/***************************************************************************

  MuFullAnisoPowderTF.cpp

  Author: Takashi U. Ito
  e-mail: tuito@post.j-parc.jp

***************************************************************************/

/***************************************************************************
 *   Copyright (C) 2026 by Takashi U. Ito                             *
 *   tuito@post.j-parc.jp                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "MuFullAnisoPowderTF.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <mutex>
#include <algorithm>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

ClassImp(MuFullAnisoPowderTF)

namespace {

  constexpr double PI = 3.141592653589793238462643383279502884;

  // Frequencies in MHz/T.
  constexpr double GAMU_MHz_T = 135.534;
  constexpr double GAE_MHz_T  = 27992.0;

  // Nested powder average.
  // Change these to adjust accuracy and speed.
  constexpr int NDIR = 300;
  constexpr int NPSI = 20;

  // Gaussian suppression of unresolved high-frequency oscillations.
  // Unit: microsecond.
  //    constexpr double SIGMA_T_US = 0.001; // no longer used as a constant

  // Ignore extremely small spectral components.
  constexpr double AMP_EPS = 1.0e-10;

  // Conservative thread-safety for static cache.
  std::mutex g_cache_mutex;
}

// Static cache definitions.
bool MuFullAnisoPowderTF::cache_valid = false;
double MuFullAnisoPowderTF::cache_Axx = 0.0;
double MuFullAnisoPowderTF::cache_Ayy = 0.0;
double MuFullAnisoPowderTF::cache_Azz = 0.0;
double MuFullAnisoPowderTF::cache_Btf = 0.0;
double MuFullAnisoPowderTF::cache_fcut = 0.0;

std::vector<double> MuFullAnisoPowderTF::cache_freq_MHz;
std::vector<double> MuFullAnisoPowderTF::cache_amp_re;
std::vector<double> MuFullAnisoPowderTF::cache_amp_im;

MuFullAnisoPowderTF::MuFullAnisoPowderTF()
{
}

MuFullAnisoPowderTF::~MuFullAnisoPowderTF()
{
  // Do not clear static cache here.
  // The cache is shared by all instances during the musrfit process lifetime.
}

Double_t MuFullAnisoPowderTF::operator()(Double_t t,
                                        const std::vector<Double_t>& par) const
{
  if (par.size() < 4) {
    std::cerr << "MuFullAnisoPowderTF: need 6 parameters: "
              << "Axx[MHz], Ayy[MHz], Azz[MHz], f_cut[MHz], t_offset [microsec], Btf[T]."
              << std::endl;
    return 0.0;
  }

  const double Axx = par[0];  // MHz
  const double Ayy = par[1];  // MHz
  const double Azz = par[2];  // MHz
  const double fcut = par[3];  // MHz
  const double t_offset = par[4];  // microsec
  const double Btf = par[5];  // T

  // positive t_offset shifts the calculated curve to later time.
  const double t_eff = t - t_offset;

  // Safe but conservative:
  // lock during cache check/rebuild/evaluation.
  std::lock_guard<std::mutex> lock(g_cache_mutex);

  if (!cache_valid || !sameParameters(Axx, Ayy, Azz, Btf, fcut)) {
    rebuildCache(Axx, Ayy, Azz, Btf, fcut);
  }

  return evalFromCache(t_eff);
}

bool MuFullAnisoPowderTF::sameParameters(double Axx,
					 double Ayy,
					 double Azz,
					 double Btf,
					 double fcut)
{
  // Exact comparison is appropriate for cache reuse within one FCN evaluation.
  // If desired, replace by tolerance comparison.
  return (Axx == cache_Axx &&
          Ayy == cache_Ayy &&
          Azz == cache_Azz &&
          Btf == cache_Btf &&
	  fcut == cache_fcut);
}

double MuFullAnisoPowderTF::evalFromCache(double t)
{
  double P = 0.0;

  const std::size_t n = cache_freq_MHz.size();

  for (std::size_t k = 0; k < n; ++k) {
    const double theta = 2.0 * PI * cache_freq_MHz[k] * t;
    const double cs = std::cos(theta);
    //    const double sn = std::sin(theta);
    //    //     real[ amp * exp(-i theta) ]
    //    //     = Re(amp) cos(theta) + Im(amp) sin(theta)
    //    if(std::abs(cache_amp_im[k]) > AMP_EPS){
    //      P += cache_amp_re[k] * cs + cache_amp_im[k] * sn;      
    //      std::cout << "c_im = " << cache_amp_im[k] << std::endl;
    //    }
    //    else{
    //      P += cache_amp_re[k] * cs;
    //    }
    P += cache_amp_re[k] * cs;  //  cache_amp_im[k]  is always zero
    
  }
  
  return P;
}

void MuFullAnisoPowderTF::rebuildCache(double Axx,
				       double Ayy,
				       double Azz,
				       double Btf,
				       double fcut)
{
  cache_freq_MHz.clear();
  cache_amp_re.clear();
  cache_amp_im.clear();

  cache_freq_MHz.reserve(static_cast<std::size_t>(NDIR) * NPSI * 16);
  cache_amp_re.reserve(static_cast<std::size_t>(NDIR) * NPSI * 16);
  cache_amp_im.reserve(static_cast<std::size_t>(NDIR) * NPSI * 16);

  std::cout << "MuFullAnisoPowderTF: rebuilding cache for "
            << "Axx=" << Axx << " MHz, "
            << "Ayy=" << Ayy << " MHz, "
            << "Azz=" << Azz << " MHz, "
            << "Btf=" << Btf << " T"
	    << "fcut=" << fcut << " MHz"
            << std::endl;

  Mat4 Sx, Sy, Sz, Ix, Iy, Iz;
  buildSpinOperators(Sx, Sy, Sz, Ix, Iy, Iz);

  Mat4 Iop[3] = {Ix, Iy, Iz};
  Mat4 Sop[3] = {Sx, Sy, Sz};

  // Precompute I_i S_j.
  Mat4 IS[3][3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      IS[i][j] = multiply(Iop[i], Sop[j]);
    }
  }

  // TF geometry:
  // B || lab z
  // P_mu(0), detection || lab x
  const Mat4 Iinit = Ix;
  const Mat4 Idet  = Ix;

  // rho0 = 1/4 * identity + 1/2 * Iinit.
  Mat4 rho0 = add(scale(eye4(), 0.25), scale(Iinit, 0.5));

  const double wm = GAMU_MHz_T * Btf;
  const double we = GAE_MHz_T  * Btf;

  const int Ntot = NDIR * NPSI;

  for (int idir = 1; idir <= NDIR; ++idir) {

    double Q0[3][3];
    makeOrientationFrame(idir, NDIR, Q0);

    for (int ipsi = 1; ipsi <= NPSI; ++ipsi) {

      const double psi =
        2.0 * PI * static_cast<double>(ipsi - 1)
        / static_cast<double>(NPSI);

      double Rz[3][3];
      double Q[3][3];

      makeRz(psi, Rz);
      mat3mul(Q0, Rz, Q);

      double Alab[3][3];
      computeAlab(Axx, Ayy, Azz, Q, Alab);

      Mat4 H = zero4();

      // Zeeman terms, B || z:
      // H = -gamma_mu B Iz + gamma_e B Sz + I^T A S.
      H = add(H, scale(Iz, -wm));
      H = add(H, scale(Sz,  we));

      // Hyperfine term.
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          H = add(H, scale(IS[i][j], Alab[i][j]));
        }
      }

      double E[4];
      cd Vcol[4][4];

      diagonalizeHermitian4(H, E, Vcol);

      // Convert eigenvectors to Mat4 V with columns as eigenvectors.
      Mat4 V = zero4();
      for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
          V.m[r][c] = Vcol[r][c];
        }
      }

      Mat4 rho_e = transformDaggerLeft(V, rho0);
      Mat4 O_e   = transformDaggerLeft(V, Idet);

      // Spectral decomposition:
      //
      // P(t) = sum_ab 2 O_e(b,a) rho_e(a,b)
      //        exp[-i 2 pi (E_a - E_b) t]
      //
      // E in MHz, t in microsecond.
      for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {

          const double f = E[a] - E[b]; // MHz

          cd amp =
            2.0 * O_e.m[b][a] * rho_e.m[a][b]
            / static_cast<double>(Ntot);

          //// Gaussian time-resolution / high-frequency suppression.
	  //  const double x = 2.0 * PI * f * sigma_t;
          //  amp *= std::exp(-0.5 * x * x);

          //// Low-pass filter
	  if(std::abs(f) > fcut){
	    continue;
	  }

          if (std::abs(amp) > AMP_EPS) {
            cache_freq_MHz.push_back(f);
            cache_amp_re.push_back(amp.real());
            cache_amp_im.push_back(amp.imag());
          }
        }
      }
    }
  }

  cache_Axx = Axx;
  cache_Ayy = Ayy;
  cache_Azz = Azz;
  cache_Btf = Btf;
  cache_fcut = fcut;
  cache_valid = true;

  std::cout << "MuFullAnisoPowderTF: cache rebuilt. "
            << "Number of spectral components = "
            << cache_freq_MHz.size()
            << std::endl;
}

// ----------------------------------------------------------------------
// 4x4 matrix helpers
// ----------------------------------------------------------------------

MuFullAnisoPowderTF::Mat4 MuFullAnisoPowderTF::zero4()
{
  Mat4 A;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      A.m[i][j] = cd(0.0, 0.0);
    }
  }
  return A;
}

MuFullAnisoPowderTF::Mat4 MuFullAnisoPowderTF::eye4()
{
  Mat4 A = zero4();
  for (int i = 0; i < 4; ++i) {
    A.m[i][i] = cd(1.0, 0.0);
  }
  return A;
}

MuFullAnisoPowderTF::Mat4
MuFullAnisoPowderTF::add(const Mat4& A, const Mat4& B)
{
  Mat4 C;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      C.m[i][j] = A.m[i][j] + B.m[i][j];
    }
  }
  return C;
}

MuFullAnisoPowderTF::Mat4
MuFullAnisoPowderTF::scale(const Mat4& A, cd c)
{
  Mat4 B;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      B.m[i][j] = c * A.m[i][j];
    }
  }
  return B;
}

MuFullAnisoPowderTF::Mat4
MuFullAnisoPowderTF::multiply(const Mat4& A, const Mat4& B)
{
  Mat4 C = zero4();

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      cd s(0.0, 0.0);
      for (int k = 0; k < 4; ++k) {
        s += A.m[i][k] * B.m[k][j];
      }
      C.m[i][j] = s;
    }
  }

  return C;
}

MuFullAnisoPowderTF::Mat4
MuFullAnisoPowderTF::dagger(const Mat4& A)
{
  Mat4 B;

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      B.m[i][j] = std::conj(A.m[j][i]);
    }
  }

  return B;
}

MuFullAnisoPowderTF::Mat4
MuFullAnisoPowderTF::transformDaggerLeft(const Mat4& V, const Mat4& A)
{
  // Return V^\dagger A V.
  Mat4 VH = dagger(V);
  return multiply(multiply(VH, A), V);
}

// ----------------------------------------------------------------------
// Spin operators
// ----------------------------------------------------------------------

MuFullAnisoPowderTF::Mat4
MuFullAnisoPowderTF::kron2(const cd A[2][2], const cd B[2][2])
{
  Mat4 K = zero4();

  for (int ae = 0; ae < 2; ++ae) {
    for (int be = 0; be < 2; ++be) {
      for (int am = 0; am < 2; ++am) {
        for (int bm = 0; bm < 2; ++bm) {
          const int row = 2 * ae + am;
          const int col = 2 * be + bm;
          K.m[row][col] = A[ae][be] * B[am][bm];
        }
      }
    }
  }

  return K;
}

void MuFullAnisoPowderTF::buildSpinOperators(Mat4& Sx, Mat4& Sy, Mat4& Sz,
                                            Mat4& Ix, Mat4& Iy, Mat4& Iz)
{
  const cd I2[2][2] = {
    {cd(1.0,0.0), cd(0.0,0.0)},
    {cd(0.0,0.0), cd(1.0,0.0)}
  };

  const cd sx[2][2] = {
    {cd(0.0,0.0), cd(1.0,0.0)},
    {cd(1.0,0.0), cd(0.0,0.0)}
  };

  const cd sy[2][2] = {
    {cd(0.0,0.0), cd(0.0,-1.0)},
    {cd(0.0,1.0), cd(0.0,0.0)}
  };

  const cd sz[2][2] = {
    {cd(1.0,0.0), cd(0.0,0.0)},
    {cd(0.0,0.0), cd(-1.0,0.0)}
  };

  // Tensor product order: electron spin outer, muon spin inner.
  Sx = scale(kron2(sx, I2), 0.5);
  Sy = scale(kron2(sy, I2), 0.5);
  Sz = scale(kron2(sz, I2), 0.5);

  Ix = scale(kron2(I2, sx), 0.5);
  Iy = scale(kron2(I2, sy), 0.5);
  Iz = scale(kron2(I2, sz), 0.5);
}

// ----------------------------------------------------------------------
// GSL diagonalization
// ----------------------------------------------------------------------

void MuFullAnisoPowderTF::diagonalizeHermitian4(const Mat4& H,
                                               double E[4],
                                               cd V[4][4])
{
  gsl_matrix_complex* Hgsl = gsl_matrix_complex_alloc(4, 4);
  gsl_vector* eval = gsl_vector_alloc(4);
  gsl_matrix_complex* evec = gsl_matrix_complex_alloc(4, 4);
  gsl_eigen_hermv_workspace* w = gsl_eigen_hermv_alloc(4);

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      gsl_complex z;
      GSL_SET_COMPLEX(&z, H.m[i][j].real(), H.m[i][j].imag());
      gsl_matrix_complex_set(Hgsl, i, j, z);
    }
  }

  gsl_eigen_hermv(Hgsl, eval, evec, w);
  gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

  for (int a = 0; a < 4; ++a) {
    E[a] = gsl_vector_get(eval, a);

    for (int r = 0; r < 4; ++r) {
      gsl_complex z = gsl_matrix_complex_get(evec, r, a);
      V[r][a] = cd(GSL_REAL(z), GSL_IMAG(z));
    }
  }

  gsl_eigen_hermv_free(w);
  gsl_matrix_complex_free(evec);
  gsl_vector_free(eval);
  gsl_matrix_complex_free(Hgsl);
}

// ----------------------------------------------------------------------
// 3D geometry helpers
// ----------------------------------------------------------------------

double MuFullAnisoPowderTF::dot3(const double a[3], const double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void MuFullAnisoPowderTF::cross3(const double a[3],
                                const double b[3],
                                double c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

void MuFullAnisoPowderTF::normalize3(double a[3])
{
  const double n = std::sqrt(dot3(a,a));

  if (n <= 0.0) {
    throw std::runtime_error("normalize3: zero vector");
  }

  a[0] /= n;
  a[1] /= n;
  a[2] /= n;
}

void MuFullAnisoPowderTF::mat3mul(const double A[3][3],
                                 const double B[3][3],
                                 double C[3][3])
{
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      double s = 0.0;
      for (int k = 0; k < 3; ++k) {
        s += A[i][k] * B[k][j];
      }
      C[i][j] = s;
    }
  }
}

void MuFullAnisoPowderTF::mat3transpose(const double A[3][3],
                                       double AT[3][3])
{
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      AT[i][j] = A[j][i];
    }
  }
}

void MuFullAnisoPowderTF::computeEulerR(double R[3][3])
{
  // Default Euler angles are zero.
  // Kept explicitly for later extension if needed.
  const double alpha = 0.0;
  const double beta  = 0.0;
  const double gamma = 0.0;

  R[0][0] =  std::cos(alpha)*std::cos(beta)*std::cos(gamma)
           - std::sin(alpha)*std::sin(gamma);

  R[0][1] = -std::cos(alpha)*std::cos(beta)*std::sin(gamma)
           - std::sin(alpha)*std::cos(gamma);

  R[0][2] =  std::cos(alpha)*std::sin(beta);

  R[1][0] =  std::sin(alpha)*std::cos(beta)*std::cos(gamma)
           + std::cos(alpha)*std::sin(gamma);

  R[1][1] = -std::sin(alpha)*std::cos(beta)*std::sin(gamma)
           + std::cos(alpha)*std::cos(gamma);

  R[1][2] =  std::sin(alpha)*std::sin(beta);

  R[2][0] = -std::sin(beta)*std::cos(gamma);
  R[2][1] =  std::sin(beta)*std::sin(gamma);
  R[2][2] =  std::cos(beta);
}

void MuFullAnisoPowderTF::computeAlab(double Axx,
                                     double Ayy,
                                     double Azz,
                                     const double Q[3][3],
                                     double Alab[3][3])
{
  // A0 = R * diag(Axx,Ayy,Azz) * R^T.
  double R[3][3];
  computeEulerR(R);

  double Aprin[3][3] = {
    {Axx, 0.0, 0.0},
    {0.0, Ayy, 0.0},
    {0.0, 0.0, Azz}
  };

  double Rt[3][3];
  double tmp[3][3];
  double A0[3][3];

  mat3transpose(R, Rt);
  mat3mul(R, Aprin, tmp);
  mat3mul(tmp, Rt, A0);

  // Alab = Q * A0 * Q^T.
  double Qt[3][3];
  mat3transpose(Q, Qt);
  mat3mul(Q, A0, tmp);
  mat3mul(tmp, Qt, Alab);
}

void MuFullAnisoPowderTF::makeOrientationFrame(int idir,
                                              int Ndir,
                                              double Q0[3][3])
{
  const double golden = (1.0 + std::sqrt(5.0)) / 2.0;

  const double zc =
    1.0 - 2.0 * (static_cast<double>(idir) - 0.5)
    / static_cast<double>(Ndir);

  const double rc = std::sqrt(std::max(0.0, 1.0 - zc*zc));

  const double phi =
    2.0 * PI * static_cast<double>(idir - 1) / golden;

  double zaxis[3] = {
    rc * std::cos(phi),
    rc * std::sin(phi),
    zc
  };

  normalize3(zaxis);

  double trial[3];

  if (std::fabs(zaxis[2]) < 0.9) {
    trial[0] = 0.0;
    trial[1] = 0.0;
    trial[2] = 1.0;
  } else {
    trial[0] = 1.0;
    trial[1] = 0.0;
    trial[2] = 0.0;
  }

  const double proj = dot3(zaxis, trial);

  double xaxis[3] = {
    trial[0] - proj*zaxis[0],
    trial[1] - proj*zaxis[1],
    trial[2] - proj*zaxis[2]
  };

  normalize3(xaxis);

  double yaxis[3];
  cross3(zaxis, xaxis, yaxis);
  normalize3(yaxis);

  // Columns of Q0 are xaxis, yaxis, zaxis.
  for (int i = 0; i < 3; ++i) {
    Q0[i][0] = xaxis[i];
    Q0[i][1] = yaxis[i];
    Q0[i][2] = zaxis[i];
  }
}

void MuFullAnisoPowderTF::makeRz(double psi, double Rz[3][3])
{
  Rz[0][0] =  std::cos(psi);
  Rz[0][1] = -std::sin(psi);
  Rz[0][2] =  0.0;

  Rz[1][0] =  std::sin(psi);
  Rz[1][1] =  std::cos(psi);
  Rz[1][2] =  0.0;

  Rz[2][0] =  0.0;
  Rz[2][1] =  0.0;
  Rz[2][2] =  1.0;
}
