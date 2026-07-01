/***************************************************************************

  MuFullAnisoPowderTF.h

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

#ifndef MU_FULL_ANISO_POWDER_TF_H
#define MU_FULL_ANISO_POWDER_TF_H

#include "PUserFcnBase.h"

#include <vector>
#include <complex>

class MuFullAnisoPowderTF : public PUserFcnBase {

public:
  MuFullAnisoPowderTF();
  virtual ~MuFullAnisoPowderTF();

  Bool_t NeedGlobalPart() const { return false; }
  void SetGlobalPart(std::vector<void *> &, UInt_t) { }
  Bool_t GlobalPartIsValid() const { return true; }

  // musrfit user function:
  // par[0] = Axx [MHz]
  // par[1] = Ayy [MHz]
  // par[2] = Azz [MHz]
  // par[3] = Btf [T]
  Double_t operator()(Double_t t, const std::vector<Double_t>& par) const;

  ClassDef(MuFullAnisoPowderTF, 1)

private:
  using cd = std::complex<double>;

  struct Mat4 {
    cd m[4][4];
  };

  // Static cache.
  static bool cache_valid;
  static double cache_Axx;
  static double cache_Ayy;
  static double cache_Azz;
  static double cache_Btf;
  static double cache_fcut;  

  static std::vector<double> cache_freq_MHz;
  static std::vector<double> cache_amp_re;
  static std::vector<double> cache_amp_im;

  static bool sameParameters(double Axx, double Ayy, double Azz, double Btf, double fcut);
  static void rebuildCache(double Axx, double Ayy, double Azz, double Btf, double fcut);
  static double evalFromCache(double t);

  // 4x4 matrix helpers.
  static Mat4 zero4();
  static Mat4 eye4();
  static Mat4 add(const Mat4& A, const Mat4& B);
  static Mat4 scale(const Mat4& A, cd c);
  static Mat4 multiply(const Mat4& A, const Mat4& B);
  static Mat4 dagger(const Mat4& A);
  static Mat4 transformDaggerLeft(const Mat4& V, const Mat4& A);

  static Mat4 kron2(const cd A[2][2], const cd B[2][2]);

  static void buildSpinOperators(Mat4& Sx, Mat4& Sy, Mat4& Sz,
                                 Mat4& Ix, Mat4& Iy, Mat4& Iz);

  static void diagonalizeHermitian4(const Mat4& H,
                                    double E[4],
                                    cd V[4][4]);

  // 3D geometry helpers.
  static double dot3(const double a[3], const double b[3]);
  static void cross3(const double a[3], const double b[3], double c[3]);
  static void normalize3(double a[3]);

  static void mat3mul(const double A[3][3],
                      const double B[3][3],
                      double C[3][3]);

  static void mat3transpose(const double A[3][3],
                            double AT[3][3]);

  static void computeEulerR(double R[3][3]);

  static void computeAlab(double Axx, double Ayy, double Azz,
                          const double Q[3][3],
                          double Alab[3][3]);

  static void makeOrientationFrame(int idir, int Ndir,
                                   double Q0[3][3]);

  static void makeRz(double psi, double Rz[3][3]);
};

#endif

