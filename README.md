# MuFullAnisoPowderTF

`MuFullAnisoPowderTF` is a plug-in user function for [musrfit](https://rmlmcfadden.github.io/musr/musrfit/) that calculates the transverse-field μSR polarization of powder-averaged muonium with a fully anisotropic hyperfine tensor characterized by three principal values [A_{xx}, A_{yy}, A_{zz}].
The function is intended for modelling non-axial Mu\(^0\) states in weak transverse-field μSR spectra.

## Author

T. U. Ito  
Japan Atomic Energy Agency

## Disclaimer

This software is provided “as is” for research purposes, with absolutely no warranty.  
Users are responsible for checking the Hamiltonian, parameter conventions, numerical settings, and physical interpretation of the calculated spectra.

## Model

The spin Hamiltonian is

\[
H =
-\gamma_\mu B I_z
+
\gamma_e B S_z
+
\sum_{i,j} A_{ij} I_i S_j .
\]

The magnetic field is fixed along the laboratory \(z\)-axis,

\[
B \parallel z_{\rm lab},
\]

and the calculated geometry is transverse field,

\[
P_\mu(0) \parallel x_{\rm lab}.
\]

The function returns the powder-averaged transverse-field muonium polarization,

\[
P_{\rm TF}(t).
\]

The gyromagnetic ratios used in the source code are

    gamma_mu = 135.534 MHz/T
    gamma_e  = 27992.0 MHz/T

## Fitting Parameters

The user function takes six parameters:

    par[0] = Axx       [MHz]
    par[1] = Ayy       [MHz]
    par[2] = Azz       [MHz]
    par[3] = f_cut     [MHz]
    par[4] = t_offset  [microsec]
    par[5] = Btf      [T]

### Axx, Ayy, Azz

Principal values of the fully anisotropic hyperfine tensor in MHz.

### Btf

External transverse magnetic field in tesla.

### f_cut

Frequency cutoff in MHz. Spectral components with

\[
|f| > f_{\rm cut}
\]

are discarded. This is used to suppress unresolved high-frequency components, such as GHz-scale hyperfine oscillations, which may otherwise cause digital aliasing artifacts.

Although `f_cut` is exposed as a fitting parameter for convenience in the `musrfit` GUI, it should usually be regarded as a numerical cutoff parameter rather than a physical parameter.

### t_offset

Time offset in microseconds. The function evaluates

\[
t_{\rm eff}=t-t_{\rm offset}.
\]

A positive `t_offset` shifts the calculated curve to later times.

## Powder Averaging

A fully random powder average is calculated by an explicit nested average over orientational degrees of freedom,

\[
S^2 \times S^1.
\]

The source code currently uses

    constexpr int NDIR = 300;
    constexpr int NPSI = 20;

where `NDIR` is the number of sampled local-axis directions and `NPSI` is the number of rotations around each local axis.

The total number of powder orientations is

\[
N_{\rm total}=N_{\rm dir}N_\psi.
\]

With the default values,

    Ntotal = 300 × 20 = 6000

orientations are used.

To change the powder-average accuracy, edit `NDIR` and `NPSI` in the source code and recompile.

## Spectral Cache

For a given parameter set,

    Axx, Ayy, Azz, Btf, f_cut

the Hamiltonian is diagonalized for all powder orientations. The result is stored in static memory as spectral components,

\[
(f_k, c_k).
\]

The polarization is evaluated as

\[
P(t)=
\sum_k
{\rm Re}
[
c_k \exp(-i2\pi f_k t)
].
\]

The cache is rebuilt only when one of

    Axx
    Ayy
    Azz
    Btf
    f_cut

changes.

`t_offset` does not trigger cache rebuilding because it only shifts the time argument.

## Dependencies

This plug-in requires:

- `musrfit`
- ROOT
- GSL

## Compilation

Edit the `Makefile` if necessary, then run

    make

To install the shared library into the ROOT/musrfit library directory, run

    make install

Check the installation path in the `Makefile` before running `make install`.

## Typical Use

The typical call through the msr-file would be
```
###############################################################
FITPARAMETER
#      Nr. Name        Value     Step      Pos_Error  Boundaries
        1 alpha    1.      0.0       none        0      2      
        2 asy      0.2     0.01      none        0.1    1    
        3 Axx      4000    0.0       none        0.1    5000
        4 Ayy      4001.4  0.1       none        0.1    5000
        5 Azz      4001.5  0.1       none        0.1    5000      
        6 f_cut    1000.0  0.0       none        0      1000
        7 t_offset 0.02    0.0       none        0      1
        8 TF       0.00023 0.0       none        0      0.01

##############################################################
THEORY
asymmetry      2
userFcn  libMuFullAnisoPowderTFLibrary.so MuFullAnisoPowderTF 3 4 5 6 7 8  (Axx Ayy Azz f_cut t_offset TF)
###############################################################
```

## Notes

- Time is assumed to be in microseconds.
- Frequencies are in MHz.
- Magnetic field is in tesla.
- Only TF geometry is implemented.
- The powder average assumes a completely random powder.
- `NDIR` and `NPSI` are compile-time constants.
- `f_cut` is exposed as a fitting parameter but mainly controls numerical high-frequency filtering.
- `t_offset` shifts the time argument and does not trigger cache rebuilding.

## Citation

If you use `MuFullAnisoPowderTF` in a scientific publication, presentation, or data analysis report, please cite this GitHub repository.


## License

Please see the `LICENSE` file.
