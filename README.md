# Advanced Magnetic Field library

See [Licence](https://github.com/Alexey-Stupishin/AMaFiL/blob/master/LICENCE.md)

<!--- [doi:10.5281/zenodo.3896222](https://zenodo.org/record/3896222) --->

The repository contains:
* sources of Magnetic Field Library (C++11)
* Windows .dll
* projects for build under Windows and Linux

## Dependencies

Wrappers for the library are implemented in:
* IDL ([idlAMaFiL](https://github.com/Alexey-Stupishin/idlAMaFiL)). idlAMaFiL used as submodule of [GX-simulator](https://github.com/Gelu-Nita/GX_SIMULATOR) package, which included in [Solar Soft](http://www.lmsal.com/solarsoft/sswdoc/sswdoc_jtop.html) environment.
* Python ([pyAMaFiL](https://github.com/Alexey-Stupishin/pyAMaFiL)). Porting to the ([pyAMPP](https://github.com/suncast-org/pyAMPP)) in progress.

This repository contains submodule [CPP-Common](https://github.com/Alexey-Stupishin/CPP-Common).

Working with submodules can be found [here](https://git-scm.com/book/en/v2/Git-Tools-Submodules).

## Routines
* _/binaries/WWNLFFFReconstruction.dll_ - Calculation library for Windows
* _/sources_ - source codes (C++11)
* _/Windows_ - MSVC project for Windows (Visual Studio 2015)
* _/Linux_ - makefile (g++) for Linux

see [Changelog](https://github.com/Alexey-Stupishin/AMaFiL/blob/master/CHANGELOG.md) for history update

## References
For NLFFF weighted method please refer to:
Wiegelmann, T. Optimization code with weighting function for the reconstruction of coronal magnetic fields. _Solar Phys_., 2004, __219__, 87108. ([doi:10.1023/B:SOLA.0000021799.39465.36](https://link.springer.com/article/10.1023/B:SOLA.0000021799.39465.36), [ADS:2004SoPh..219...87W](https://ui.adsabs.harvard.edu/abs/2004SoPh..219...87W/abstract), [arXiv:0802.0124](https://arxiv.org/abs/0802.0124))

Some proves and using of this library can be found in:

Fleishman, G., Anfinogentov, S., Loukitcheva, M., Mysh'yakov, I., Stupishin, A. Casting the Coronal Magnetic Field Reconstruction Tools in 3D Using the MHD Bifrost Model. _ApJ_, 2017, __839__, 30 ([doi:10.3847/1538-4357/aa6840](https://iopscience.iop.org/article/10.3847/1538-4357/aa6840), [ADS:2017ApJ...839...30F](https://ui.adsabs.harvard.edu/abs/2017ApJ...839...30F/abstract), [arXiv:1703.06360](https://arxiv.org/abs/1703.06360))

Stupishin, A., Kaltman, T., Bogod, V., Yasnov, L. Modeling of Solar Atmosphere Parameters Above Sunspots Using RATAN-600 Microwave Observations. _Solar Phys_, 2018, __293__, 13 ([doi:10.1007/s11207-017-1228-7](https://link.springer.com/article/10.1007/s11207-017-1228-7), [ADS:2018SoPh..293...13S](https://ui.adsabs.harvard.edu/abs/2018SoPh..293...13S/abstract))

Anfinogentov, S., Stupishin, A., Myshyakov, I., Fleishman, G. Record-breaking Coronal Magnetic Field in Solar Active Region 12673. _ApJL_, 2019, __880__, L29 ([doi:10.3847/2041-8213/ab3042](https://iopscience.iop.org/article/10.3847/2041-8213/ab3042), [ADS:2019ApJ...880L..29A](https://ui.adsabs.harvard.edu/abs/2019ApJ...880L..29A/abstract), [arXiv:1907.06398](https://arxiv.org/abs/1907.06398))

Fleishman, G., Myshyakov, I., Stupishin, A., Loukitcheva, M., Anfinogentov, S. Force-free Field Reconstructions Enhanced by Chromospheric Magnetic Field Data. _ApJ_, 2019, __870__, 101 ([doi:10.3847/1538-4357/aaf384](https://iopscience.iop.org/article/10.3847/1538-4357/aaf384), [ADS:2019ApJ...870..101F](https://ui.adsabs.harvard.edu/abs/2019ApJ...870..101F/abstract), [arXiv:1811.02093](https://arxiv.org/abs/1811.02093))

Fleishman, G., Anfinogentov, Stupishin, A., Kuznetsov, A., Nita, G. Coronal Heating Law Constrained by Microwave Gyroresonant Emission. _ApJ_, 2021, __909__, 89 ([doi:10.3847/1538-4357/abdab](https://iopscience.iop.org/article/10.3847/1538-4357/abdab1), [ADS:2021ApJ...909...89F](https://ui.adsabs.harvard.edu/abs/2021ApJ...909...89F/abstract), [arXiv:2101.03651](https://arxiv.org/abs/2101.03651))

Gelu M. Nita, Gregory D. Fleishman, Alexey A. Kuznetsov, Sergey A. Anfinogentov, Alexey G. Stupishin, Eduard P. Kontar, Samuel J. Schonfeld, James A. Klimchuk, and Dale E. Gary. Data-constrained Solar Modeling with GX Simulator. _ApJSS_, 2023, __267__, 6 ([doi:10.3847/1538-4365/acd343](https://link.springer.com/article/10.1023/B:SOLA.0000021799.39465.36), [ADS:2023ApJS..267....6N](https://ui.adsabs.harvard.edu/abs/2023ApJS..267....6N/abstract), [arXiv:2301.00795](https://arxiv.org/abs/2301.00795))
