# HPB_tools

Various bits of code written to carry out numerical calculations of atomic energy levels, state admixtures, and more for alkali-metal species in the hyperfine Paschen-Back regime.

Currently tested for rubidium (n = 5) in experiments realised in QLM @ Durham University.

`atomdata`: Simple database of atomic structure constants for use in calculations.

`uncoupledbasis`: Code for numerically calculating the atomic energy levels, as well as the state admixtures, at a given magnetic field strength. Also used for generating energy-level diagrams as a function of magnetic field.

`durhamcolours`: Custom colour codes for use in plotting.

`BIG_diagram`: Code for plotting a "big" diagram composed of: 1) atomic energy levels as a function of magnetic field strength, 2) calculated atomic absorption spectrum at maximum value of magnetic field and 3) dipole-allowed transitions from ground-state to excited-state energy levels ($$\pi,\sigma^{\pm}$$).
*Note: requires installation of [ElecSus](https://github.com/jameskeaveney/ElecSus) for calculation of absorption spectra.*

*Additional details coming soon...*
