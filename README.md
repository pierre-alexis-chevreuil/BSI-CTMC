

This set of codes computes individual electron trajectories in the Barrier-Suppression Ionisation (BSI) regime.
It is entirely based on Matlab codes. However, because of the large number of electrons needed to obtain statistical convergence of the results, the codes using the most CPU time have been optimised to be converted into binaries (MEX-files), see https://ch.mathworks.com/help/matlab/call-mex-file-functions.html.
The "Functions" folder contains all the necessary functions to execute the main code. All functions having "%#codgen" are meant to be converted to MEX-files. Switching between non-compiled and MEX files can be done by changing the variable "params.use_mex" from 0 to 1.
The "Examples" folder is meant to show how the different sub-functions (electric field propagation in waveplates, ionisationa and electron trajectories propagation) can be used. It is also to be meant for the creation of the individual MEX files.

Shield: [![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[![CC BY-SA 4.0][cc-by-sa-image]][cc-by-sa]

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg
