This package contains routines that are being developed for full-spectrum (ion line & plasma line) ISR analysis. Most things do not work yet...


The priori models apriorimodel_plasmaline and apriorimodel_bafim_flipchem_plasmaline must be used together with modified versions of guisdap mrqmn.m and dirthe.m. Perhaps the simplest way is to include an addpath command that points to the location of these files in the "special" box  of the guisdap GUI. To verify that these versions are actually used, check that GUISDAP prints occassionally "MATLAB mrqmn forced and dirthe modified for plasma line input -- IV 2026" on the screen.

The modified versions should not be used with anything else than the plasma line prior models!


apriorimodel_bafim_flipchem_plasmaline can use both plasma line data and the Flipchem chemistry model. As this combination is probably not needed in many cases, the flipchem model can be controlled via the global variable bafim_flipchem_model_error. Negative values mean that the model is not used. For example
global bafim_flipchem_model_error
bafim_flipchem_model_error=-1
(for example in the "special" box) will disable the Flipchem part.

Also some experimental routines for plasma line detection and Gaussian spectrum fits are included. See comments in the files...


IV 2026
