# Falaise_Modules

Reconstruction and Analysis modules for the SuperNEMO Experiment


## FLMANU modules library

#### FLMANU_CD

Basic inspection of "CD" bank (calorimeter) with ROOT output file.

#### FLMANU_GENBB

Basic inspection of GENBB generator and particle spectra from "SD" bank

#### FLMANU_NeMg

Basic inspection of "PTD" banks splitted in topologies of N electron(s) with N gamma(s) with ROOT output file.

#### FLMANU_Test

Basic module doing a removal of zero calorimeter hit.


## ZeroHitCleaner library

Short module removing events based on (customizable) number of hits in tracker and calorimeter according to SD bank. Usefull to be run after flsimulate in order to reduce output file and speed up processing of reconstruction and analysis.
