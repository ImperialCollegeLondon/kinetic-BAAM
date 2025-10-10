

    +=============================================================================+
    |                                                                             |
    |                                                                             |
    |                                                                             |
    |                                                                             |
    |   `7MM            `7MM"""Yp,       db            db      `7MMM.     ,MMF'   |
    |     MM              MM    Yb      ;MM:          ;MM:       MMMb    dPMM     |
    |     MM  ,MP'        MM    dP     ,V^MM.        ,V^MM.      M YM   ,M MM     |
    |     MM ;Y           MM"""bg.    ,M  `MM       ,M  `MM      M  Mb  M' MM     |
    |     MM;Mm   mmmmm   MM    `Y    AbmmmqMA      AbmmmqMA     M  YM.P'  MM     |
    |     MM `Mb.         MM    ,9   A'     VML    A'     VML    M  `YM'   MM     |
    |   .JMML. YA.      .JMMmmmd9  .AMA.   .AMMA..AMA.   .AMMA..JML. `'  .JMML.   |
    |                                                                             |
    |                                                                             |
    |                                                                             |
    |                                                                             |
    +=============================================================================+
				
						 kinetic-Batch Adsorber Analogue Model     
						    Copyright (C) 2025 Hassan Azzan
						 
						Imperial College London, United Kingdom
						     Multiphase Systems Laboratory
						
					 Authors:  Hassan Azzan (HA), Ayca Yilmaz (AY)
		

## INTRODUCTION
The repository contains the MATLAB (and equivalent python) scripts and functions of a model for a zero-dimensional time variant kinetics-based cyclic adsorption process. Contains routines for forward simulations and multi-objective optimization of PVSA, and VSA processes.

![alt text](https://github.com/ImperialCollegeLondon/kinetic-BAAM/blob/main/rawData/optExample.png?raw=true)

## INSTALLATION

### Dependencies

The following dependencies are required for the proper execution of this program using MATLAB or python.

1. MATLAB 2020b or later 
    - Global Optimization Toolbox *[required for optimization]*
2. python 3.12 or later 
    - numpy>=1.26
    - scipy>=1.13
    - pymoo>=0.6.0 *[required for optimization]*
    - matplotlib>=3.8
    - matplotlib-inline>=0.1.6
    - pandas>=2.2

### Installation instructions

1. Clone the full software package from the GitHub server into the preferred installation directory (e.g. Desktop, My Documents). The command is as follows:
```sh
git clone https://github.com/ImperialCollegeLondon/kinetic-BAAM.git
```

## REFERENCES
 
## AUTHORS

### Maintainers of the repository
* Dr. Hassan Azzan (hassan.azzan15@imperial.ac.uk)

### Project contributors
* Ms. Ayca Yilmaz
* Prof. Ronny Pini (r.pini@imperial.ac.uk)


## Change Log
All notable changes to this project will be documented in this file.
 
###  Create python version of the platform - 2025-10-10
#### Added
- kBAAM_ODEs_nonIsothermal_ND.py
- kBAAM_Outputs_nonIsothermal.py
- LDFCoefficient.py
- DSL.py
- run_NSGA.py
- kBAAM_run_NSGA.py
- pythonPort/


###  Initial creation of the repository - 2025-10-07
#### Added
- kBAAM_ODEs_nonIsothermal_ND.m
- kBAAM_Outputs_nonIsothermal.m
- LDFCoefficient.m
- DSL.m
- run_NSGA.m
- kBAAM_run_NSGA.m
- plotFunctions/
- rawData/
- matFiles/
