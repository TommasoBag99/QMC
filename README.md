# QMC

This repository contains the Fortran 90 source code for performing Quantum
Monte Carlo (QMC) simulations to calculate the energies of chemical systems,
primarily composed of hydrogen and helium. 

The code supports both Variational Monte Carlo (VMC) and Pure Diffusion Monte
Carlo (PDMC) techniques. For comprehensive details on these algorithms and 
additional information, please refer to the `Documentation.pdf` file located in 
this folder.

## Repository Contents

- `License`: GNU General Public License.
- `Documentation.pdf`: Provides a concise overview of the Monte Carlo methods used,
                       including Variational and Pure Diffusion Monte Carlo, detailing
                       wavefunction definitions and system energy calculations.
- `QMC_main.F90`: Implements Variational and Pure Diffusion Monte Carlo methods.
- `QMC_energy.F90`: Calculates total energy as the sum of potential and kinetic
                    energies, and manages wavefunction computations and its first
                    and second derivative.
- `QMC_utilities.F90`: Provides utilities for random electron positioning,
                       average energy and acceptance rate calculation, and
                       results display.
- `QMC_input`: Input file for configuring the computations.
- `file.xyz`: Coordinates file for the atoms in the simulated chemical system.
- `Makefile`: For compiling the source code and generating the executable.


## Prerequisites

The software compilation requires the `gfortran` compiler. Install it on
Debian-based systems with the following command:
```bash
sudo apt-get install gfortran
```

## Compilation and Execution
Before you can compile and execute the Quantum Monte Carlo simulations, ensure that you 
have downloaded the project folder to your local machine. You can download the entire 
repository as a ZIP file by clicking 
[here](https://github.com/TommasoBag99/QMC/archive/refs/heads/main.zip).

Once the folder is downloaded, follow these steps to compile and run the simulations:

1. **Navigate to the Project Directory:**
   
   Open your terminal and change to the directory containing the project files.
   ```bash
   cd path/to/your/project/QMC
   ```   

2. **Compile the Code:**
   
   Execute the `make` command to compile the source files using the provided `Makefile`.
   This will automatically compile the files `QMC_main.F90`, `QMC_energy.F90`, and
   `QMC_utilities.F90` with gfortran to generate the executable `QMC_run`.
   ```bash
   make
   ```
      
3. **Run the Simulation:**

   Ensure that the input files `QMC_input` and `file.xyz` are correctly configured with all
   necessary parameters before starting the simulation.
   
   Start the simulation by running the generated executable:
   ```bash
   ./QMC_run
   ```
