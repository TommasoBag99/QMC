# QMC

This repository contains the Fortran 90 source code for performing Quantum
Monte Carlo (QMC) simulations to calculate the energies of chemical systems.

## Repository Contents

- `License`: GNU General Public License.
- `Makefile`: For compiling the source code and generating the executable.
- `QMC_main.F90`: Implements Variational and Pure Diffusion Monte Carlo methods.
- `QMC_energy.F90`: Calculates total energy as the sum of potential and kinetic
                    energies, and manages wavefunction computations and its first
                    and second derivative.
- `QMC_utilities.F90`: Provides utilities for random electron positioning,
                       average energy and acceptance rate calculation, and
                       results display.
- `QMC_input`: Input file for configuring the computations.
- `file.xyz`: Coordinates file for the atoms in the simulated chemical system.

## Prerequisites

The software compilation requires the `gfortran` compiler. Install it on
Debian-based systems with the following command:
```bash
sudo apt-get install gfortran

## Compilation and Execution

To compile and execute the Quantum Monte Carlo simulations, follow these simple steps:

1. **Navigate to the Project Directory:**
   Open your terminal and change to the directory containing the project files.
   ```bash
   cd path/to/your/project
