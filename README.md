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

## Input Files Description

The simulation uses two main input files to configure and execute the Quantum Monte Carlo 
calculations:

### `QMC_input`

This file contains the parameters needed to set up and control the simulation method and 
environment. It specifies the Monte Carlo method (variational or diffusion), the system 
geometry (`file.xyz`), the number of electrons, and various simulation parameters such as
the Jastrow factor, time step, number of steps, number of walkers, reference energy, and 
projection time (the last two are specific to the Pure Diffusion method).

**Example of `QMC_input`:**
```bash
dif       ! Method var/dif     
h3.xyz    ! Geometry file        
2         ! Number of electrons  
1.2       ! Jastrow factor (a)

0.05      ! Time step            
100000    ! Number of steps
50        ! Number of walkers

-1.342    ! Reference energy !! ONLY FOR PURE DIFFUSION
100       ! Projection time  !! ONLY FOR PURE DIFFUSION
```

### `file.xyz`

This file defines the positions of the atoms in the chemical system being simulated, with
each coordinate given in Angstroms. Each line after the first specifies an atom type followed
by its x, y, and z coordinates in space.

**Example of `file.xyz`:**
```bash
3 # number of atoms

H -0.049222 0.000000 -0.085255 
H -0.049222 0.000000  0.785255 
H  0.704662 0.000000  0.350000
```

## Output Description

The output of the simulation is displayed directly in the command prompt as a simple table. 
Below is a breakdown of the information you will receive:

### System Coordinates

This table displays the coordinates of the atoms in the system under simulation, listing each 
atom type along with its respective x, y, and z coordinates.
```bash
+--------+-----------------------------------+
| SYSTEM |      X          Y          Z      |
+--------+-----------------------------------+
|   H    | -0.0930161  0.0000000 -0.1611086  |
|   H    | -0.0930161  0.0000000  1.4839168  |
|   H    |  1.3316181  0.0000000  0.6614041  |
+--------+-----------------------------------+
```

### Charge

Displays the total charge of the system.
```bash
+--------------------------------------------+
| CHARGE =  1                                |
+--------------------------------------------+
```

### Quantum Monte Carlo Results

This section summarizes the results of the Quantum Monte Carlo simulation, including the reference 
energy, the calculated Quantum Monte Carlo energy with its standard error, and the acceptance rate 
with its standard error.
```bash
+--------------------------------------------+
|     PURE DIFFUSION QUANTUM MONTE CARLO     |
+--------------------------------------------+
| ENERGY REF = -1.342000                     |
| ENERGY QMC = -1.342193 +/- 0.002652        |
| ACCEPT QMC =  0.967137 +/- 0.000136        |
+--------------------------------------------+
```

These tables provide a concise summary of the simulation results, useful for quick analysis and 
verification of the Quantum Monte Carlo calculations.
