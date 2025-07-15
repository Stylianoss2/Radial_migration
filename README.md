# Radial Migration with Spiral Arm Perturbation

This project simulates **stellar radial migration** in a galactic disc under the influence of **logarithmic potentials** and **rotating spiral arm perturbations**. 
It models orbital dynamics using Leapfrog integration in both the inertial and rotating frames, supporting 2D and 3D configurations.

## Features

- Galactic potential modeled via a **logarithmic potential**.
- Inclusion of **spiral arm perturbations** with constant pattern speed (ωₚ).
- Leapfrog integrator adapted for use in a **rotating frame**.
- Computation of conserved quantities: energy, angular momentum, and Jacobi integral.
- Support for:
  - 2D and 3D simulations.
  - Output of full orbital trajectories, accelerations, and derived parameters.

## Requirements

- C++17 or later
- CMake (>= 3.15)
- GCC or Clang
- MPI (optional, if using parallel versions)

## Building the Project
cmake -B build
cmake --build build
./build/bin/LogPotExec -no_pert // for the logarithmic potential only
./build/bin/LogPotExec -pert // To include a perturbation imposed by a spiral arm pattern

