# SiC-Laser-Annealing-Simulator

A MATLAB-based simulation toolkit for modeling thermal profiles during green laser annealing of implanted Silicon Carbide (SiC) with various capping materials.

## Overview

This project provides a comprehensive simulation framework for studying laser-material interactions in semiconductor manufacturing, specifically focused on the activation of implanted dopants in SiC through laser annealing. The simulator models heat transfer dynamics, temperature profiles, and energy absorption characteristics to optimize annealing parameters.

## Features

- **Multi-layer Material Simulation**: Model complex structures with different material layers including caps like Carbon, AlN, Si₃N₄, TiN, TaN, BN, SiO₂, and Al₂O₃
- **Customizable Laser Parameters**: Adjust wavelength, pulse width, energy density, and timing
- **Dual Laser Support**: Simulate multiple laser exposures with configurable delay times
- **Thermal Analysis**: Generate temperature profiles as a function of depth and time
- **Visualization Tools**: Plot laser power profiles, maximum temperature distributions, and absorption curves
- **Data Export**: Save simulation results to Excel for further analysis

## Getting Started

1. Clone this repository to your local MATLAB environment
2. Run `main.m` to execute a simulation with default parameters
3. Modify material and laser properties in the main function to test different configurations

## Core Functions

- `Calculation.m`: Main computation engine for temperature profiles
- `MatLasProtertyExtraction.m`: Material and laser property extraction
- `LasPowFrofAtTimeExtraction.m`: Laser power profile generation
- `xAxisExtraction.m` and `tAxisExtraction.m`: Mesh generation utilities
- `Compress.m`: Data compression for large outputs

## Usage Example

```matlab
% Define material properties
% [thickness-um, Capacity-J/gC, Density-g/cm3, thermal conductivity-W/cmK]
MatPro = [0.1, 0.97, 2.05, 2.05     % C-film
          0.05, 1.3, 3.21, 0.011    % amf-SiC
          100.0, 1.3, 3.21, 2.8];   % crs-SiC

% Define laser properties
% [wavelength-nm, pulsewidth-us, energydensity-J/cm2, delay-us, mode]
LasPro = [527, 0.25, 1.00, 0.00, 1  % First laser (Gaussian)
          527, 0.25, 1.00, 1.00, 1]; % Second laser (Gaussian)

% Run simulation and analyze results
[TMax, Tatx] = Calculation(MatPro, LasPro, RefAbs, xMesh, tMesh, ProfilAtDepth);
```

## Applications

- Optimizing laser parameters for SiC dopant activation
- Testing different capping materials for enhanced annealing
- Studying thermal profiles at material interfaces
- Investigating dual-laser approaches for controlled heating

## Research Context

This simulator supports research in green laser annealing technology for SiC semiconductor manufacturing, addressing challenges in:

- Finding optimal capping materials with favorable thermal and optical properties
- Determining ideal laser parameters (energy density, pulse width, delay between dual lasers)
- Analyzing temperature profiles at material interfaces
- Comparing different annealing conditions for dopant activation

