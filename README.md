# PANDA – Predicting Angle from Nanoscale Density Analysis

![Logo](logo.png)

This repository provides tools for calculating the contact angle of a droplet using a one-dimensional density profile. The code requires several input files generated with GROMACS. Below is an overview of the process and the necessary files.

## Requirements

### Input Files

To perform the contact angle calculations, the following GROMACS-generated files are needed:

- .gro - A GROMACS structure file containing the molecular coordinates and box dimensions of the system **after** simulation.
- .xtc - A trajectory file containing the time evolution of the system's coordinates. This file is essential for density profile calculations across different time steps.
- .tpr - A portable binary run file that contains all information about the simulation (topology, parameters, etc.), which allows GROMACS to analyze the trajectory properly.

### Example GROMACS Density Profile Command

To calculate the density profile, use the built-in `gmx density` command. Below is an example command to compute the density profile along the z-axis:
```
gmx density -s *topology*.tpr -f *trajectory*.xtc -o *density_profile*.xvg -sl 200 -d Z -b *initial_time*
```

Here:
- -s specifies the input topology (`.tpr`) file.
- -f specifies the input trajectory (`.xtc`) file.
- -o specifies the output file for the density profile data in .xvg format.
- -sl 200 sets the number of slices to 200 along the chosen axis (for better spatial resolution).
- -d Z indicates that the density profile is calculated along the z-axis.
- -b sets the starting time (in ps) for analysis, allowing you to exclude the initial equilibration phase from the calculation. In practice, it is worth choosing this time so that averaging over the last 5 ns or ~2500 frames is performed.


## Contact Angle Calculation

Two methods are available for calculating the contact angle of the droplet:

- `profile_approx`: Standard method for calculating contact angle from one-dimensional density profile based on the algorithm presented in the original paper
- `profile_approx_modified`: An enhanced version of profile_approx that incorporates calculation of pore parameters from the molecular coordinatese and the addition of another hyperparameter – offset, resulting in a more accurate result when calculating the angle

The choice of method depends on the purpose. The first method is verified and has been published in a journal. However, due to the non-ideal model, the contact angle value may not be accurate. The second method was developed empirically and does not fit the idea of the whole project. However, until general formulas are derived that take into account the non-zero wetting layer, this algorithm is more accurate. Later, after further research, the second method will be retracted

## Examples

In the examples folder, you will find a Jupyter Notebook demonstrating a complete pipeline for contact angle calculations. The notebook includes:

1. Examples of typical constants that are used in calculations.
2. Using two provided methods to calculate the contact angle.
3. Example of density profile visualization and its approximation.
