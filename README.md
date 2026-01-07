# Shimazoe_2026
Source Code for the simulations of the paper 'Linker histone H1 functions as a liquid-like glue to organize chromatin in living human cells'

You can find here our model for chromatin phase separation and analysis scripts with the community. Please use it freely and cite our paper: [link after acceptance]

For questions, contact Jan Huertas: jh2366[at]cam.ac.uk

## Installation Guide
For compiling the chromatin multiscale model, please follow the instructions in https://github.com/CollepardoLab/CollepardoLab_Chromatin_Model/tree/plugin-dev

## Usage
The files needed to run a typical simulation from the paper are in the demo folder. For the simulations of 108 nucleosome fibers, check the demo_108 folder. Please download the whole folder, and run using at least 16 cores, for example:

mpirun -np 16 ./lmp -in run.in

## Visualization
We recommend using OVITO (https://www.ovito.org/) or VMD (https://www.ks.uiuc.edu/Research/vmd) for visualising the trajectory.
