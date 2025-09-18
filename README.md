# graphene
Spectral calculations for multi-layer graphene as described in [QAE Phases in gated rhobohedral graphene](https://arxiv.org/abs/2509.05439). 

## spectrum_calc_scan.m
Use this script to calculate spectra. Requires MATLAB parallel computing toolbox as written. To remove this toolbox requirement change the annotated lines in "residue_map.m" and "grid_adapt.m". 

## graphene_eps.m, graphene_om.m
Functions which define the appropriate interface operator (operator $A$ as described). Change parameter values for number of layers, $\gamma$, and $u$ as described in above reference inside of these functions ($\gamma=$ eps, $u$ = om). 

## make_H.m
Use to generate and diagonalize bulk (constant coefficient) Hamiltonians of the system.

## test_modes_vec.m
Can be used to extract eigenvectors at a particular $(k_x, E)$ value.
