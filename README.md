# Optimal Ground-based Anchor Placement

This repository contains the MATLAB scripts used to reproduce the figures and compute optimal anchor placements in the paper:

J. Ko et al., “Optimal Ground-based Anchor Placements for Least-Square Multilateration,” Electronics Letters, 2025.

## Repository Structure

- README.md
- gen_fig.m
  - Generates the figures (Fig 1-a, 1-b, 2-a, 2-b) shown in the paper. Figures 3 and 4 are omitted because 'BoundaryConstraintOptimization.m' can replace the role.
- NormConstraintOptimization.m
  - Computes the optimal anchor positions under the norm-constraint case
- BoundaryConstraintOptimization.m
  - Computes the optimal anchor positions under the boundary-constraint case

## Prerequisites

- In this paper, we use matlab R2024a. 
- To run the code, you need matlab toolbox whose name is'Optimization Toolbox'.

## Usage
1. Clone the repository
```bash
git clone https://github.com/KO-JUNSUNG/OAP.git
cd OAP
```

2. You can change the number of anchors by changing the 'N'. 
3. If you want to change the norm or radius constraint then change the constraint function.