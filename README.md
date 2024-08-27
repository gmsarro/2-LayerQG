# Reproducible code for "Sarro et al. 2024"
# Varying Blocking Dynamics with Increased Latent Heating in Idealized 2-Layer QG Models

This repository contains the scripts used in the methods section of the paper "Sarro et al. 2024". The scripts are divided into two main categories:

1. **Local Wave Activity and Budget Calculation for a 2-Layer QG Model (with Latent Heating)**
2. **Eddy Growth Rates Calculation for a 2-Layer QG Model**

## Repository Structure

```bash
2-LayerQG/
│
│
│   ├── Local_Wave_Activity/       # Scripts for local wave activity and budget calculation
│   │   ├── Make_LWA/              # Scripts for local wave activity calculation + reference PV
│   │   ├── Make_LWA_budget/       # Scripts for reference wind + reference temperature + local wave activity budget calculation
│   ├── Eddy_Growth/   # Scripts for eddy growth rates calculation
```

## Scripts Overview
# 1. **Local Wave Activity and Budget Calculation
This category includes scripts for calculating the local wave activity and its budget within a 2-layer quasi-geostrophic (QG) model, with the additional consideration of latent heating effects. These scripts are used to reproduce the results presented in Section [X] of Sarro et al. (2024).

# Key Subdirectories:
Make_LWA/: Contains scripts for calculating local wave activity and reference potential vorticity (PV).
Make_LWA_budget/: Includes scripts for calculating reference wind, reference temperature, and the local wave activity budget.
# 2. Eddy Growth Rates Calculation
This category includes scripts to compute the eddy growth rates for the 2-layer QG model. These scripts correspond to the analysis discussed in Section [Y] of Sarro et al. (2024).

# Key Subdirectory:
Eddy_Growth/: Contains scripts for determining the eddy growth rates for the 2-layer QG model.
