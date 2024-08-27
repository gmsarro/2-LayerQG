# Code for Sarro et al. (2024); Varying Blocking Dynamics with Increased Latent Heating in Idealized 2-Layer QG Models

This repository contains the scripts described in the methods section of Sarro et al. (2024). The scripts are divided into two main categories:

1. **Local Wave Activity and Budget Calculation for a 2-Layer QG Model (with Latent Heating)**
2. **Eddy Growth Rates Calculation for a 2-Layer QG Model**

The scripts are made for the [Lutsko and Hell (2021) model](https://github.com/nicklutsko/moist_QG_channel/tree/main), but can be modified to work on any 2 layer QG model.
The LWA budget here presented work for 2-layer QG models. To compute the LWA budget in more complex datasets, such as reanalysis, see the [falwa page](https://github.com/csyhuang/hn2016_falwa).



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

# Scripts Overview
 1. **Local Wave Activity and Budget Calculation
This category includes scripts for calculating the local wave activity and its budget within a 2-layer quasi-geostrophic (QG) model, with the additional consideration of latent heating effects. These scripts are used to reproduce the results presented in Section 2.3 of Sarro et al. (2024).

# Key Subdirectories:
Make_LWA/: Contains scripts for calculating local wave activity and reference potential vorticity (PV).
Make_LWA_budget/: Includes scripts for calculating reference wind, reference temperature, and the local wave activity budget.

# To get the LWA budget to work compute in this order:
**1) 2-LayerQG/Local_Wave_Activity/Make_LWA/Calculate_LWA.py 
**2) 2-LayerQG/Local_Wave_Activity/Make_LWA_budget/UREF_make.py
**2) 2-LayerQG/Local_Wave_Activity/Make_LWA_budget/LP_budget_make.py or 2-LayerQG/Local_Wave_Activity/Make_LWA_budget/LWA_budget_make.py

 2. ** Eddy Growth Rates Calculation
This category includes scripts to compute the eddy growth rates for the 2-layer QG model. These scripts correspond to the analysis discussed in Section 2.6 of Sarro et al. (2024).
