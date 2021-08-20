# Gear restrictions create conservation and fisheries trade-offs for management

This repository accompanies the paper: Carvalho, PG and Humphries, AT. 2021. Gear restrictions create conservation and fisheries trade-offs for management. _Fish and Fisheries_. [10.1111/faf.12607](https://doi.org/10.1111/faf.12607)

## Abstract

Gear-based management for coral reef fisheries is often overlooked in the scientific literature. Empirical studies have demonstrated the conservation benefits of gear-restricted areas (i.e. prohibiting fishing gears), which can support greater biomass than unrestricted areas and protect species that play key functional roles. However, population dynamics of functional feeding groups of reef fishes under specific gear-restriction regimes remains uncertain. Here, we constructed a multi-species, length-based fisheries model to observe relative biomass and catch of reef fishes under various gear-restriction management scenarios. We used fishery-dependent and fishery-independent data to determine the catchability of functional groups and selectivity of size classes for hook-and-line, net and spear fishing, which are widely used gear types on coral reefs globally. Our model revealed trade-offs involved with gear-restriction management such that no single management strategy was able to maximize biomass or catch of all functional groups simultaneously. Also, we found that spear fishing (i.e. prohibiting hook-and-line and net fishing) maintained the highest total biomass summed across functional groups, whilst hook-and-line fishing (i.e. prohibiting net and spear fishing) and a ban on spears maintained the lowest biomass. However, hook-and-line fishing generated the highest catch-per-unit-effort. Our model results were primarily driven by differential growth rates, maximum per capita production of recruits, and catchability of functional groups targeted by each fishing gear. We demonstrate that gear restrictions can be a critical management tool for maintaining biomass and catch of certain functional groups but will likely require additional management to protect all key functional feeding groups of coral reef fishes.

## Usage Notes

### gear_reef_fisheries_model.R

Data: This script reads .csv files provided in the ["data"](https://github.com/paulcarvalho/faf_gear_restrictions_paper/tree/main/data) folder

Supporting files: 
_model_functions.R_ contains functions that are called in _gear_reef_fisheries_model.R_. Documentation (e.g., description, usage, etc.) for functions can be viewed in R using the command docstring("function name"). Requires the [docstring package](https://cran.r-project.org/web/packages/docstring/vignettes/docstring_intro.html).

_run_model.R_