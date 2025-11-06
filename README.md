
# Overview

This repository contains and data and scripts for reproducing the results accompanying the manuscript  

### Chronic infections can generate SARS-CoV-2-like bursts of viral evolution without epistasis
Edwin Rodríguez-Horta<sup>1,2</sup>, John Strahan<sup>3</sup>, Aaron R. Dinner<sup>3,#</sup> and John P. Barton<sup>1,#</sup>

<sup>1</sup> Department of Computational and Systems Biology, University of Pittsburgh School of Medicine  
<sup>2</sup> Group of Complex Systems and Statistical Physics, Department of Theoretical Physics, Physics Faculty, University of Havana    
<sup>3</sup> Department of Chemistry and James Franck Institute, University of Chicago  

<sup>#</sup> correspondence to [dinner@uchicago.edu](mailto:dinner@uchicago.edu) and [jpbarton@pitt.edu](mailto:jpbarton@pitt.edu)  

The preprint is available at __LINK PENDING__.


# Contents
- `src/` — Julia scripts implementing the evolutionary model and analysis functions  
- `notebooks/Run_Test.ipynb` — Demonstration notebook for running the simulation and generating example plots  
- `data/` — Simulation results and input files for reproducing figures in the paper  

### Software dependencies

- **Julia 1.10 or higher** (https://julialang.org/)  
- **IJulia** to run the notebooks  
- **Python 3** for figures
  
# Installation and Environment Setup

1. Clone this repository and navigate into it:
```bash
git clone https://github.com/bartonlab/paper-SARS-CoV-2-evolution.git
cd paper-SARS-CoV-2-evolution
```

2. Start a Julia session in the project folder. You can either:
```bash
julia --project=.
```
or start Julia normally and activate the environment:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
This will download all dependencies specified in `Manifest.toml`.

3. To run the notebooks:
```julia
using IJulia
IJulia.notebook()
```
This opens Jupyter in your browser. Navigate to `notebooks/` and open `Run_Test.ipynb` to run an example simulation.

# Usage

This repository provides **two main functions**:

1. **`evolve_pop_EpiModel(case0::individual, tEnd::Int64, N::Int64, params::Param; chronical_cases=true)`**  
   Simulates population evolution over time in an epidemiological model that combines within-host replication and between-host transmission, accounting for chronic cases.
   
   - **Arguments:**  
     - `case0::individual`: initial infected individual  
     - `tEnd::Int64`: duration of simulation  
     - `N::Int64`: intra-host population size  
     - `params::Param`: model parameters (mutation rates, selection coefficients, generation times, probability of chronic infection)  
     - `chronical_cases::Bool=true`: include chronic infections (default)  
   - **Returns:** `(pop_final, num_sick, rare_events, chronical_times)`  
     - `pop_final`: population history over time  
     - `num_sick`: infected individuals per time step  
     - `rare_events`: chronic cases per time step  
     - `chronical_times`: times of chronic infections

3. **`av_number_of_mut_variant(pop_final::Vector{Vector{individual}})`**  
   Computes the **average number of mutations per variant** from a population returned by `evolve_pop_EpiModel`. Returns an array of mutation fractions per variant over time.

**Other useful functions:**
- `slope_change_score(fraction_mut_variant_res::Array)` – detects slope changes in mutation accumulation  
- `init_individual(...)` – initializes the first infected individual with given parameters  
- `savitzky_golay(y::Array, window::Int, order::Int)` – smooths time series data  

---

### Example workflow
```julia
# Load selection coefficients and fit distributions
selection_coeff = readdlm("data/selection_coeff_SC2.txt")[:]
bf_effect_dist = fit(LogNormal, selection_coeff[selection_coeff .> 0.02])
nt_effect_dist = fit(Normal, selection_coeff[-0.02 .< selection_coeff .< 0.02])

# Define parameters
parameters = Param(1e-3, 1e-4, bf_effect_dist, nt_effect_dist, 2.0, 4.9, 0.5, 1.6e-4)

# Initialize first individual
case0 = init_individual(1, 0.0, parameters.t_s, 1000, 0)

# Simulate population
pop_final, num_sick, rare_events, chronical_times = evolve_pop_EpiModel(case0, 1000, 1000, parameters, chronical_cases=true)

# Compute mutation fractions and slopes
fraction_mut_variant_res = av_number_of_mut_variant(pop_final)
slope_vs_t, score = slope_change_score(fraction_mut_variant_res[1:end])
```

This workflow generates the population dynamics, tracks chronic infections, and allows further analysis or plotting.



# License

This repository is dual licensed as [GPL-3.0](LICENSE-GPL) (source code) and [CC0 1.0](LICENSE-CC0) (figures, documentation, and our presentation of the data).
