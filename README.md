# Replication of Lee et al. (2026)
Coursework for PSZB17-214 course (Research Field Practice in Cognitive/Affective Psychology). 

Materials, data and code used in the original article are available at [https://doi.org/10.17605/OSF.IO/VUQFN](https://doi.org/10.17605/OSF.IO/VUQFN). 

The original article: 
*Lee, PS., Ballard, T. & Sewell, D.K. Conflict Dynamics Within and Across Trials: The Effect of Distractor-Target Stimulus Onset Asynchrony. Comput Brain Behav (2026). https://doi.org/10.1007/s42113-026-00330-7.*

This repository currently contains the R analysis scripts that test the presence of the CSE. I am re-analyzing the data available on OSF, for the aforementioned course. 

The implementation of the revised diffusion model for conflict tasks (Lee and Sewell, 2024) in R is a work in progress, with availability expected in the future.
*(Lee, P.-S., & Sewell, D. K. (2024). A revised diffusion model for conflict tasks. Psychonomic Bulletin & Review, 31(1), 1–31. https://doi.org/10.3758/s13423-023-02288-0)*

### Cargo crates: 
- `rayon`: parallel trial simulation
- `rand` + `rand_distr`: RNG and normal draws
- `serde`, `toml`, `csv`: configs and I/O
- `clap`: command_line arguments
- `anyhow`: error handling
- `argmin`: optional, for Nelder-Mead

### Filesystem:
`rdmc-fit/
├── Cargo.toml
├── .gitignore            # target/, results/
├── data/                 # observed RT quantiles, error rates (CSV)
├── configs/              # starting values, bounds, n_trials (TOML)
├── results/              # fitted params, simulated CDFs, logs
├── src/
│   ├── lib.rs
│   ├── main.rs           # CLI entry point
│   ├── model.rs          # RDMC simulator (Euler steps)
│   ├── data.rs           # load data, compute quantile bins
│   ├── objective.rs      # G² / chi-square over 6 RT bins + error bin
│   └── optim.rs          # Nelder-Mead (SIMPLEX) loop
├── benches/              # criterion
└── tests/                # e.g. checks, tests of parameters`
