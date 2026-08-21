# STTE-MDR-Code

Implementation and simulation-analysis code for comparing inverse-probability, doubly robust, model-based density-ratio (MDR), trial-emulation, and longitudinal TMLE estimators of the sustained-treatment effect (STTE).

**Authors:** Zern Ke (1), Mingshi Cui (1), Feng Dai (2), Birol Emir (2), Javier Cabrera (1), and Demissie Alemayehu (2)

**Affiliations:**

1. Department of Statistics, Rutgers University, New Brunswick
2. Pfizer Inc, New York, NY, USA

## Estimators

The simulation compares the following ten estimators:

| Code name | Display name |
|---|---|
| `ipcw_pool` | IPCW-Entry |
| `dr_ipcw_pool` | DR-IPCW-Entry |
| `ipcw_state_gcomp` | IPCW-State |
| `dr_ipcw_state_gcomp` | DR-IPCW-State |
| `mdr_state_gcomp_follow_balanced_cell_normalized` | MDR-State |
| `mdr_pool_follow_balanced_cell_normalized` | MDR-Entry |
| `dr_mdr_state_gcomp_follow_balanced_cell_normalized` | DR-MDR-State |
| `dr_mdr_pool_follow_balanced_cell_normalized` | DR-MDR-Entry |
| `trial_emulation` | TrialEmulation |
| `ltmle_tmle_glm` | LTMLE-GLM |

All four MDR estimators use the follow-up-balanced, cell-normalized MDR construction reported in the paper.

## Repository contents

| File | Purpose |
|---|---|
| `stte_simulation_all_parallel.R` | Main parallel driver for the paper simulation. |
| `stte_simulation_all.R` | Coordinates the ten estimators and provides serial simulation helpers. |
| `stte_per_protocol_simulation.R` | Data-generating process, target-trial construction, non-MDR estimators, summaries, and plotting helpers. |
| `mdr_weight.R` | MDR weight estimation used by the main driver. |
| `stte_doubly_robust_helpers.R` | Doubly robust estimation helpers. |
| `ltmle.R` | Longitudinal TMLE data preparation and estimation. |
| `save_observed_data.R` | Generates and saves replication-level observed datasets for the scenario-diagnostics workflow. |
| `analyze_simulation_performance.R` | Produces the paper's performance tables and figures from a separately generated `result500.rds`. |
| `analyze_simulation_scenario_diagnostics.R` | Produces scenario-diagnostic tables and figures from `observed_data500.rds`. |

## Software requirements

R 4.4 or later is recommended. Install the required contributed packages with:

```r
install.packages(c(
  "data.table",
  "dplyr",
  "ggplot2",
  "ltmle",
  "patchwork",
  "randomForest",
  "ranger",
  "scales",
  "tidyr",
  "tidyselect",
  "TrialEmulation"
))
```

The code was checked with the following package versions:

| Package | Version |
|---|---:|
| data.table | 1.17.8 |
| dplyr | 1.1.4 |
| ggplot2 | 4.0.3 |
| ltmle | 1.3-0 |
| patchwork | 1.3.1 |
| randomForest | 4.7-1.2 |
| ranger | 0.17.0 |
| scales | 1.4.0 |
| tidyr | 1.3.1 |
| tidyselect | 1.2.1 |
| TrialEmulation | 0.0.4.11 |

## Paper simulation settings

The defaults in `stte_simulation_all_parallel.R` reproduce the paper configuration:

- Seven scenarios: base, small sample, strong confounding, poor positivity, poor adherence, rare event, and nonlinear.
- 500 replications per scenario.
- Cohort size `n = 2000`, except `n = 100` in the small-sample scenario.
- Truth calculation based on `n_truth = 50000`.
- Random seed `20260630`.
- Simulated MDR target sample ten times the number of eligible person-trials.
- Strategy-specific `ranger` probability forests with 500 trees.
- Up to ten parallel workers, subject to the available cores.

## Reproducing the simulation

Run R from the repository directory and execute:

```r
source("stte_simulation_all_parallel.R")

reproduced_results <-
  run_all_methods_mdr_solutions_simulation_parallel()

saveRDS(reproduced_results, "result500.rds")
```

The complete run evaluates 3,500 simulated cohorts and all ten estimators. It is computationally intensive and requires substantially more time and memory than a smoke test.

The large simulation-output file `result500.rds` is not included in this repository. Running the code above creates it with 35,000 scenario-method-replication rows: seven scenarios, ten estimators, and 500 replications per combination. Estimator failures are retained in the result object and summarized through `n_success` and `failure_rate`.

## Reproducing the performance tables and figures

After generating or separately obtaining `result500.rds`, place it in the repository directory and run:

```sh
Rscript analyze_simulation_performance.R
```

The script writes CSV tables and PDF/PNG figures beside the analysis script, including:

- `table_simulation_performance.csv`
- `figure_simulation_performance_heatmap.pdf`
- `table_simulation_weight_diagnostics.csv`
- `figure_simulation_weight_diagnostics.pdf`
- Weight-maximum and weight-mean distribution tables and figures

To read and write these files in another directory, set the `STTE_MDR_RESULT_DIR` environment variable to the directory containing `result500.rds`.

## Scenario diagnostics and observed data

`observed_data500.rds` is not included in this repository because the file is approximately 1.11 GB. The file contains 500 observed datasets for each of the seven scenarios and is only required to reproduce the scenario-diagnostics analysis.

The main paper simulation uses master seed `20260630`. The observed data used for the scenario diagnostics were generated in six independent batches. Seeds `20260630` and `42` generated 50 replications each; seeds `20260731`, `666666`, `88667`, and `688202` generated 100 replications each. The batches were combined in that order to obtain replications 1--500 within every scenario.

`save_observed_data.R` implements the replication-level observed-data generation step and saves each dataset as an RDS file. Its executable block generates the final 100-replication batch for all seven scenarios using seed `688202`. Other batches can be generated by sourcing the script and calling `save_generated_observed_data_all_methods_mdr_solutions()` with the corresponding seed, replication count, and output directory.

If `observed_data500.rds` is obtained separately, place it beside `analyze_simulation_scenario_diagnostics.R` and run:

```sh
Rscript analyze_simulation_scenario_diagnostics.R
```

Outputs are written to the `scenario_diagnostics/` directory.

## Citation

If you use this code, please cite both the accompanying manuscript and the software repository. Replace the remaining bracketed fields below with the final publication and repository information.

### Manuscript

> Ke, Z., Cui, M., Dai, F., Emir, B., Cabrera, J., & Alemayehu, D. (2026). *From Cumulative Weights to Marginal Density Ratios: Per-Protocol Estimation in Sequential Target Trials*. [Journal, preprint server, or other publication information].

```bibtex
@unpublished{ke2026marginaldensityratios,
  author = {Ke, Zern and Cui, Mingshi and Dai, Feng and Emir, Birol and Cabrera, Javier and Alemayehu, Demissie},
  title  = {From Cumulative Weights to Marginal Density Ratios: Per-Protocol Estimation in Sequential Target Trial Emulation},
  year   = {2026},
  note   = {[Manuscript status, journal, or preprint information]},
  doi    = {[DOI, if available]}
}
```

### Software repository

> Ke, Z., Cui, M., Dai, F., Emir, B., Cabrera, J., & Alemayehu, D. (2026). *STTE-MDR-Code: Implementation and simulation-analysis code for “From Cumulative Weights to Marginal Density Ratios”* [Computer software]. GitHub. [Repository URL]

```bibtex
@software{ke2026sttemdrcode,
  author  = {Ke, Zern and Cui, Mingshi and Dai, Feng and Emir, Birol and Cabrera, Javier and Alemayehu, Demissie},
  title   = {{STTE-MDR-Code}: Implementation and Simulation-Analysis Code for ``From Cumulative Weights to Marginal Density Ratios: Per-Protocol Estimation in Sequential Target Trial Emulation''},
  year    = {2026},
  url     = {[Repository URL]},
  version = {[Release or commit identifier]}
}
```

## License

This repository is released under the [MIT License](LICENSE). Copyright is held by the authors listed above.

## Reproducibility notes

- Run scripts from the repository directory or use their absolute paths.
- The small-sample cohort size is applied automatically by the `small_size` scenario definition.
- The main paper implementation is the parallel driver; `stte_simulation_all.R` is sourced as a dependency and also provides serial helpers.
- Temporary TrialEmulation model files are written under the R session's temporary directory unless another model directory is supplied.
