# Replication Package — Asymmetric Income-Emissions Adjustment in a Global Panel

Replication code and data for an applied panel-econometrics study examining
whether carbon dioxide emissions respond asymmetrically to rising versus
falling income, whether this asymmetry depends on which fossil fuel
dominates a country's energy supply (carbon lock-in), and whether it holds
uniformly across income levels and is robust to small-sample bias
correction.

This repository is anonymised for peer review. Author name, affiliation, and
paper title are withheld until after acceptance.

---

## Software requirements

- R version 4.3.3 or later
- Packages: `tidyverse` (2.0.0+), `readxl` (1.5.0+), `openxlsx` (4.2.8.1+)

Install with:
```r
install.packages(c("tidyverse", "readxl", "openxlsx"))
```

All required packages are also auto-installed by `code/01_packages_data.R` if
missing.

---

## Data

Balanced panel of 153 economies, annual observations 1990-2022 (T = 33).
See `data/README_data.md` for the full variable list, units, and original
data sources (Our World in Data CO2 and Energy datasets). The processed,
analysis-ready panel is provided directly as
`data/panel_data_income_emissions.xlsx` — this is not the raw downloaded
data, but the cleaned and merged panel actually used in estimation. All
original variables (levels, before the natural-log transformation used in
estimation) are retained in the file for transparency.

---

## Repository structure

```
.
├── README.md                                    (this file)
├── LICENSE
├── data/
│   ├── README_data.md
│   └── panel_data_income_emissions.xlsx
├── code/
│   ├── 00_run_all.R                              master runner
│   ├── 01_packages_data.R                        loads packages and data
│   ├── 02_estimator_functions.R                  estimator library (sourced by 03-10)
│   ├── 03_descriptive_statistics_table1.R        -> Table 1
│   ├── 04_preliminary_tests_table2.R             -> Table 2 (4 panels: CSD, unit root, structural-break unit root, cointegration)
│   ├── 05_main_estimation_table3.R               -> Table 3, Panel A + Wald/residual-CD rows
│   ├── 06_supplementary_dl_and_pmg.R             -> Table 3, Panels B and C (all 4 samples)
│   ├── 07_subsample_sensitivity_table4.R         -> Table 4 (income-tercile sub-panels)
│   ├── 08_temporal_stability_check.R             diagnostic only, not a manuscript table
│   ├── 09_csa_joint_test_table3.R                -> Table 3, "CSA terms" row (all 4 samples)
│   └── 10_bias_corrected_table5.R                -> Table 5 (RMA/HPJ bias-corrected estimates)
└── output/
    └── (created by the scripts; one file per manuscript table)
```

---

## Replication instructions

**Option 1 — run everything at once:**
```r
setwd("path/to/this/repository")
source("code/00_run_all.R")
```
Total runtime: approximately 6-7 minutes on a standard laptop (tested on R 4.3.3,
Ubuntu 24.04).

**Option 2 — run scripts individually**, in numeric order, from the repository
root:
```r
setwd("path/to/this/repository")
source("code/01_packages_data.R")
source("code/02_estimator_functions.R")
source("code/03_descriptive_statistics_table1.R")
source("code/04_preliminary_tests_table2.R")
source("code/05_main_estimation_table3.R")
source("code/06_supplementary_dl_and_pmg.R")
source("code/07_subsample_sensitivity_table4.R")
source("code/08_temporal_stability_check.R")
source("code/09_csa_joint_test_table3.R")
source("code/10_bias_corrected_table5.R")
```
Scripts 03 and onward each call `source("code/01_packages_data.R")`
internally, so they can also be run standalone in any order after 01 and 02
have been sourced once.

No random-number generation is used anywhere in this pipeline (all
estimators are closed-form OLS or a deterministic numerical optimiser for the
pooled mean group estimator), so no seed is required for exact
reproducibility.

---

## Output files → manuscript tables

| Script | Output file(s) | Manuscript table |
|---|---|---|
| `03_descriptive_statistics_table1.R` | `output/table1_descriptive_statistics.csv` | Table 1 |
| `04_preliminary_tests_table2.R` | `output/table2_prelim_tests_full.json`, `output/table2_prelim_tests_full.rds` | Table 2 (four panels: cross-sectional dependence; first/second-generation unit root; third-generation unit root with structural breaks; cointegration and slope homogeneity) |
| `05_main_estimation_table3.R` | `output/table3a-c_*.csv`, `output/table3d_pmg_by_panel_detail.csv` | Table 3, Panel A (all 4 samples) + Wald symmetry/residual CD rows |
| `06_supplementary_dl_and_pmg.R` | `output/table3_panelB_supplementary_dl.csv`, `output/table3_panelB_wald_symmetry.csv`, `output/table3_panelC_pmg.csv` | Table 3, Panels B and C (all four samples, including each panel's own Wald symmetry test) |
| `07_subsample_sensitivity_table4.R` | `output/table4_income_terciles.csv`, `output/table4b_income_terciles_csa_joint_f.csv`, `output/table4c_income_terciles_wald_symmetry.csv`, `output/table4d_income_terciles_wald_symmetry_dl.csv`, `output/table4_income_group_assignment.csv` | Table 4 |
| `09_csa_joint_test_table3.R` | `output/table3e_csa_joint_f_main_panels.csv` | Table 3, "CSA terms" row (all 4 samples) |
| `10_bias_corrected_table5.R` | `output/table5_rma_results.csv`, `output/table5_hpj_results.csv`, `output/table5_wald_rma.csv`, `output/table5_wald_hpj.csv` | Table 5 |
| `08_temporal_stability_check.R` | `output/supplementary_temporal_stability_not_a_table.csv` | Not a manuscript table — referenced only in the Conclusion's Limitations discussion |

---

## Estimators

| Code label | Estimator | Key citation |
|---|---|---|
| CS-ARDL | Cross-sectionally augmented ARDL | Chudik & Pesaran (2015) |
| CS-NARDL | Asymmetric CS-ARDL (partial-sum decomposition) | Shin, Yu & Greenwood-Nimmo (2014); cross-sectional averages built from the partial sums per Kapetanios, Mitchell & Shin (2014) |
| FCS-NARDL | CS-NARDL with low-frequency Fourier augmentation | Enders & Lee (2012a) |
| CS-DL | Cross-sectionally augmented distributed lag (symmetric, no lagged dependent variable) | Chudik, Mohaddes, Pesaran & Raissi (2016) |
| CS-NDL | Asymmetric CS-DL | Combines Shin et al. (2014) with Chudik et al. (2016) |
| FCS-NDL | Fourier-augmented CS-NDL | — |
| PMG-ARDL | Pooled mean group (symmetric) | Pesaran, Shin & Smith (1999); estimated here by concentrated maximum likelihood, following Blackburne & Frank (2007) |
| PMG-NARDL | Asymmetric PMG | — |
| FPMG-NARDL | Fourier-augmented asymmetric PMG | — |
| CS-ARDL-RMA / CS-NARDL-RMA / FCS-NARDL-RMA | CS-ARDL group with recursive mean adjustment | Bias-correction transform demeans the dependent variable by its own expanding-window mean before estimation |
| CS-ARDL-HPJ / CS-NARDL-HPJ / FCS-NARDL-HPJ | CS-ARDL group with half-panel jackknife | Bias correction: 2 x full-sample estimate minus half the sum of the two half-sample estimates, applied to long-run coefficients only |

CS-DL, CS-NDL, and FCS-NDL are reported only as supplementary robustness
checks and are not proposed as a contribution of this study. PMG-ARDL,
PMG-NARDL, and FPMG-NARDL serve as a slope-homogeneity benchmark against the
heterogeneous, cross-sectionally augmented mean-group estimators that form
the primary contribution. The RMA and HPJ variants in Table 5 test whether
the paper's headline findings on energy-use dominance and long-run income
symmetry are sensitive to Nickell-type small-sample bias in the lagged
dependent variable.

---

## Notes

- Tested on R 4.3.3, Ubuntu 24.04.
- Estimated total runtime for the full pipeline: approximately 6-7 minutes, most of which is Table 2's two bootstrap-based third-generation unit root tests (CADF-SB and BCIPS-BCADF, B=99 replications each, run for all four variables).
- No simulation or bootstrap procedure is used; all reported estimates are
  exactly reproducible from the provided data with no random seed required.
- The carbon lock-in classification (coal-, gas-, oil-dominant, or
  unclassified) is included in the processed data file and is not
  recomputed by this code; see `data/README_data.md` for its construction.
- The income-group classification (low-, middle-, high-income terciles,
  used in Table 4) is computed directly by `code/07_subsample_sensitivity_table4.R`
  from each country's time-averaged GDP per capita; it is also included in
  the processed data file's `income_group` column and the
  `Income_Group_Classification` sheet for transparency, but the code does
  not read either of these back in — it recomputes the classification
  itself so the table is reproducible from the underlying GDP data alone.
