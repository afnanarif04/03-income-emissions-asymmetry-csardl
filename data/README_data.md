# Data Documentation

## File: `panel_data_income_emissions.xlsx`

Three sheets: `Panel_Data` (the estimation panel), `Carbon_Lockin_Classification`
(the country-level carbon lock-in group assignment, used in Table 3's
sub-panels), and `Income_Group_Classification` (the country-level income
tercile assignment, used in Table 4's sub-panels).

### Sheet: `Panel_Data`

| Column | Description | Unit | Source |
|---|---|---|---|
| `country_iso3` | ISO 3166-1 alpha-3 country code | — | Our World in Data |
| `country_name` | Country name | — | Our World in Data |
| `year` | Calendar year | — | — |
| `co2_per_capita` | CO2 emissions per capita | tonnes | Our World in Data CO2 dataset |
| `energy_per_capita` | Primary energy consumption per capita | kWh | Our World in Data Energy dataset |
| `gdp_current_usd` | Gross domestic product, current US dollars | USD | Our World in Data (Maddison Project / World Bank via OWID) |
| `population` | Total population | persons | Our World in Data |
| `gdp_per_capita` | GDP per capita, current US dollars | USD | Derived: `gdp_current_usd / population` |
| `coal_share_energy_pct` | Share of primary energy from coal | % | Our World in Data Energy dataset |
| `gas_share_energy_pct` | Share of primary energy from natural gas | % | Our World in Data Energy dataset |
| `oil_share_energy_pct` | Share of primary energy from oil | % | Our World in Data Energy dataset |
| `fossil_share_energy_pct` | Share of primary energy from all fossil sources | % | Our World in Data Energy dataset |
| `renewables_share_energy_pct` | Share of primary energy from renewables | % | Our World in Data Energy dataset |
| `ln_co2_per_capita` | Natural log of `co2_per_capita` | — | Derived |
| `ln_gdp_per_capita` | Natural log of `gdp_per_capita` | — | Derived |
| `ln_energy_per_capita` | Natural log of `energy_per_capita` | — | Derived |
| `ln_population` | Natural log of `population` | — | Derived |
| `carbon_lockin_group` | Carbon lock-in classification (see below) | — | Derived |
| `income_group` | Income tercile classification (see below) | — | Derived |

The four `ln_*` columns are the variables actually used in estimation
(`lco2`, `lgdp`, `lenergy`, `lpop` in the code). They have zero missing
values across the full balanced panel. The energy-mix share columns have
missing values for some country-years where Our World in Data does not
report a full mix breakdown; these are not used directly in estimation and
only inform the carbon lock-in classification, which is itself computed
without missing values (see below).

The `income_group` column is provided for transparency only; the estimation
code (`code/07_subsample_sensitivity_table4.R`) recomputes this
classification directly from `gdp_per_capita` and does not read this column
back in, so Table 4 remains reproducible from the underlying GDP data alone
even if this column were absent.

### Sheet: `Carbon_Lockin_Classification`

One row per country. `carbon_lockin_group` takes one of four values:
`coal`, `gas`, `oil`, or `unclassified`. A country is classified as
dominated by a given fossil fuel if that fuel's time-averaged share of its
energy mix exceeds a threshold; countries that do not meet the threshold
for any single fuel are `unclassified` and are retained only in the
full-panel ("All countries") estimation, not in any carbon lock-in
sub-panel.

| Column | Description |
|---|---|
| `country_iso3` | ISO 3166-1 alpha-3 country code |
| `country_name` | Country name |
| `coal_share_avg_pct` | Time-averaged coal share of energy mix |
| `gas_share_avg_pct` | Time-averaged gas share of energy mix |
| `oil_share_avg_pct` | Time-averaged oil share of energy mix |
| `n_energy_sources_used` | Number of energy sources with non-missing data used in the classification |
| `carbon_lockin_group` | Final classification: coal / gas / oil / unclassified |

### Sheet: `Income_Group_Classification`

One row per country, sorted by time-averaged GDP per capita. Countries are
split into three equally sized groups (terciles) on `avg_gdp_per_capita`:
the bottom third are `low`, the middle third are `middle`, and the top
third are `high`.

| Column | Description |
|---|---|
| `country_iso3` | ISO 3166-1 alpha-3 country code |
| `country_name` | Country name |
| `avg_gdp_per_capita` | Time-averaged GDP per capita, current USD (1990-2022 mean) |
| `income_group` | Final classification: low / middle / high |

This sheet is provided for transparency only (see the note under
`Panel_Data` above); `code/07_subsample_sensitivity_table4.R` recomputes the
same classification directly from `Panel_Data`'s `gdp_per_capita` column.

---

## Sample

- **N (countries):** 153
- **T (years):** 33 (1990-2022)
- **Panel:** balanced (5,049 country-year observations)
- **Carbon lock-in sub-panels** (used in Table 3): coal-dominant = 10
  countries, gas-dominant = 18 countries, oil-dominant = 48 countries,
  unclassified = 77 countries
- **Income tercile sub-panels** (used in Table 4): low-income = 51
  countries, middle-income = 51 countries, high-income = 51 countries,
  split at the 33rd and 66th percentiles of time-averaged GDP per capita
  (cutoffs approximately USD 4,802 and USD 16,004)

## Data sources

| Source | URL | Variables | Notes |
|---|---|---|---|
| Our World in Data — CO2 and Greenhouse Gas Emissions | https://github.com/owid/co2-data | `co2_per_capita` | Downloaded as `owid_co2_data.csv` in the original project; merged and cleaned into the panel provided here |
| Our World in Data — Energy | https://github.com/owid/energy-data | `energy_per_capita`, `gdp_current_usd`, `population`, energy-mix shares | Downloaded as `owid_energy_data.csv` in the original project; merged and cleaned into the panel provided here |

The raw downloaded CSV files are not included in this repository; the
processed panel (`panel_data_income_emissions.xlsx`) is the analysis-ready
dataset actually used to produce every result in the paper, and is
sufficient on its own to run every script in `code/`.

## Notes on missing values and panel balance

- The core estimation variables (`ln_co2_per_capita`, `ln_gdp_per_capita`,
  `ln_energy_per_capita`, `ln_population`) have no missing values; the panel
  is fully balanced across all 153 countries and all 33 years for these
  variables.
- The energy-mix share columns (`coal_share_energy_pct` etc.) have missing
  values for some country-years, reflecting gaps in the underlying Our World
  in Data energy-mix reporting. These columns are not used directly in any
  estimator; they only informed the (separately provided) carbon lock-in
  classification.
