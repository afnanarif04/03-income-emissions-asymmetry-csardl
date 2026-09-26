# =====================================================================
# 01_packages_data.R
# Always runs first. Installs/loads required packages, reads the
# processed panel data from data/panel_data_income_emissions.xlsx, and
# builds the shared panel objects used by every downstream script.
#
# Set the working directory to the repository root before running this
# script (or run via 00_run_all.R, which does this automatically).
# =====================================================================

required_pkgs <- c("tidyverse", "readxl", "openxlsx")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(required_pkgs, library, character.only = TRUE))

select   <- dplyr::select
filter   <- dplyr::filter
rename   <- dplyr::rename
mutate   <- dplyr::mutate
arrange  <- dplyr::arrange
summarise <- dplyr::summarise

# ---------------------------------------------------------------------
# Read the processed panel data (already cleaned; see data/README_data.md
# for variable definitions and original sources).
# ---------------------------------------------------------------------
panel <- read_excel("data/panel_data_income_emissions.xlsx", sheet = "Panel_Data")
lockin_map <- read_excel("data/panel_data_income_emissions.xlsx",
                         sheet = "Carbon_Lockin_Classification")

panel <- panel %>%
  arrange(country_iso3, year) %>%
  mutate(id = as.integer(factor(country_iso3, levels = unique(country_iso3))))

cat("Panel loaded:", n_distinct(panel$id), "countries,",
    min(panel$year), "-", max(panel$year), "(T =",
    max(panel$year) - min(panel$year) + 1, ")\n")
cat("Carbon lock-in groups:", paste(names(table(panel$carbon_lockin_group)),
    table(panel$carbon_lockin_group), sep = "=", collapse = ", "), "\n")

# Rename to the short internal variable names used throughout the paper's
# estimation code, for consistency with the manuscript's notation.
panel <- panel %>%
  rename(iso3 = country_iso3, country = country_name,
         lco2 = ln_co2_per_capita, lgdp = ln_gdp_per_capita,
         lenergy = ln_energy_per_capita, lpop = ln_population,
         lockin = carbon_lockin_group)

KEY  <- "lgdp"
CTRL <- c("lenergy", "lpop")

# Sub-panel index vectors used by 06_subsample_sensitivity.R
oecd_iso3 <- c("AUS","AUT","BEL","CAN","CHL","COL","CRI","CZE","DNK","EST",
               "FIN","FRA","DEU","GRC","HUN","ISL","IRL","ISR","ITA","JPN",
               "KOR","LVA","LTU","LUX","MEX","NLD","NZL","NOR","POL","PRT",
               "SVK","SVN","ESP","SWE","CHE","TUR","GBR","USA")

dir.create("output", showWarnings = FALSE)

cat("01_packages_data.R complete\n\n")
