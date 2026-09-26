# =====================================================================
# 03_descriptive_statistics_table1.R
# Produces output/table1_descriptive_statistics.csv (manuscript Table 1).
# =====================================================================
source("code/01_packages_data.R")

desc <- panel %>%
  summarise(
    lco2_mean = mean(lco2, na.rm = TRUE), lco2_sd = sd(lco2, na.rm = TRUE),
    lco2_min = min(lco2, na.rm = TRUE), lco2_max = max(lco2, na.rm = TRUE),
    lgdp_mean = mean(lgdp, na.rm = TRUE), lgdp_sd = sd(lgdp, na.rm = TRUE),
    lgdp_min = min(lgdp, na.rm = TRUE), lgdp_max = max(lgdp, na.rm = TRUE),
    lenergy_mean = mean(lenergy, na.rm = TRUE), lenergy_sd = sd(lenergy, na.rm = TRUE),
    lenergy_min = min(lenergy, na.rm = TRUE), lenergy_max = max(lenergy, na.rm = TRUE),
    lpop_mean = mean(lpop, na.rm = TRUE), lpop_sd = sd(lpop, na.rm = TRUE),
    lpop_min = min(lpop, na.rm = TRUE), lpop_max = max(lpop, na.rm = TRUE)
  )

table1 <- tibble(
  Variable = c("ln(CO2 per capita)", "ln(GDP per capita)",
               "ln(Energy per capita)", "ln(Population)"),
  Mean = c(desc$lco2_mean, desc$lgdp_mean, desc$lenergy_mean, desc$lpop_mean),
  `Std. Dev.` = c(desc$lco2_sd, desc$lgdp_sd, desc$lenergy_sd, desc$lpop_sd),
  Minimum = c(desc$lco2_min, desc$lgdp_min, desc$lenergy_min, desc$lpop_min),
  Maximum = c(desc$lco2_max, desc$lgdp_max, desc$lenergy_max, desc$lpop_max)
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("\n=== Manuscript Table 1 ===\n")
print(table1)
write_csv(table1, "output/table1_descriptive_statistics.csv")

cat("\n03_descriptive_statistics_table1.R complete\n\n")
