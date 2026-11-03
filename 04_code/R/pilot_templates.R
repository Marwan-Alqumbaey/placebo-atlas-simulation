# Table structures from the pilot products (structure only; no outcome data).
# Arm sizes and row-inclusion rules of each pilot product's main adverse event table.
# Placebo-controlled: Beyonttra, Leqembi, Pyrukynd, Ruconest. Active-controlled tables
# (Xgeva, Onerji) are used only as size/rule templates. Drovelis (single-arm) and the
# three products without a usable table (Hexacima, Stelara, Peyona) give no template.
pilot_templates <- data.frame(
  product = c("Beyonttra", "Leqembi", "Pyrukynd", "Ruconest", "Xgeva", "Onerji"),
  ND = c(421, 898, 40, 44, 2841, 128),
  NP = c(211, 897, 39, 29, 2836, 131),
  rule = c("either", "either", "none", "none", "either", "either"),
  cut = c(10, 5, NA, NA, 5, 5),
  ge = c(FALSE, TRUE, NA, NA, TRUE, TRUE),
  source = c("AR Table 23: >10% in either group", "AR Table 38: at least 5% in any treatment group",
             "AR Table 70: summary table, no cutoff stated", "AR Table 26: summary, no cutoff",
             "AR Table 35: >=5% in either group (active control)", "AR Table 30: at least 5% (active control)"))

# Number of side effects labelled very common or common in section 4.8 of each pilot label
# (structure count only, made without opening any assessment report).
pilot_rows_per_product <- c(Hexacima = 12, Onerji = 31, Drovelis = 12, Stelara = 16, Beyonttra = 2,
                            Pyrukynd = 6, Peyona = 4, Leqembi = 14, Xgeva = 9, Ruconest = 1)
