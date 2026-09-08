# Appendix B. Taxonomic abundance data -----------------------------------


# Table 3: Major insect groups by sample

table_3 <- taxa_raw %>%
  filter(Class == "Insecta") %>%
  group_by(Sample_ID, Order) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Order,
    values_from = Abundance,
    values_fill = 0
  ) %>%
  mutate(
    Total_Insect_Abundance = rowSums(
      across(where(is.numeric))
    )
  ) %>%
  arrange(Sample_ID)

print(
  table_3,
  n = 48,
  width = Inf
)


# Table 4: Hemiptera abundance by major taxonomic group

hemiptera_groups <- c(
  "Auchenorrhyncha",
  "Heteroptera",
  "Sternorrhyncha",
  "Unresolved Hemiptera"
)

table_4_long <- taxa_raw %>%
  filter(Order == "Hemiptera") %>%
  mutate(
    Hemiptera_Group = case_when(
      Suborder == "Auchenorrhyncha" ~ "Auchenorrhyncha",
      Suborder == "Heteroptera" ~ "Heteroptera",
      Suborder == "Sternorrhyncha" ~ "Sternorrhyncha",
      TRUE ~ "Unresolved Hemiptera"
    )
  ) %>%
  group_by(
    Sample_ID,
    Hemiptera_Group
  ) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

table_4 <- expand_grid(
  Sample_ID = sample_data$Sample_ID,
  Hemiptera_Group = hemiptera_groups
) %>%
  left_join(
    table_4_long,
    by = c(
      "Sample_ID",
      "Hemiptera_Group"
    )
  ) %>%
  mutate(
    Abundance = replace_na(Abundance, 0)
  ) %>%
  pivot_wider(
    names_from = Hemiptera_Group,
    values_from = Abundance
  ) %>%
  left_join(
    sample_data %>%
      select(
        Sample_ID,
        Mowing_Regime
      ),
    by = "Sample_ID"
  ) %>%
  mutate(
    Total_Hemiptera =
      Auchenorrhyncha +
      Heteroptera +
      Sternorrhyncha +
      `Unresolved Hemiptera`
  ) %>%
  select(
    Sample_ID,
    Mowing_Regime,
    Auchenorrhyncha,
    Heteroptera,
    Sternorrhyncha,
    `Unresolved Hemiptera`,
    Total_Hemiptera
  ) %>%
  arrange(Sample_ID)

print(
  table_4,
  n = 48,
  width = Inf
)


# Table 5: Coleoptera family abundance by sample

coleoptera_families <- c(
  "Oedemeridae",
  "Coccinellidae",
  "Brentidae",
  "Curculionidae",
  "Chrysomelidae",
  "Cantharidae",
  "Staphylinidae",
  "Dermestidae",
  "Limnichidae",
  "Tenebrionidae",
  "Elateridae",
  "Melyridae",
  "Phalacridae",
  "Scarabaeidae"
)

table_5_long <- taxa_raw %>%
  filter(
    Order == "Coleoptera",
    !is.na(Family),
    Family != ""
  ) %>%
  group_by(
    Sample_ID,
    Family
  ) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

table_5_total <- taxa_raw %>%
  filter(Order == "Coleoptera") %>%
  group_by(Sample_ID) %>%
  summarise(
    Total_Coleoptera = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

table_5 <- expand_grid(
  Sample_ID = sample_data$Sample_ID,
  Family = coleoptera_families
) %>%
  left_join(
    table_5_long,
    by = c(
      "Sample_ID",
      "Family"
    )
  ) %>%
  mutate(
    Abundance = replace_na(Abundance, 0)
  ) %>%
  pivot_wider(
    names_from = Family,
    values_from = Abundance
  ) %>%
  mutate(
    Coleoptera_Family_Richness = rowSums(
      across(
        all_of(coleoptera_families),
        ~ .x > 0
      )
    ),
    Total_Family_Resolved_Coleoptera = rowSums(
      across(
        all_of(coleoptera_families)
      )
    )
  ) %>%
  left_join(
    sample_data %>%
      select(
        Sample_ID,
        Mowing_Regime
      ),
    by = "Sample_ID"
  ) %>%
  left_join(
    table_5_total,
    by = "Sample_ID"
  ) %>%
  mutate(
    Total_Coleoptera = replace_na(
      Total_Coleoptera,
      0
    )
  ) %>%
  select(
    Sample_ID,
    Mowing_Regime,
    all_of(coleoptera_families),
    Coleoptera_Family_Richness,
    Total_Family_Resolved_Coleoptera,
    Total_Coleoptera
  ) %>%
  arrange(Sample_ID)

print(
  table_5,
  n = 48,
  width = Inf
)