# Supplementary taxonomic summaries --------------------------------------


# Packages

library(readxl)
library(dplyr)
library(tidyr)


# Import and prepare data

file_path <- file.path(
  "data",
  "Sample_Metadata.xlsx"
)

metadata_raw <- read_excel(
  file_path,
  sheet = "Sample_Metadata"
)

taxa_raw <- read_excel(
  file_path,
  sheet = "Identification_Sheet"
)

taxa_raw <- taxa_raw %>%
  mutate(
    Family = if_else(
      Family == "Apionidae",
      "Brentidae",
      Family
    )
  )

sample_data <- metadata_raw %>%
  mutate(
    Site = recode(
      Site,
      "Sparks" = "SPARK"
    ),
    
    Mowing_Regime = case_when(
      Site_Code %in% c("MB", "AB") ~ "Frequently Mown",
      Site_Code == "AS" & Transect == 1 ~ "Frequently Mown",
      Site_Code == "AS" & Transect == 2 ~ "Low Mowing",
      Site_Code %in% c("BB", "TP", "RW") ~ "Low Mowing",
      Site_Code %in% c("SP", "CU") ~ "No Mowing"
    ),
    
    Mowing_Regime = factor(
      Mowing_Regime,
      levels = c(
        "Frequently Mown",
        "Low Mowing",
        "No Mowing"
      )
    ),
    
    Month = factor(
      Month,
      levels = c(
        "May",
        "June",
        "July"
      )
    )
  )


# Hemiptera suborders

hemiptera_suborder_totals <- taxa_raw %>%
  filter(
    Order == "Hemiptera",
    Suborder %in% c(
      "Auchenorrhyncha",
      "Heteroptera",
      "Sternorrhyncha"
    )
  ) %>%
  group_by(Suborder) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

hemiptera_suborder_totals


hemiptera_suborder_regime <- taxa_raw %>%
  filter(
    Order == "Hemiptera",
    Suborder %in% c(
      "Auchenorrhyncha",
      "Heteroptera",
      "Sternorrhyncha"
    )
  ) %>%
  left_join(
    sample_data %>%
      select(Sample_ID, Mowing_Regime),
    by = "Sample_ID"
  ) %>%
  group_by(
    Suborder,
    Mowing_Regime
  ) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Mowing_Regime,
    values_from = Abundance,
    values_fill = 0
  )

hemiptera_suborder_regime


# Coleoptera families

coleoptera_family_regime <- taxa_raw %>%
  filter(
    Order == "Coleoptera",
    !is.na(Family),
    Family != ""
  ) %>%
  left_join(
    sample_data %>%
      select(Sample_ID, Mowing_Regime),
    by = "Sample_ID"
  ) %>%
  group_by(
    Family,
    Mowing_Regime
  ) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Mowing_Regime,
    values_from = Abundance,
    values_fill = 0
  )

coleoptera_family_regime


coleoptera_no_mowing_composition <- taxa_raw %>%
  filter(
    Order == "Coleoptera",
    !is.na(Family),
    Family != ""
  ) %>%
  left_join(
    sample_data %>%
      select(Sample_ID, Mowing_Regime),
    by = "Sample_ID"
  ) %>%
  filter(
    Mowing_Regime == "No Mowing"
  ) %>%
  group_by(Family) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Percentage = 100 * Abundance / sum(Abundance)
  ) %>%
  arrange(
    desc(Abundance)
  )

coleoptera_no_mowing_composition


# Coleoptera family richness by site

coleoptera_richness <- taxa_raw %>%
  filter(
    Order == "Coleoptera",
    !is.na(Family),
    Family != ""
  ) %>%
  group_by(Sample_ID) %>%
  summarise(
    Coleoptera_Richness = n_distinct(Family),
    .groups = "drop"
  )

richness_data <- sample_data %>%
  left_join(
    coleoptera_richness,
    by = "Sample_ID"
  ) %>%
  mutate(
    Coleoptera_Richness = coalesce(
      Coleoptera_Richness,
      0L
    )
  )


coleoptera_richness_site <- richness_data %>%
  group_by(
    Site,
    Mowing_Regime
  ) %>%
  summarise(
    n = n(),
    Mean = mean(Coleoptera_Richness),
    Median = median(Coleoptera_Richness),
    SD = sd(Coleoptera_Richness),
    Minimum = min(Coleoptera_Richness),
    Maximum = max(Coleoptera_Richness),
    .groups = "drop"
  ) %>%
  arrange(
    Mowing_Regime,
    desc(Mean)
  )

coleoptera_richness_site


coleoptera_richness_no_mowing <- richness_data %>%
  filter(
    Mowing_Regime == "No Mowing"
  ) %>%
  select(
    Sample_ID,
    Site,
    Month,
    Transect,
    Coleoptera_Richness
  ) %>%
  arrange(
    Site,
    Month,
    Transect
  )

coleoptera_richness_no_mowing
