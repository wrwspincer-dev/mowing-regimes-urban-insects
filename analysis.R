# MSc Dissertation: Mowing Regimes and Urban Insect Communities
# Statistical analysis and figures


# 1. Packages ------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(glmmTMB)
library(DHARMa)
library(emmeans)


# 2. Import and prepare data ---------------------------------------------

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


## Standardise taxonomy --------------------------------------------------

taxa_raw <- taxa_raw %>%
  mutate(
    Family = if_else(
      Family == "Apionidae",
      "Brentidae",
      Family
    )
  )


## Prepare sample metadata ----------------------------------------------

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
    ),
    
    Transect_ID = interaction(
      Site_Code,
      Transect,
      drop = TRUE
    )
  )


# 3. Figure settings -----------------------------------------------------

site_colours <- c(
  "Aberdare" = "#0072B2",
  "Main Building" = "#56B4E9",
  "Bute Building" = "#006D5B",
  "Redwood" = "#009E73",
  "ASSL" = "#4DB6A3",
  "Temple of Peace" = "#80CDC1",
  "SPARK" = "#6A51A3",
  "CUBRIC" = "#9E9AC8"
)

regime_colours <- c(
  "Frequently Mown" = "#0072B2",
  "Low Mowing" = "#009E73",
  "No Mowing" = "#6A51A3"
)

# 4. Hypothesis 1: Total insect abundance --------------------------------


# Prepare response variable

insect_abundance <- taxa_raw %>%
  filter(Class == "Insecta") %>%
  group_by(Sample_ID) %>%
  summarise(
    Insect_Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

sample_data <- sample_data %>%
  left_join(insect_abundance, by = "Sample_ID") %>%
  mutate(Insect_Abundance = coalesce(Insect_Abundance, 0))


# Descriptive statistics

h1_regime_summary <- sample_data %>%
  group_by(Mowing_Regime) %>%
  summarise(
    n = n(),
    Mean = mean(Insect_Abundance),
    Median = median(Insect_Abundance),
    SD = sd(Insect_Abundance),
    Minimum = min(Insect_Abundance),
    Maximum = max(Insect_Abundance),
    .groups = "drop"
  )

h1_regime_summary

sample_data %>%
  summarise(
    Mean = mean(Insect_Abundance),
    Variance = var(Insect_Abundance),
    Variance_Mean_Ratio = Variance / Mean,
    Zero_Counts = sum(Insect_Abundance == 0)
  )


# Model fitting and diagnostics

h1_model_hierarchical <- glmmTMB(
  Insect_Abundance ~ Mowing_Regime + Month +
    (1 | Site) + (1 | Transect_ID),
  family = nbinom2,
  data = sample_data
)

h1_model <- glmmTMB(
  Insect_Abundance ~ Mowing_Regime + Month,
  family = nbinom2,
  data = sample_data
)

VarCorr(h1_model_hierarchical)
AIC(h1_model_hierarchical, h1_model)

set.seed(123)

h1_residuals <- simulateResiduals(
  h1_model,
  n = 1000
)

plot(h1_residuals)
testUniformity(h1_residuals)
testDispersion(h1_residuals)
testOutliers(h1_residuals)


# Model results

drop1(h1_model, test = "Chisq")

h1_emmeans <- emmeans(
  h1_model,
  ~ Mowing_Regime,
  type = "response"
)

h1_emmeans

summary(
  pairs(h1_emmeans, adjust = "tukey"),
  infer = c(TRUE, TRUE),
  type = "response"
)

h1_month_emmeans <- emmeans(
  h1_model,
  ~ Month,
  type = "response"
)

h1_month_emmeans

summary(
  pairs(h1_month_emmeans, adjust = "tukey"),
  infer = c(TRUE, TRUE),
  type = "response"
)


# Figure 1

h1_plot_estimates <- as.data.frame(h1_emmeans)

figure_1 <- ggplot(
  sample_data,
  aes(x = Mowing_Regime, y = Insect_Abundance)
) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.35,
    colour = "grey40"
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    fill = NA,
    colour = "grey40"
  ) +
  geom_point(
    aes(colour = Site),
    size = 2.5,
    alpha = 0.8,
    position = position_jitter(
      width = 0.15,
      height = 0,
      seed = 123
    )
  ) +
  geom_errorbar(
    data = h1_plot_estimates,
    aes(
      x = Mowing_Regime,
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 1
  ) +
  geom_point(
    data = h1_plot_estimates,
    aes(x = Mowing_Regime, y = response),
    inherit.aes = FALSE,
    shape = 18,
    size = 4.5
  ) +
  scale_colour_manual(values = site_colours) +
  labs(
    x = "Mowing Regime",
    y = "Total Insect Abundance",
    colour = "Site"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

figure_1

# 5. Hypothesis 2: Hemiptera and Coleoptera abundance --------------------


# Prepare abundance data

h2_abundance <- taxa_raw %>%
  filter(Order %in% c("Hemiptera", "Coleoptera")) %>%
  group_by(Sample_ID, Order) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Order,
    values_from = Abundance,
    values_fill = 0
  )

sample_data <- sample_data %>%
  left_join(h2_abundance, by = "Sample_ID") %>%
  mutate(
    Hemiptera = coalesce(Hemiptera, 0),
    Coleoptera = coalesce(Coleoptera, 0)
  )


# Descriptive statistics

h2_regime_summary <- sample_data %>%
  group_by(Mowing_Regime) %>%
  summarise(
    n = n(),
    Hemiptera_Mean = mean(Hemiptera),
    Hemiptera_Median = median(Hemiptera),
    Hemiptera_SD = sd(Hemiptera),
    Hemiptera_Minimum = min(Hemiptera),
    Hemiptera_Maximum = max(Hemiptera),
    Coleoptera_Mean = mean(Coleoptera),
    Coleoptera_Median = median(Coleoptera),
    Coleoptera_SD = sd(Coleoptera),
    Coleoptera_Minimum = min(Coleoptera),
    Coleoptera_Maximum = max(Coleoptera),
    .groups = "drop"
  )

h2_regime_summary

sample_data %>%
  summarise(
    Hemiptera_Mean = mean(Hemiptera),
    Hemiptera_Median = median(Hemiptera),
    Hemiptera_Variance = var(Hemiptera),
    Hemiptera_VMR = Hemiptera_Variance / Hemiptera_Mean,
    Hemiptera_Zeros = sum(Hemiptera == 0),
    Coleoptera_Mean = mean(Coleoptera),
    Coleoptera_Median = median(Coleoptera),
    Coleoptera_Variance = var(Coleoptera),
    Coleoptera_VMR = Coleoptera_Variance / Coleoptera_Mean,
    Coleoptera_Zeros = sum(Coleoptera == 0)
  )


# Order-specific models

h2_hemiptera_hierarchical <- glmmTMB(
  Hemiptera ~ Mowing_Regime + Month +
    (1 | Site) + (1 | Transect_ID),
  family = nbinom2,
  data = sample_data
)

h2_hemiptera_model <- glmmTMB(
  Hemiptera ~ Mowing_Regime + Month,
  family = nbinom2,
  data = sample_data
)

VarCorr(h2_hemiptera_hierarchical)
AIC(h2_hemiptera_hierarchical, h2_hemiptera_model)
drop1(h2_hemiptera_model, test = "Chisq")

h2_hemiptera_emmeans <- emmeans(
  h2_hemiptera_model,
  ~ Mowing_Regime,
  type = "response"
)

h2_hemiptera_emmeans

summary(
  pairs(h2_hemiptera_emmeans, adjust = "tukey"),
  infer = c(TRUE, TRUE),
  type = "response"
)

emmeans(
  h2_hemiptera_model,
  ~ Month,
  type = "response"
)


h2_coleoptera_hierarchical <- glmmTMB(
  Coleoptera ~ Mowing_Regime + Month +
    (1 | Site) + (1 | Transect_ID),
  family = nbinom2,
  data = sample_data
)

h2_coleoptera_model <- glmmTMB(
  Coleoptera ~ Mowing_Regime + Month,
  family = nbinom2,
  data = sample_data
)

VarCorr(h2_coleoptera_hierarchical)
AIC(h2_coleoptera_hierarchical, h2_coleoptera_model)
drop1(h2_coleoptera_model, test = "Chisq")

h2_coleoptera_emmeans <- emmeans(
  h2_coleoptera_model,
  ~ Mowing_Regime,
  type = "response"
)

h2_coleoptera_emmeans

summary(
  pairs(h2_coleoptera_emmeans, adjust = "tukey"),
  infer = c(TRUE, TRUE),
  type = "response"
)

emmeans(
  h2_coleoptera_model,
  ~ Month,
  type = "response"
)


# Compare mowing responses between orders

h2_long <- sample_data %>%
  select(
    Sample_ID,
    Site,
    Mowing_Regime,
    Month,
    Transect_ID,
    Hemiptera,
    Coleoptera
  ) %>%
  pivot_longer(
    cols = c(Hemiptera, Coleoptera),
    names_to = "Taxon",
    values_to = "Abundance"
  ) %>%
  mutate(
    Taxon = factor(
      Taxon,
      levels = c("Hemiptera", "Coleoptera")
    )
  )

h2_interaction_hierarchical <- glmmTMB(
  Abundance ~ Taxon * Mowing_Regime + Month +
    (1 | Site) + (1 | Transect_ID) + (1 | Sample_ID),
  family = nbinom2,
  dispformula = ~ Taxon,
  data = h2_long
)

h2_interaction_model <- glmmTMB(
  Abundance ~ Taxon * Mowing_Regime + Month +
    (1 | Sample_ID),
  family = nbinom2,
  dispformula = ~ Taxon,
  data = h2_long
)

VarCorr(h2_interaction_hierarchical)
AIC(h2_interaction_hierarchical, h2_interaction_model)

h2_no_interaction <- glmmTMB(
  Abundance ~ Taxon + Mowing_Regime + Month +
    (1 | Sample_ID),
  family = nbinom2,
  dispformula = ~ Taxon,
  data = h2_long
)

anova(h2_no_interaction, h2_interaction_model)

emmeans(
  h2_interaction_model,
  ~ Mowing_Regime | Taxon,
  type = "response"
)

h2_joint_emmeans <- emmeans(
  h2_interaction_model,
  ~ Taxon * Mowing_Regime
)

summary(
  contrast(
    h2_joint_emmeans,
    interaction = "pairwise"
  ),
  infer = c(TRUE, TRUE),
  type = "response"
)


# Model diagnostics

set.seed(123)

h2_hemiptera_residuals <- simulateResiduals(
  h2_hemiptera_model,
  n = 1000
)

plot(h2_hemiptera_residuals)
testUniformity(h2_hemiptera_residuals)
testDispersion(h2_hemiptera_residuals)
testOutliers(h2_hemiptera_residuals)

h2_coleoptera_residuals <- simulateResiduals(
  h2_coleoptera_model,
  n = 1000
)

plot(h2_coleoptera_residuals)
testUniformity(h2_coleoptera_residuals)
testDispersion(h2_coleoptera_residuals)
testOutliers(h2_coleoptera_residuals)

h2_interaction_residuals <- simulateResiduals(
  h2_interaction_model,
  n = 1000
)

plot(h2_interaction_residuals)
testUniformity(h2_interaction_residuals)
testDispersion(h2_interaction_residuals)
testOutliers(h2_interaction_residuals)


# Figure 2

h2_figure_estimates <- bind_rows(
  as.data.frame(h2_hemiptera_emmeans) %>%
    mutate(Taxon = "Hemiptera"),
  as.data.frame(h2_coleoptera_emmeans) %>%
    mutate(Taxon = "Coleoptera")
) %>%
  mutate(
    Taxon = factor(
      Taxon,
      levels = c("Hemiptera", "Coleoptera")
    )
  )

figure2_panel_labels <- data.frame(
  Taxon = factor(
    c("Hemiptera", "Coleoptera"),
    levels = c("Hemiptera", "Coleoptera")
  ),
  label = c("a", "b")
)

figure_2 <- ggplot(
  h2_long,
  aes(x = Mowing_Regime, y = Abundance)
) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.35,
    colour = "grey40"
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    fill = NA,
    colour = "grey40"
  ) +
  geom_point(
    aes(colour = Site),
    size = 2.5,
    alpha = 0.8,
    position = position_jitter(
      width = 0.15,
      height = 0,
      seed = 123
    )
  ) +
  geom_errorbar(
    data = h2_figure_estimates,
    aes(
      x = Mowing_Regime,
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 1
  ) +
  geom_point(
    data = h2_figure_estimates,
    aes(x = Mowing_Regime, y = response),
    inherit.aes = FALSE,
    shape = 18,
    size = 4.5
  ) +
  geom_text(
    data = figure2_panel_labels,
    aes(
      x = "Frequently Mown",
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1.8,
    vjust = 1.4,
    size = 5,
    fontface = "bold"
  ) +
  facet_wrap(
    ~ Taxon,
    scales = "free_y"
  ) +
  scale_colour_manual(values = site_colours) +
  labs(
    x = "Mowing Regime",
    y = "Abundance",
    colour = "Site"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 13),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

figure_2

# 6. Hypothesis 3: Coleoptera family richness ----------------------------


# Prepare richness data

h3_richness <- taxa_raw %>%
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

sample_data <- sample_data %>%
  left_join(h3_richness, by = "Sample_ID") %>%
  mutate(
    Coleoptera_Richness = coalesce(
      Coleoptera_Richness,
      0L
    )
  )


# Descriptive statistics

h3_regime_summary <- sample_data %>%
  group_by(Mowing_Regime) %>%
  summarise(
    n = n(),
    Mean = mean(Coleoptera_Richness),
    Median = median(Coleoptera_Richness),
    SD = sd(Coleoptera_Richness),
    Minimum = min(Coleoptera_Richness),
    Maximum = max(Coleoptera_Richness),
    Zero_Samples = sum(Coleoptera_Richness == 0),
    .groups = "drop"
  )

h3_regime_summary

sample_data %>%
  summarise(
    Mean = mean(Coleoptera_Richness),
    Median = median(Coleoptera_Richness),
    Variance = var(Coleoptera_Richness),
    Minimum = min(Coleoptera_Richness),
    Maximum = max(Coleoptera_Richness),
    Zero_Samples = sum(Coleoptera_Richness == 0)
  )


# Model fitting and diagnostics

h3_model_hierarchical <- glmmTMB(
  Coleoptera_Richness ~ Mowing_Regime + Month +
    (1 | Site) + (1 | Transect_ID),
  family = poisson,
  data = sample_data
)

h3_model <- glmmTMB(
  Coleoptera_Richness ~ Mowing_Regime + Month,
  family = poisson,
  data = sample_data
)

VarCorr(h3_model_hierarchical)
AIC(h3_model_hierarchical, h3_model)

set.seed(123)

h3_residuals <- simulateResiduals(
  h3_model,
  n = 1000
)

plot(h3_residuals)
testUniformity(h3_residuals)
testDispersion(h3_residuals)
testOutliers(h3_residuals)
testQuantiles(h3_residuals)


# Model results

drop1(h3_model, test = "Chisq")

h3_emmeans <- emmeans(
  h3_model,
  ~ Mowing_Regime,
  type = "response"
)

h3_emmeans

summary(
  pairs(h3_emmeans, adjust = "tukey"),
  infer = c(TRUE, TRUE),
  type = "response"
)

h3_month_emmeans <- emmeans(
  h3_model,
  ~ Month,
  type = "response"
)

h3_month_emmeans

summary(
  pairs(h3_month_emmeans, adjust = "tukey"),
  infer = c(TRUE, TRUE),
  type = "response"
)


# Figure 3

h3_plot_estimates <- as.data.frame(h3_emmeans)

figure_3 <- ggplot(
  sample_data,
  aes(
    x = Mowing_Regime,
    y = Coleoptera_Richness
  )
) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.35,
    colour = "grey40"
  ) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    fill = NA,
    colour = "grey40"
  ) +
  geom_point(
    aes(colour = Site),
    size = 2.5,
    alpha = 0.8,
    position = position_jitter(
      width = 0.15,
      height = 0,
      seed = 123
    )
  ) +
  geom_errorbar(
    data = h3_plot_estimates,
    aes(
      x = Mowing_Regime,
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    inherit.aes = FALSE,
    width = 0.12,
    linewidth = 1
  ) +
  geom_point(
    data = h3_plot_estimates,
    aes(
      x = Mowing_Regime,
      y = rate
    ),
    inherit.aes = FALSE,
    shape = 18,
    size = 4.5
  ) +
  scale_colour_manual(values = site_colours) +
  scale_y_continuous(
    breaks = 0:7,
    limits = c(0, 7.3),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Mowing Regime",
    y = "Coleoptera Family Richness",
    colour = "Site"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 13),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

figure_3

# 7. Hypothesis 4: Coleoptera family composition -------------------------


# Build sample-level community matrix

h4_family_long <- taxa_raw %>%
  filter(
    Order == "Coleoptera",
    !is.na(Family),
    Family != ""
  ) %>%
  group_by(Sample_ID, Family) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

h4_comm_full <- h4_family_long %>%
  complete(
    Sample_ID = sample_data$Sample_ID,
    Family,
    fill = list(Abundance = 0)
  ) %>%
  pivot_wider(
    names_from = Family,
    values_from = Abundance,
    values_fill = 0
  )

h4_comm_nonempty <- h4_comm_full %>%
  filter(
    rowSums(across(-Sample_ID)) > 0
  )

h4_sample_ids <- h4_comm_nonempty$Sample_ID

h4_comm_matrix <- h4_comm_nonempty %>%
  select(-Sample_ID) %>%
  as.data.frame()

rownames(h4_comm_matrix) <- h4_sample_ids

h4_metadata <- sample_data %>%
  filter(Sample_ID %in% h4_sample_ids) %>%
  slice(match(h4_sample_ids, Sample_ID))

dim(h4_comm_matrix)
sum(rowSums(h4_comm_matrix) == 0)
all(rownames(h4_comm_matrix) == h4_metadata$Sample_ID)


# NMDS ordination

h4_family_frequency <- data.frame(
  Family = colnames(h4_comm_matrix),
  Samples_Present = colSums(h4_comm_matrix > 0),
  Total_Abundance = colSums(h4_comm_matrix)
) %>%
  arrange(Samples_Present, Total_Abundance)

h4_family_frequency

set.seed(123)

h4_nmds_full <- metaMDS(
  h4_comm_matrix,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

h4_nmds_full$stress

# Singleton families were removed for ordination only because the
# full matrix produced a near-zero, poorly informative NMDS solution.

h4_singleton_families <- h4_family_frequency %>%
  filter(Samples_Present == 1) %>%
  pull(Family)

h4_comm_filtered <- h4_comm_matrix[
  ,
  !colnames(h4_comm_matrix) %in% h4_singleton_families,
  drop = FALSE
]

h4_comm_filtered <- h4_comm_filtered[
  rowSums(h4_comm_filtered) > 0,
  ,
  drop = FALSE
]

set.seed(123)

h4_nmds <- metaMDS(
  h4_comm_filtered,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

h4_nmds$stress

h4_nmds_scores <- as.data.frame(
  scores(h4_nmds, display = "sites")
)

h4_nmds_scores$Sample_ID <- rownames(h4_nmds_scores)

h4_nmds_scores <- h4_nmds_scores %>%
  left_join(
    sample_data %>%
      select(
        Sample_ID,
        Mowing_Regime,
        Site,
        Month
      ),
    by = "Sample_ID"
  )


# Mowing-regime effect

# Samples were pooled across months within fixed transects before
# testing mowing regime to avoid treating repeated samples as independent.

h4_transect_long <- h4_family_long %>%
  left_join(
    sample_data %>%
      select(
        Sample_ID,
        Transect_ID,
        Mowing_Regime
      ),
    by = "Sample_ID"
  ) %>%
  group_by(
    Transect_ID,
    Mowing_Regime,
    Family
  ) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

h4_transect_comm <- h4_transect_long %>%
  select(
    Transect_ID,
    Family,
    Abundance
  ) %>%
  pivot_wider(
    names_from = Family,
    values_from = Abundance,
    values_fill = 0
  )

h4_transect_metadata <- h4_transect_long %>%
  distinct(
    Transect_ID,
    Mowing_Regime
  ) %>%
  slice(
    match(
      h4_transect_comm$Transect_ID,
      Transect_ID
    )
  )

h4_transect_matrix <- h4_transect_comm %>%
  select(-Transect_ID) %>%
  as.data.frame()

rownames(h4_transect_matrix) <- h4_transect_comm$Transect_ID

all(
  rownames(h4_transect_matrix) ==
    h4_transect_metadata$Transect_ID
)

table(h4_transect_metadata$Mowing_Regime)

set.seed(123)

h4_permanova <- adonis2(
  h4_transect_matrix ~ Mowing_Regime,
  data = h4_transect_metadata,
  method = "bray",
  permutations = 9999
)

h4_permanova


# Multivariate dispersion

h4_bray <- vegdist(
  h4_transect_matrix,
  method = "bray"
)

h4_permdisp <- betadisper(
  h4_bray,
  group = h4_transect_metadata$Mowing_Regime,
  bias.adjust = TRUE
)

set.seed(123)

permutest(
  h4_permdisp,
  permutations = 9999
)

h4_dispersion_summary <- data.frame(
  Mowing_Regime = h4_transect_metadata$Mowing_Regime,
  Distance = h4_permdisp$distances
) %>%
  group_by(Mowing_Regime) %>%
  summarise(
    n = n(),
    Mean = mean(Distance),
    SD = sd(Distance),
    Median = median(Distance),
    .groups = "drop"
  )

h4_dispersion_summary
TukeyHSD(h4_permdisp)


# Pairwise mowing-regime PERMANOVAs

h4_regime_pairs <- combn(
  levels(h4_transect_metadata$Mowing_Regime),
  2,
  simplify = FALSE
)

h4_pairwise_permanova <- lapply(
  h4_regime_pairs,
  function(regime_pair) {
    
    keep <- h4_transect_metadata$Mowing_Regime %in% regime_pair
    
    pair_matrix <- h4_transect_matrix[
      keep,
      ,
      drop = FALSE
    ]
    
    pair_metadata <- droplevels(
      h4_transect_metadata[keep, ]
    )
    
    set.seed(123)
    
    result <- adonis2(
      pair_matrix ~ Mowing_Regime,
      data = pair_metadata,
      method = "bray",
      permutations = 9999
    )
    
    data.frame(
      Contrast = paste(regime_pair, collapse = " vs "),
      F = result$F[1],
      R2 = result$R2[1],
      P = result$`Pr(>F)`[1]
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    P_Adjusted = p.adjust(
      P,
      method = "holm"
    )
  )

h4_pairwise_permanova


# Temporal variation in family composition

set.seed(123)

h4_permanova_month <- adonis2(
  h4_comm_matrix ~ Mowing_Regime + Month,
  data = h4_metadata,
  method = "bray",
  permutations = 9999,
  strata = h4_metadata$Transect_ID,
  by = "terms"
)

h4_permanova_month

h4_month_pairs <- combn(
  levels(h4_metadata$Month),
  2,
  simplify = FALSE
)

h4_pairwise_month <- lapply(
  h4_month_pairs,
  function(month_pair) {
    
    pair_metadata <- h4_metadata %>%
      filter(Month %in% month_pair) %>%
      group_by(Transect_ID) %>%
      filter(n_distinct(Month) == 2) %>%
      ungroup() %>%
      droplevels()
    
    pair_matrix <- h4_comm_matrix[
      pair_metadata$Sample_ID,
      ,
      drop = FALSE
    ]
    
    set.seed(123)
    
    result <- adonis2(
      pair_matrix ~ Month,
      data = pair_metadata,
      method = "bray",
      permutations = 9999,
      strata = pair_metadata$Transect_ID
    )
    
    data.frame(
      Contrast = paste(month_pair, collapse = " vs "),
      Transects = n_distinct(pair_metadata$Transect_ID),
      F = result$F[1],
      R2 = result$R2[1],
      P = result$`Pr(>F)`[1]
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    P_Adjusted = p.adjust(
      P,
      method = "holm"
    )
  )

h4_pairwise_month


# Sensitivity analysis without singleton families

h4_transect_matrix_filtered <- h4_transect_matrix[
  ,
  !colnames(h4_transect_matrix) %in% h4_singleton_families,
  drop = FALSE
]

set.seed(123)

adonis2(
  h4_transect_matrix_filtered ~ Mowing_Regime,
  data = h4_transect_metadata,
  method = "bray",
  permutations = 9999
)

h4_bray_filtered <- vegdist(
  h4_transect_matrix_filtered,
  method = "bray"
)

h4_permdisp_filtered <- betadisper(
  h4_bray_filtered,
  group = h4_transect_metadata$Mowing_Regime,
  bias.adjust = TRUE
)

set.seed(123)

permutest(
  h4_permdisp_filtered,
  permutations = 9999
)

h4_pairwise_permanova_filtered <- lapply(
  h4_regime_pairs,
  function(regime_pair) {
    
    keep <- h4_transect_metadata$Mowing_Regime %in% regime_pair
    
    pair_matrix <- h4_transect_matrix_filtered[
      keep,
      ,
      drop = FALSE
    ]
    
    pair_metadata <- droplevels(
      h4_transect_metadata[keep, ]
    )
    
    set.seed(123)
    
    result <- adonis2(
      pair_matrix ~ Mowing_Regime,
      data = pair_metadata,
      method = "bray",
      permutations = 9999
    )
    
    data.frame(
      Contrast = paste(regime_pair, collapse = " vs "),
      F = result$F[1],
      R2 = result$R2[1],
      P = result$`Pr(>F)`[1]
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    P_Adjusted = p.adjust(
      P,
      method = "holm"
    )
  )

h4_pairwise_permanova_filtered


# Figure 4

figure_4 <- ggplot(
  h4_nmds_scores,
  aes(
    x = NMDS1,
    y = NMDS2,
    colour = Mowing_Regime
  )
) +
  stat_ellipse(
    aes(group = Mowing_Regime),
    level = 0.95,
    linewidth = 1,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3.5,
    alpha = 0.85
  ) +
  scale_colour_manual(
    values = regime_colours
  ) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    colour = "Mowing Regime"
  ) +
  coord_equal() +
  theme_classic() +
  theme(
    text = element_text(size = 13),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

figure_4

# 8. Site-level temporal variation ---------------------------------------


# Site-level model

h5_site_model <- glmmTMB(
  Insect_Abundance ~ Site * Month + (1 | Transect_ID),
  family = nbinom2,
  data = sample_data
)

h5_site_model_additive <- glmmTMB(
  Insect_Abundance ~ Site + Month + (1 | Transect_ID),
  family = nbinom2,
  data = sample_data
)

anova(
  h5_site_model_additive,
  h5_site_model
)

set.seed(123)

h5_residuals <- simulateResiduals(
  h5_site_model,
  n = 1000
)

plot(h5_residuals)
testUniformity(h5_residuals)
testDispersion(h5_residuals)
testOutliers(h5_residuals)


# Within-site month comparisons

h5_emmeans <- emmeans(
  h5_site_model,
  ~ Month | Site,
  type = "response"
)

h5_emmeans

h5_month_contrasts <- pairs(
  h5_emmeans,
  adjust = "tukey"
)

summary(
  h5_month_contrasts,
  infer = c(TRUE, TRUE),
  type = "response"
)


# Monthly site means

h5_site_month <- sample_data %>%
  group_by(
    Site,
    Month
  ) %>%
  summarise(
    Mean_Abundance = mean(Insect_Abundance),
    .groups = "drop"
  )


# Compact-letter groups from Tukey-adjusted contrasts

h5_significance <- h5_site_month %>%
  mutate(
    Significance = case_when(
      Site == "Aberdare" ~ "a",
      
      Site == "ASSL" & Month == "May" ~ "a",
      Site == "ASSL" & Month %in% c("June", "July") ~ "b",
      
      Site == "Bute Building" & Month %in% c("May", "June") ~ "a",
      Site == "Bute Building" & Month == "July" ~ "b",
      
      Site == "CUBRIC" & Month == "May" ~ "a",
      Site == "CUBRIC" & Month == "June" ~ "ab",
      Site == "CUBRIC" & Month == "July" ~ "b",
      
      Site == "Main Building" & Month == "May" ~ "a",
      Site == "Main Building" & Month %in% c("June", "July") ~ "b",
      
      Site == "Redwood" & Month == "May" ~ "b",
      Site == "Redwood" & Month == "June" ~ "a",
      Site == "Redwood" & Month == "July" ~ "b",
      
      Site == "SPARK" ~ "a",
      
      Site == "Temple of Peace" & Month == "May" ~ "a",
      Site == "Temple of Peace" & Month == "June" ~ "ab",
      Site == "Temple of Peace" & Month == "July" ~ "b"
    ),
    Label_Y = Mean_Abundance + 18
  )


# Figure 5

figure_5 <- ggplot() +
  geom_line(
    data = sample_data,
    aes(
      x = Month,
      y = Insect_Abundance,
      group = Transect_ID,
      colour = Site
    ),
    linewidth = 0.45,
    alpha = 0.25
  ) +
  geom_point(
    data = sample_data,
    aes(
      x = Month,
      y = Insect_Abundance,
      colour = Site
    ),
    size = 1.4,
    alpha = 0.35
  ) +
  geom_line(
    data = h5_site_month,
    aes(
      x = Month,
      y = Mean_Abundance,
      group = Site,
      colour = Site
    ),
    linewidth = 1
  ) +
  geom_point(
    data = h5_site_month,
    aes(
      x = Month,
      y = Mean_Abundance,
      colour = Site
    ),
    size = 2.6
  ) +
  geom_text(
    data = h5_significance,
    aes(
      x = Month,
      y = Label_Y,
      label = Significance
    ),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  facet_wrap(
    ~ Site,
    ncol = 4
  ) +
  scale_colour_manual(
    values = site_colours
  ) +
  labs(
    x = "Month",
    y = "Insect abundance"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    panel.background = element_rect(
      fill = "grey98",
      colour = "grey80",
      linewidth = 0.6
    ),
    panel.grid = element_blank(),
    panel.spacing = unit(
      0.8,
      "lines"
    ),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

figure_5

# 9. Reproducibility ------------------------------------------------------

sessionInfo()
