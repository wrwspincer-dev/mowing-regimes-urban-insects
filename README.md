# Effects of mowing regimes on urban insect communities

Code and data supporting the MSc dissertation:

"Effects of mowing regimes on insect communities in urban grasslands"

MSc Global Ecology and Conservation
Cardiff University, 2026

## Study overview

This study examined insect communities across urban grasslands in Cardiff
managed under three mowing regimes: Frequently Mown, Low Mowing and No Mowing.

Sweep-net sampling was conducted along fixed transects during May, June and
July 2026. Analyses focused on total insect abundance, Hemiptera and Coleoptera
abundance, Coleoptera family richness, and Coleoptera family composition.

## Repository contents

- `analysis.R` - statistical analyses for the four study hypotheses and Figures 1-5
- `supplementary_summaries.R` - additional taxonomic summaries used in interpretation
- `appendix_tables.R` - code used to produce dissertation appendix tables
- `data/` - sample metadata and taxonomic identification data
- `figures/` - final figures produced from the analyses

## Software

Analyses were conducted in R 4.5.2.

Major packages used include:

- glmmTMB
- DHARMa
- emmeans
- vegan
- ggplot2
- dplyr
- tidyr
- readxl

Full session information can be reproduced by running the final section of
`analysis.R`.

## Running the analysis

Clone or download the repository and open the project directory in R.

The analysis assumes the input workbook is located at:

`data/Sample_Metadata.xlsx`

Run `analysis.R` from the repository root directory.

## Notes

Mowing regime was assigned at transect level. At ASSL, Transect 1 was
classified as Frequently Mown and Transect 2 as Low Mowing.

Coleoptera records assigned to Apionidae were standardised to Brentidae prior
to analysis.
