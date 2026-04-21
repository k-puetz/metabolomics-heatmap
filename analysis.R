# Load libraries
library(dplyr)
library(readxl)
library(ComplexHeatmap)
library(RColorBrewer)
library(tibble)

# Read data
metabolomics <- suppressWarnings(read_excel('metabolomics.xlsx', sheet = 2))
annotation_col <- read.csv( "sample_annotations.csv",
    stringsAsFactors = FALSE, check.names = FALSE)

# Set column names
metabolomics <- metabolomics %>%
  mutate(Name_unique = make.unique(as.character(Name))) %>%
  tibble::column_to_rownames("Name_unique")
colnames(annotation_col)[1] = "SampleID"

# Remove samples with high QC variance
metabolomics_filtered_qc_na <- metabolomics %>%
  filter(`Group CV [%]: QC` <= 30)

# Select columns with metabolomic measurements
area_cols <- grep("^Area_", colnames(metabolomics_filtered_qc_na), value = TRUE)
area_cols <- area_cols[!grepl("QC|Verd|Blank|zero", area_cols)]

# Create matrix, normalise data, log-transform and scale row wise
mat_qc_na <- metabolomics_filtered_qc_na[, area_cols]
mat_norm <- sweep(mat_qc_na, 2, apply(mat_qc_na, 2, median, na.rm = TRUE), FUN = "/")
mat_norm <- log10(mat_norm + 1)
mat_scaled_qc_na_norm <- t(scale(t(mat_norm)))

# Matching metadata to metabolomics
colnames(mat_scaled_qc_na_norm) <- sub(".*-", "", colnames(mat_scaled_qc_na_norm))
colnames(mat_scaled_qc_na_norm) <- sub("^0+", "", colnames(mat_scaled_qc_na_norm))

annotation_col <- annotation_col[colnames(mat_scaled_qc_na_norm), , drop = FALSE]

# Define colors

fam_colors <- c(
  "Akkermansiaceae"   = "#1b9e77",
  "Bacteroidaceae"    = "#54478c",
  "Bifidobacteriaceae"= "#7570b3",
  "Enterobacteriaceae"= "#bbdefb",
  "Enterococcaceae"   = "#ffee33",
  "Lachnospiraceae"   = "#e6ab02",
  "Prevotellaceae"    = "#a6761d",
  "Rhizobiaceae"      = "#66a61e",
  "Ruminococcaceae"   = "#1f78b4",
  "Tannerellaceae"    = "#b2df8a",
  "other"             = "#fb9a99"
)

ann_colors <- list(
  "shotgun" = c('0' = "green", '1' = "orange", '2' = "red"),
  "timepoint" = c(a = "blue", c = "purple", e = "yellow"),
  "entity" = c('myelom' = "blue", 'lymphom' = "purple"),
  "survival" = c('0' = "white", '1' = "black"),
  "InvSimpson" = c('Low'= "green", 'Medium' = "orange", 'High' = "red"),
  "Dominant Family" = fam_colors
)

# Plot

pdf("heatmap_metabolomics.pdf", width = 25, height = 20)

ht <- Heatmap(
  mat_scaled_qc_na_norm,
  name = "Expression",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_km = 4,
  show_row_names = FALSE,
  row_km = 4,
  top_annotation = HeatmapAnnotation(
  df = annotation_col,
  col = ann_colors
  )

)

draw(ht, padding = unit(c(2, 4, 2, 4), "cm"))
dev.off()
