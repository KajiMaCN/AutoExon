library(data.table)
library(dplyr)
library(ggplot2)
library(jsonlite)
library(magrittr)
library(PepMapViz)


# Require: sample_id, class, timepoint, batch(?)
# Protein.Group, Genes, Precursor.Charge, Precursor.Normalised
file <- '/workspace/datasets/Match_result-X401SC24071912_Z01_F001_B1_43/report-matched-M0.csv'
report <- read.csv(
  file, header = TRUE,
  stringsAsFactors = FALSE, check.names = FALSE
)
data <- as.data.table(report)

# Check data
head(data)
dim(data)
length(unique(data$Polypeptide.Novogene.ID))

data$label[
  data$label %in% c("maintain", "early_remit", "late_remit", "relapse", "remit")
] <- "Non-convert"
data$label[data$label == "convert"] <- "Convert"

# Exon information


file <- '/workspace/results/alb-1.json'
itih4_1 <- fromJSON(file)
print(itih4_1)


# Skipped functions: strip_sequence_DIANN, obtain_mod_DIANN 
# No need to strip sequence as we have that already
# No need to obtain PTM because we do not want to plot

gene_symbol <- "ITIH4"
print(gene_symbol %in% unique(data$Genes))

# Filter data for just peptides mapping to gene of interest
data_fltr <- data[data$Genes == gene_symbol, ]
head(data_fltr)
dim(data_fltr)

whole_seq <- data.frame(
  Region_Sequence = c(itih4_1$sequence),
  Genes = c(gene_symbol)
)
print(whole_seq)

# Matches all peptide sequences to the reference sequence provided 
matching_result <- match_and_calculate_positions(
  data_fltr,
  column = "Stripped.Sequence", # name of column with sequence
  whole_seq,
  match_columns = c("Genes"), # match column from data with whole_seq
  column_keep = c("Polypeptide.Novogene.ID", "label") # information to keep 
)
# colnames(matching_result)[1] <- "Sequence" # To avoid error
dim(matching_result)

# TODO: Raise issue - Allow matching_result to not have reps and PTM_position
matching_result$reps <- NA
matching_result$PTM_position <- NA
data_with_quantification <- peptide_quantification(
  whole_seq,
  matching_result,
  matching_columns = c("Genes"), # identifying columns for sequences in whole_seq
  distinct_columns = c("Polypeptide.Novogene.ID", "label"), # not aggregated for these conditions
  quantify_method = "PSM", # PSM or Area
  # area_column = "Precursor.Normalised",
  with_PTM = FALSE,
  reps = FALSE
)
dim(data_with_quantification)

# file <- 'tmp/astral/data_with_quantification.csv'
# write.csv(data_with_quantification, file, row.names = TRUE)

# Average PSMs across samples from the same label
data_avg <- data_with_quantification %>%
  select(-Genes, -Polypeptide.Novogene.ID) %>%
  group_by(Character, Position, label) %>%
  summarise(PSM = mean(PSM, na.rm = TRUE), .groups = "keep") %>%
  arrange(label, Position) %>%
  as.data.table()

# Domain information
domain <- itih4_1$exons
domain_color <- c(
  rep("#F8766D", 14),
  rep("#B79F00", 2),
  # rep("#00BA38", 1),
  # rep("#B79F00", 1),
  rep("#F8766D", 8)
  # rep("#00BFC4", 1),
  # rep("#619CFF", 1),
  # rep("#F564E3", 15)
)
domain$domain_color <- domain_color
head(domain)

# For samples grouped by label. It will plot the last sample label by means
# of geom_raster, and does not aggregate across samples.
p <- create_peptide_plot(
  data_avg,
  y_axis_vars = c("label"),
  x_axis_vars = NULL, # Facets according to region.
  y_expand = c(0.8, 0.8),
  x_expand = c(0, 0),
  theme_options = list(),
  labs_options = list(title = "ITIH4-1", x = "Position"),
  color_fill_column = 'PSM',
  fill_gradient_options = list(),
  label_size = 3,
  add_domain = TRUE, # TODO: Even if FALSE, domain will be plotted
  domain = domain, # TODO: ERROR if domain not provided
  domain_start_column = "start",
  domain_end_column = "end",
  domain_type_column = "exon",
  domain_fill_color_column = "domain_color",
  domain_label_size = 2.5, 
  # PTM = FALSE,
  # PTM_type_column = "PTM_type",
  # PTM_color = PTM_color,
  add_label = FALSE,
  # label_column = "Character",
  # column_order = list(Region_1 = 'VH, VL')
)
file = "/workspace/results/peptides-ITIH4-1-PSM.pdf"
ggsave(p, filename = file, width = 12, height = 2.5)


# # Count number of labels
# data %>%
#   select(Polypeptide.Novogene.ID, label) %>%
#   unique() %>%
#   group_by(label) %>%
#   summarise(count = n())
