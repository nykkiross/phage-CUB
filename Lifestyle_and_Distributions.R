# =============================================================
# Lifestyle assignment and distribution of lifestyle and tRNAs
# =============================================================

# Use: when working with INPHARED data, lifestyle is not included
# this script assigns a phage as "putative virulent" or "putative temperate"
# distribution of lifestyle within a set of phages, as well as distribution
# of tRNAs is then calculated

# =========
# PACKAGES
# =========
req_pkgs <- c(
  "readxl",
  "tidyverse",
  "rentrez",
  "ggplot2"
)

not_installed <- req_pkgs[!req_pkgs %in% installed.packages()[, "Package"]]

if (length(not_installed) > 0) {
  install.packages(not_installed, dependencies = TRUE)
}

invisible(lapply(req_pkgs, library, character.only = TRUE))


# ========================================
# LIFESTYLE ASSIGNMENT terms and function
# ========================================
temperate_terms <- c(
  "integrase",
  "excisionase",
  "lysogen",
  "lysogenic"
)

classify_lifestyle_from_accession <- function(accession, sleep_time = 0.34) {
  
  Sys.sleep(sleep_time)
  
  gb_text <- tryCatch(
    entrez_fetch(
      db = "nuccore",
      id = accession,
      rettype = "gb",
      retmode = "text"
    ),
    error = function(e) NA_character_
  )
  
  if (is.na(gb_text)) {
    return(tibble(
      Accession = accession,
      lifestyle = "fetch_failed",
      matched_terms = NA_character_
    ))
  }
  
  gb_lower <- tolower(gb_text)
  
  found_terms <- temperate_terms[
    str_detect(gb_lower, temperate_terms)
  ]
  
  lifestyle <- case_when(
    length(found_terms) > 0 ~ "putative_temperate",
    TRUE ~ "putative_virulent"
  )
  
  tibble(
    Accession = accession,
    lifestyle = lifestyle,
    matched_terms = paste(found_terms, collapse = "; ")
  )
}

# ---------------------------------
# Load data and test with ONE host
# ---------------------------------
# metadata was formatted as in the INPHARED database for the generation of this script
# but should be able to run with any formatting so long as Accession # is included
metadata <- read_csv("/file/path/to/sample_phage_metadata.csv")

lifestyle_results <- metadata %>%
  distinct(Accession) %>%
  pull(Accession) %>%
  map_dfr(classify_lifestyle_from_accession)

metadata_with_lifestyle <- metadata %>%
  left_join(lifestyle_results, by = "Accession")

write_csv(metadata_with_lifestyle, "sample_phage_lifestyle.csv")

# -------------------------------------------
# Run lifestyle assignment on all phage sets
# -------------------------------------------
metadata_files <- list.files(
  "/file/path/to/all_metadata_files",
  pattern = "\\.csv$",
  full.names = TRUE
)

dir.create("lifestyle_outputs", showWarnings = FALSE)

for (fp in metadata_files) {
  
  message("Processing: ", basename(fp))
  
  metadata <- read_csv(fp, show_col_types = FALSE)
  
  lifestyle_results <- metadata %>%
    distinct(Accession) %>%
    pull(Accession) %>%
    map_dfr(classify_lifestyle_from_accession)
  
  output <- metadata %>%
    left_join(lifestyle_results, by = "Accession")
  
  out_name <- str_replace(
    basename(fp),
    "\\.csv$",
    "_lifestyle.csv"
  )
  
  write_csv(output, file.path("lifestyle_outputs", out_name))
}


# =====================================
# DISTRIBUTIONS for tRNA and lifestyle
# =====================================

# --------------------------------
# Merge all CSVs into one dataset
# --------------------------------
files <- list.files(
  "/file/path/to/all_metadata_files/lifestyle_outputs",
  pattern = "_lifestyle\\.csv$",
  full.names = TRUE
)

all_phages <- map_dfr(
  files,
  read_csv,
  show_col_types = FALSE
)

# make sure important columns are correct
all_phages <- all_phages %>%
  mutate(
    Host = as.factor(Host),
    lifestyle = as.factor(lifestyle),
    tRNAs = as.numeric(tRNAs)
  )

all_phages$Host <- factor(
  all_phages$Host,
  levels = c(
    "Set",
    "Custom",
    "Host",
    "Order",
    "If",
    "Desired"
  )
)

# ---------------------------------
# Summary statistics for lifestyle
# ---------------------------------
lifestyle_summary <- all_phages %>%
  count(Host, lifestyle) %>%
  group_by(Host) %>%
  mutate(
    proportion = n / sum(n)
  )

lifestyle_summary

write_csv(
  lifestyle_summary,
  "lifestyle_distribution_by_host.csv"
)

# lifestyle proportion plot
lifestyle_summary_plot <- lifestyle_summary %>%
  filter(lifestyle != "fetch_failed")

ggplot(lifestyle_summary_plot,
       aes(x = Host,
           y = proportion,
           fill = lifestyle)) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Lifestyle Distribution Across Host Genera",
    y = "Proportion of Phages"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# --------------------------------------------------------------------
# Chi-squared test - is lifestyle distribution related to host genus?
# --------------------------------------------------------------------

lifestyle_table <- table(
  all_phages$Host,
  all_phages$lifestyle
)

chisq.test(lifestyle_table)

# ------------------------------
# tRNA presence/absence summary
# ------------------------------
trna_presence_summary <- all_phages %>%
  group_by(Host) %>%
  summarise(
    total_phages = n(),
    phages_with_tRNAs = sum(tRNAs > 0, na.rm = TRUE),
    phages_without_tRNAs = sum(tRNAs == 0, na.rm = TRUE),
    proportion_with_tRNAs =
      phages_with_tRNAs / total_phages,
    percent_with_tRNAs =
      proportion_with_tRNAs * 100
  )

trna_presence_summary

write_csv(
  trna_presence_summary,
  "tRNA_presence_summary_by_host.csv"
)

# plot proportions
ggplot(
  trna_presence_summary,
  aes(x = Host,
      y = proportion_with_tRNAs,
      fill = Host)
) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Proportion of Phages Encoding tRNAs",
    y = "Proportion with tRNAs"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# stacked proportions plot
trna_presence_plot <- all_phages %>%
  mutate(
    tRNA_presence = ifelse(
      tRNAs > 0,
      "Has tRNAs",
      "No tRNAs"
    )
  ) %>%
  count(Host, tRNA_presence) %>%
  group_by(Host) %>%
  mutate(
    proportion = n / sum(n)
  )

ggplot(
  trna_presence_plot,
  aes(x = Host,
      y = proportion,
      fill = tRNA_presence)
) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Distribution of tRNA-Encoding Phages",
    y = "Proportion"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# -----------------------------------------------------------
# Chi-squared test - is tRNA presence related to host genus?
# -----------------------------------------------------------
trna_table <- table(
  all_phages$Host,
  all_phages$tRNAs > 0
)

chisq.test(trna_table)

# -----------------------------------
# Summary statistics for tRNA counts
# -----------------------------------
trna_summary <- all_phages %>%
  group_by(Host) %>%
  summarise(
    n_phages = n(),
    mean_tRNAs = mean(tRNAs, na.rm = TRUE),
    median_tRNAs = median(tRNAs, na.rm = TRUE),
    sd_tRNAs = sd(tRNAs, na.rm = TRUE),
    min_tRNAs = min(tRNAs, na.rm = TRUE),
    max_tRNAs = max(tRNAs, na.rm = TRUE)
  )

trna_summary

write_csv(
  trna_summary,
  "tRNA_distribution_summary_by_host.csv"
)

# plot it
ggplot(all_phages,
       aes(x = Host,
           y = tRNAs,
           fill = Host)) +
  geom_violin(trim = FALSE,
              alpha = 0.7) +
  geom_boxplot(width = 0.1,
               outlier.shape = NA) +
  theme_bw() +
  labs(
    title = "Distribution of tRNA Genes Across Phage Host Genera",
    y = "Number of tRNA Genes"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# -------------------------------------------------------
# Kruskal-Wallis test on tRNA counts across host genera
# -------------------------------------------------------
kruskal.test(tRNAs ~ Host,
             data = all_phages)

# optional: pairwise comparisons
pairwise.wilcox.test(
  all_phages$tRNAs,
  all_phages$Host,
  p.adjust.method = "BH"
)

# -------------------------------------------------
# tRNA distributions for tRNA-encoding phages ONLY
# -------------------------------------------------
trna_positive_phages <- all_phages %>%
  filter(tRNAs > 0)

trna_positive_summary <- trna_positive_phages %>%
  group_by(Host) %>%
  summarise(
    n_trna_phages = n(),
    mean_tRNAs = mean(tRNAs, na.rm = TRUE),
    median_tRNAs = median(tRNAs, na.rm = TRUE),
    sd_tRNAs = sd(tRNAs, na.rm = TRUE),
    min_tRNAs = min(tRNAs, na.rm = TRUE),
    max_tRNAs = max(tRNAs, na.rm = TRUE)
  )

trna_positive_summary

write_csv(
  trna_positive_summary,
  "tRNA_positive_summary_by_host.csv"
)

ggplot(
  trna_positive_phages,
  aes(x = Host,
      y = tRNAs,
      fill = Host)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(
    width = 0.1,
    outlier.shape = NA
  ) +
  theme_bw() +
  labs(
    title = "tRNA Counts Among tRNA-Encoding Phages",
    y = "Number of tRNA Genes"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

