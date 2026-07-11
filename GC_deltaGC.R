# ====================
# GC AND ∆GC ANALYSIS
# ====================

# ----------------
# Settings/inputs
# ----------------
metadata_xlsx <- "/your/file/path/sample_phages.xlsx"

out_dir <- "/your/directory/path/Results/GC_deltaGC"

# Folder where downloaded CDS FASTA files will be cached
cds_cache_dir <- file.path(out_dir, "downloaded_CDS_fastas")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cds_cache_dir, recursive = TRUE, showWarnings = FALSE)

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Optional: set your preferred genus order
genus_order <- c(
  "Custom",
  "Host",
  "Order",
  "Optional"
)

# Optional colors
set_colors <- c(
  "Host" = "#F8766D",
  "Phage" = "#9F79EE"
)

trna_presence_colors <- c(
  "Without tRNAs" = "#999999",
  "With tRNAs" = "#66CDAA"
)

lifestyle_colors <- c(
  "virulent" = "#EE30A7",
  "temperate" = "#00B0F6"
)

# =========
# PACKAGES
# =========
req_pkgs <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "readr",
  "seqinr",
  "rentrez",
  "ggplot2",
  "scales",
  "broom"
)

not_installed <- req_pkgs[!req_pkgs %in% installed.packages()[, "Package"]]

if (length(not_installed) > 0) {
  install.packages(not_installed, dependencies = TRUE)
}

invisible(lapply(req_pkgs, library, character.only = TRUE))

# ======================================================
# Helper functions for accessions and statistical tests
# ======================================================
standardize_text <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish()
}

clean_accession <- function(x) {
  x <- x %>%
    as.character() %>%
    stringr::str_squish() %>%
    na_if("") %>%
    na_if("NA")
  # Fix accidental spaces in RefSeq accessions, e.g. "NC 004586" -> "NC_004586"
  x <- stringr::str_replace(x, "^(NC|NZ|NM|NR|XM|XR|YP|XP)\\s+([0-9]+)", "\\1_\\2")
  x
}

clean_name_for_file <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9_.-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

stop_if_missing_cols <- function(df, cols, df_name = "metadata") {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      df_name,
      " is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

safe_wilcox <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }
  tryCatch(
    wilcox.test(x, y)$p.value,
    error = function(e) NA_real_
  )
}

safe_spearman <- function(x, y) {
  x <- x[!is.na(x) & !is.na(y)]
  y <- y[!is.na(x) & !is.na(y)]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(
      estimate = NA_real_,
      p.value = NA_real_,
      n = length(x)
    ))
  }
  out <- suppressWarnings(cor.test(x, y, method = "spearman"))
  tibble(
    estimate = unname(out$estimate),
    p.value = out$p.value,
    n = length(x)
  )
}

# ======================================
# Function to fetch CDS FASTA from NCBI
# ======================================
fetch_cds_fasta <- function(accession, cache_dir = cds_cache_dir) {
  
  accession <- clean_accession(accession)
  
  if (is.na(accession)) {
    return(NA_character_)
  }
  
  safe_acc <- clean_name_for_file(accession)
  out_file <- file.path(cache_dir, paste0(safe_acc, "_CDS.fasta"))
  
  if (file.exists(out_file) && file.info(out_file)$size > 0) {
    first_line <- tryCatch(
      readLines(out_file, n = 1, warn = FALSE),
      error = function(e) ""
    )
    
    if (length(first_line) > 0 && stringr::str_starts(first_line[1], ">")) {
      return(out_file)
    } else {
      warning("Cached file is not valid FASTA, deleting: ", out_file)
      file.remove(out_file)
    }
  }
  message("Fetching CDS FASTA from NCBI: ", accession)
  # --------------------------
  # Direct fetch by accession
  # --------------------------
  fasta_txt <- tryCatch(
    rentrez::entrez_fetch(
      db = "nuccore",
      id = accession,
      rettype = "fasta_cds_na",
      retmode = "text"
    ),
    error = function(e) {
      NA_character_
    }
  )
  # -------------------------------------------------
  # Resolve accession to NCBI UID, then fetch by UID
  # -------------------------------------------------
  if (is.na(fasta_txt) || !nzchar(fasta_txt) ||
      !stringr::str_starts(strsplit(fasta_txt, "\n")[[1]][1], ">")) {
    
    message("Direct fetch failed; trying Entrez search for: ", accession)
    
    search_res <- tryCatch(
      rentrez::entrez_search(
        db = "nuccore",
        term = paste0(accession, "[Accession]"),
        retmax = 5
      ),
      error = function(e) NULL
    )
    
    if (!is.null(search_res) && length(search_res$ids) > 0) {
      uid <- search_res$ids[[1]]
      message("Resolved ", accession, " to NCBI UID: ", uid)
      fasta_txt <- tryCatch(
        rentrez::entrez_fetch(
          db = "nuccore",
          id = uid,
          rettype = "fasta_cds_na",
          retmode = "text"
        ),
        error = function(e) {
          warning("Could not fetch CDS FASTA for accession/UID: ", accession, " / ", uid)
          NA_character_
        }
      )
    }
  }
  # -------------------------
  # Validate returned FASTA
  # -------------------------
  if (is.na(fasta_txt) || !nzchar(fasta_txt)) {
    warning("Empty CDS FASTA returned for accession: ", accession)
    return(NA_character_)
  }
  first_line <- strsplit(fasta_txt, "\n")[[1]][1]
  if (!stringr::str_starts(first_line, ">")) {
    warning("NCBI did not return valid FASTA for accession: ", accession)
    return(NA_character_)
  }
  writeLines(fasta_txt, out_file)
  out_file
}

# =========================
# GC calculation functions
# =========================
calculate_gc_content <- function(dna_seq) {
  c(
    Total_GC = seqinr::GC(dna_seq),
    GC1      = seqinr::GC1(dna_seq),
    GC2      = seqinr::GC2(dna_seq),
    GC3      = seqinr::GC3(dna_seq)
  )
}

calc_gc_from_cds_fasta <- function(fasta_file, sample_id, accession, sample_type) {
  
  if (is.na(fasta_file) || !file.exists(fasta_file)) {
    return(tibble())
  }
  
  x <- tryCatch(
    seqinr::read.fasta(file = fasta_file, seqtype = "DNA"),
    error = function(e) {
      warning("Could not read FASTA file: ", fasta_file)
      return(NULL)
    }
  )
  
  if (is.null(x) || length(x) == 0) {
    return(tibble())
  }
  
  gc <- lapply(x, calculate_gc_content)
  
  df <- as.data.frame(do.call(rbind, gc)) %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(
      sample_id = sample_id,
      accession = accession,
      sample_type = sample_type
    ) %>%
    select(
      Gene,
      sample_id,
      accession,
      sample_type,
      Total_GC,
      GC1,
      GC2,
      GC3
    )
  df
}

gc_means <- function(df) {
  df %>%
    summarise(
      Total_GC = mean(Total_GC, na.rm = TRUE),
      GC1      = mean(GC1, na.rm = TRUE),
      GC2      = mean(GC2, na.rm = TRUE),
      GC3      = mean(GC3, na.rm = TRUE),
      n_genes  = n(),
      .groups = "drop"
    )
}

# ----------------------------------------------
# Read metadata from multi-sheet Excel workbook
# ----------------------------------------------
sheet_names <- readxl::excel_sheets(metadata_xlsx)

meta_raw <- purrr::map_dfr(
  sheet_names,
  function(sh) {
    readxl::read_excel(metadata_xlsx, sheet = sh) %>%
      mutate(source_sheet = sh)
  }
)

# columns based on INPHARED database columns
stop_if_missing_cols(
  meta_raw,
  c(
    "Accession",
    "Phage_ID",
    "Lifestyle",
    "tRNAs",
    "Host_Genus",
    "Host_Species",
    "Host_Accession"
  ),
  df_name = "metadata"
)

# this script DOES allow for organisms with multiple chromosomes and has been run successfully on Vibrio species
# Add Host_Accession2 as NA if that column is absent from sheets with only one chromosome/accession
if (!("Host_Accession2" %in% names(meta_raw))) {
  meta_raw$Host_Accession2 <- NA_character_
}

meta <- meta_raw %>%
  mutate(
    across(
      .cols = where(is.character) &
        !any_of(c("Accession", "Host_Accession", "Host_Accession2")),
      .fns = standardize_text
    ),
    phage_accession = clean_accession(Accession),
    host_accession  = clean_accession(Host_Accession),
    host_accession2 = clean_accession(Host_Accession2),
    phage_id = as.character(Phage_ID),
    host_genus = as.character(Host_Genus),
    host_species = as.character(Host_Species),
    host_key = dplyr::case_when(
      !is.na(host_accession2) ~ paste(host_accession, host_accession2, sep = " + "),
      TRUE ~ host_accession
    ),
    lifestyle = tolower(as.character(Lifestyle)),
    lifestyle = case_when(
      lifestyle %in% c("temperate", "temp") ~ "temperate",
      lifestyle %in% c("virulent", "vir") ~ "virulent",
      TRUE ~ lifestyle
    ),
    tRNAs = as.numeric(tRNAs),
    has_tRNAs = case_when(
      is.na(tRNAs) ~ NA_character_,
      tRNAs > 0 ~ "With tRNAs",
      tRNAs == 0 ~ "Without tRNAs"
    ),
    has_tRNAs = factor(
      has_tRNAs,
      levels = c("Without tRNAs", "With tRNAs")
    ),
    host_genus = factor(host_genus, levels = genus_order)
  ) %>%
  filter(
    !is.na(phage_accession),
    !is.na(host_accession),
    !is.na(host_genus)
  )

meta <- meta %>%
  mutate(host_genus = factor(as.character(host_genus), levels = genus_order_final))

readr::write_csv(meta, file.path(out_dir, "metadata_cleaned_accession_based.csv"))

# ----------------------------------------
# Build organism table for GC calculation
# ----------------------------------------
phage_orgs <- meta %>%
  transmute(
    sample_id = phage_id,
    host_key = host_key,
    accession = phage_accession,
    sample_type = "phage",
    chromosome = NA_character_
  )

host_orgs_chr1 <- meta %>%
  transmute(
    sample_id = host_accession,
    host_key = host_key,
    accession = host_accession,
    sample_type = "host",
    chromosome = "chromosome_1"
  )

host_orgs_chr2 <- meta %>%
  filter(!is.na(host_accession2)) %>%
  transmute(
    sample_id = host_accession2,
    host_key = host_key,
    accession = host_accession2,
    sample_type = "host",
    chromosome = "chromosome_2"
  )

host_orgs <- bind_rows(host_orgs_chr1, host_orgs_chr2)

organisms <- bind_rows(phage_orgs, host_orgs) %>%
  distinct(accession, sample_type, .keep_all = TRUE) %>%
  filter(!is.na(accession))

# ----------------------------------------------
# Fetch FASTAs and calculate per-CDS GC content
# ----------------------------------------------
gc_out_dir <- file.path(out_dir, "GC")
dir.create(gc_out_dir, recursive = TRUE, showWarnings = FALSE)

per_gene_gc <- organisms %>%
  mutate(
    cds_fasta = purrr::map_chr(accession, fetch_cds_fasta)
  ) %>%
  mutate(
    gc_data = purrr::pmap(
      list(cds_fasta, sample_id, accession, sample_type),
      calc_gc_from_cds_fasta
    )
  ) %>%
  select(host_key, chromosome, gc_data) %>%
  tidyr::unnest(gc_data)

readr::write_csv(
  per_gene_gc,
  file.path(gc_out_dir, "ALL_organisms_per_CDS_GC.csv")
)

per_sample_mean_gc <- per_gene_gc %>%
  group_by(sample_id, accession, sample_type) %>%
  gc_means() %>%
  ungroup()

readr::write_csv(
  per_sample_mean_gc,
  file.path(gc_out_dir, "ALL_organisms_mean_CDS_GC.csv")
)

# ------------------------------------------------------
# Attach the GC data for phage and host to the metadata
# ------------------------------------------------------
gc_types <- c("GC1", "GC2", "GC3", "Total_GC")

phage_gc <- per_gene_gc %>%
  dplyr::filter(sample_type == "phage") %>%
  dplyr::left_join(
    meta %>%
      dplyr::select(
        phage_id,
        phage_accession,
        host_accession,
        host_accession2,
        host_genus,
        host_species,
        lifestyle,
        tRNAs,
        has_tRNAs,
        Family
      ),
    by = c(
      "sample_id" = "phage_id",
      "accession" = "phage_accession"
    )
  )

host_gc <- per_gene_gc %>%
  filter(sample_type == "host")

readr::write_csv(
  phage_gc,
  file.path(gc_out_dir, "ALL_phages_per_CDS_GC_with_metadata.csv")
)

readr::write_csv(
  host_gc,
  file.path(gc_out_dir, "ALL_hosts_per_CDS_GC.csv")
)

# --------------------------------------------------
# Calculate ∆GC (GC phage - GC host) for all phages
# --------------------------------------------------
delta_gc <- phage_gc %>%
  dplyr::filter(!is.na(host_key)) %>%
  dplyr::group_by(
    sample_id,
    host_key,
    host_accession,
    host_accession2,
    host_genus,
    host_species,
    lifestyle,
    tRNAs,
    has_tRNAs
  ) %>%
  dplyr::group_modify(~{
    phage_df <- .x
    keys <- .y
    matched_host <- host_gc %>%
      dplyr::filter(host_key == keys$host_key[[1]])
    
    if (nrow(matched_host) == 0) {
      return(tibble(
        GC_Type = gc_types,
        Host_Mean_GC = NA_real_,
        Phage_Mean_GC = sapply(gc_types, function(gt) mean(phage_df[[gt]], na.rm = TRUE)),
        Delta_GC = NA_real_,
        p_value = NA_real_,
        n_host_genes = 0L,
        n_phage_genes = nrow(phage_df)
      ))
    }
    purrr::map_dfr(gc_types, function(gt) {
      host_vals <- matched_host[[gt]]
      phage_vals <- phage_df[[gt]]
      host_mean <- mean(host_vals, na.rm = TRUE)
      phage_mean <- mean(phage_vals, na.rm = TRUE)
      tibble(
        GC_Type = gt,
        Host_Mean_GC = host_mean,
        Phage_Mean_GC = phage_mean,
        Delta_GC = phage_mean - host_mean,
        p_value = safe_wilcox(phage_vals, host_vals),
        n_host_genes = sum(!is.na(host_vals)),
        n_phage_genes = sum(!is.na(phage_vals))
      )
    })
  }) %>%
  ungroup() %>%
  rename(Phage = sample_id) %>%
  mutate(
    host_genus = factor(as.character(host_genus), levels = genus_order_final),
    GC_Type = factor(GC_Type, levels = c("GC1", "GC2", "GC3", "Total_GC")),
    lifestyle = factor(lifestyle, levels = c("virulent", "temperate")),
    has_tRNAs = factor(has_tRNAs, levels = c("Without tRNAs", "With tRNAs")),
    p_adj_BH = p.adjust(p_value, method = "BH")
  )

readr::write_csv(
  delta_gc,
  file.path(out_dir, "deltaGC_per_phage_vs_matched_host.csv")
)

# ================================
# Statistical summaries and tests
# ================================

# ----------
# Summaries
# ----------
delta_gc_genus_summary <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(host_genus)) %>%
  group_by(host_genus, GC_Type) %>%
  summarise(
    n_phages = n_distinct(Phage),
    mean_deltaGC = mean(Delta_GC, na.rm = TRUE),
    median_deltaGC = median(Delta_GC, na.rm = TRUE),
    sd_deltaGC = sd(Delta_GC, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  delta_gc_genus_summary,
  file.path(out_dir, "deltaGC_summary_by_genus.csv")
)

# lifestyle summary
delta_gc_lifestyle_summary <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(lifestyle), !is.na(host_genus)) %>%
  group_by(host_genus, lifestyle, GC_Type) %>%
  summarise(
    n_phages = n_distinct(Phage),
    mean_deltaGC = mean(Delta_GC, na.rm = TRUE),
    median_deltaGC = median(Delta_GC, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  delta_gc_lifestyle_summary,
  file.path(out_dir, "deltaGC_summary_lifestyle.csv")
)

# tRNA presence summary
delta_gc_trna_presence_summary <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(has_tRNAs), !is.na(host_genus)) %>%
  group_by(host_genus, has_tRNAs, GC_Type) %>%
  summarise(
    n_phages = n_distinct(Phage),
    mean_deltaGC = mean(Delta_GC, na.rm = TRUE),
    median_deltaGC = median(Delta_GC, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  delta_gc_trna_presence_summary,
  file.path(out_dir, "deltaGC_summary_tRNA_presence.csv")
)

# tRNA count summary
delta_gc_trna_numeric_summary <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(tRNAs), !is.na(host_genus)) %>%
  group_by(host_genus, GC_Type) %>%
  summarise(
    n_phages = n_distinct(Phage),
    spearman_rho = safe_spearman(tRNAs, Delta_GC)$estimate,
    spearman_p = safe_spearman(tRNAs, Delta_GC)$p.value,
    .groups = "drop"
  )

readr::write_csv(
  delta_gc_trna_numeric_summary,
  file.path(out_dir, "deltaGC_spearman_tRNA_number.csv")
)

# -----------------------------------------------
# Statistical comparisons of lifestyle and tRNAs
# -----------------------------------------------
# lifestyle comparison for all GC types
lifestyle_tests <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(lifestyle)) %>%
  group_by(host_genus, GC_Type) %>%
  summarise(
    n_virulent = n_distinct(Phage[lifestyle == "virulent"]),
    n_temperate = n_distinct(Phage[lifestyle == "temperate"]),
    p_value = safe_wilcox(
      Delta_GC[lifestyle == "virulent"],
      Delta_GC[lifestyle == "temperate"]
    ),
    mean_virulent = mean(Delta_GC[lifestyle == "virulent"], na.rm = TRUE),
    mean_temperate = mean(Delta_GC[lifestyle == "temperate"], na.rm = TRUE),
    median_virulent = median(Delta_GC[lifestyle == "virulent"], na.rm = TRUE),
    median_temperate = median(Delta_GC[lifestyle == "temperate"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  lifestyle_tests,
  file.path(out_dir, "deltaGC_lifestyle_wilcox.csv")
)

# tRNA presence/absence
trna_presence_tests <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(has_tRNAs)) %>%
  group_by(host_genus, GC_Type) %>%
  summarise(
    n_without = n_distinct(Phage[has_tRNAs == "Without tRNAs"]),
    n_with = n_distinct(Phage[has_tRNAs == "With tRNAs"]),
    p_value = safe_wilcox(
      Delta_GC[has_tRNAs == "Without tRNAs"],
      Delta_GC[has_tRNAs == "With tRNAs"]
    ),
    mean_without = mean(Delta_GC[has_tRNAs == "Without tRNAs"], na.rm = TRUE),
    mean_with = mean(Delta_GC[has_tRNAs == "With tRNAs"], na.rm = TRUE),
    median_without = median(Delta_GC[has_tRNAs == "Without tRNAs"], na.rm = TRUE),
    median_with = median(Delta_GC[has_tRNAs == "With tRNAs"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  trna_presence_tests,
  file.path(out_dir, "deltaGC_tRNA_presence_wilcox.csv")
)

# tRNA numeric counts
trna_numeric_tests <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(tRNAs)) %>%
  group_by(host_genus, GC_Type) %>%
  group_modify(~{
    safe_spearman(.x$tRNAs, .x$Delta_GC)
  }) %>%
  ungroup() %>%
  rename(
    spearman_rho = estimate,
    p_value = p.value
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  trna_numeric_tests,
  file.path(out_dir, "deltaGC_numeric_tRNA_spearman.csv")
)

# =======
# PLOTS
# =======
heat_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(size = 14, face = "bold")
  )

# ---------
# Heatmaps
# ---------
# lifestyle ∆GC
p_lifestyle_heat <- ggplot(
  delta_gc_lifestyle_summary,
  aes(x = lifestyle, y = host_genus, fill = mean_deltaGC)
) +
  geom_tile(color = "white") +
  facet_wrap(~ GC_Type, nrow = 1) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    name = "Mean ΔGC"
  ) +
  labs(
    x = "Lifestyle",
    y = "Host genus",
    title = "Mean ΔGC by lifestyle and host genus"
  ) +
  heat_theme


ggsave(
  file.path(plot_dir, "deltaGC_heatmap_lifestyle.png"),
  p_lifestyle_heat,
  width = 12,
  height = 5,
  dpi = 300
)

# tRNA presence/absence ∆GC
p_trna_presence_heat <- ggplot(
  delta_gc_trna_presence_summary,
  aes(x = has_tRNAs, y = host_genus, fill = mean_deltaGC)
) +
  geom_tile(color = "white") +
  facet_wrap(~ GC_Type, nrow = 1) +
  scale_fill_gradient2(
    low = "#3B82F6",
    mid = "white",
    high = "#EF4444",
    midpoint = 0,
    name = "Mean ΔGC"
  ) +
  labs(
    x = "tRNA presence",
    y = "Host genus",
    title = "Mean ΔGC by tRNA presence and host genus"
  ) +
  heat_theme


ggsave(
  file.path(plot_dir, "deltaGC_heatmap_tRNA_presence.png"),
  p_trna_presence_heat,
  width = 12,
  height = 5,
  dpi = 300
)

# ------------------------------
# Violin plot host vs. phage GC
# ------------------------------
# bind host and phage GC values and write to new file
host_phage_gc_long <- bind_rows(
  host_gc %>%
    left_join(
      meta %>%
        distinct(
          host_key,
          host_genus,
          host_species,
          host_accession,
          host_accession2
        ),
      by = "host_key"
    ) %>%
    mutate(
      comparison_group = "Host"
    ),
  phage_gc %>%
    mutate(
      comparison_group = "Phage"
    )
) %>%
  filter(!is.na(host_genus)) %>%
  select(
    Gene,
    sample_id,
    accession,
    sample_type,
    comparison_group,
    host_genus,
    GC1,
    GC2,
    GC3,
    Total_GC
  ) %>%
  pivot_longer(
    cols = all_of(gc_types),
    names_to = "GC_Type",
    values_to = "GC"
  ) %>%
  mutate(
    comparison_group = factor(comparison_group, levels = c("Host", "Phage")),
    host_genus = factor(as.character(host_genus), levels = genus_order_final),
    GC_Type = factor(GC_Type, levels = c("GC1", "GC2", "GC3", "Total_GC"))
  )

readr::write_csv(
  host_phage_gc_long,
  file.path(gc_out_dir, "ALL_hosts_phages_per_CDS_GC_long_for_plotting.csv")
)

# plot in a violin plot comparing phages with their hosts
p_host_phage_gc_violin <- host_phage_gc_long %>%
  filter(!is.na(GC), !is.na(comparison_group)) %>%
  ggplot(aes(x = comparison_group, y = GC, fill = comparison_group)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, color = "black") +
  facet_grid(GC_Type ~ host_genus) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = set_colors, drop = FALSE) +
  labs(
    x = "",
    y = "CDS GC percentage",
    title = "Host vs phage CDS GC percentage by host genus"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(plot_dir, "host_vs_phage_GC_violin.png"),
  p_host_phage_gc_violin,
  width = 18,
  height = 9,
  dpi = 300
)

# -------------------
# Other violin plots
# -------------------
# lifestyle ∆GC
p_lifestyle_violin <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(lifestyle)) %>%
  ggplot(aes(x = lifestyle, y = Delta_GC, fill = lifestyle)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, color = "black") +
  facet_grid(GC_Type ~ host_genus) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = lifestyle_colors, drop = FALSE) +
  labs(
    x = "",
    y = "ΔGC: phage − matched host",
    title = "∆GC by lifestyle"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(plot_dir, "deltaGC_violin_lifestyle.png"),
  p_lifestyle_violin,
  width = 18,
  height = 9,
  dpi = 300
)

# tRNA presence/absence ∆GC
p_trna_presence_violin <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(has_tRNAs)) %>%
  ggplot(aes(x = has_tRNAs, y = Delta_GC, fill = has_tRNAs)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, color = "black") +
  facet_grid(GC_Type ~ host_genus) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = trna_presence_colors, drop = FALSE) +
  labs(
    x = "",
    y = "ΔGC: phage − matched host",
    title = "∆GC by tRNA presence"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(plot_dir, "deltaGC_violin_tRNA_presence.png"),
  p_trna_presence_violin,
  width = 18,
  height = 9,
  dpi = 300
)

# -------------------------------------------------
# Scatterplot for tRNA number correlation with ∆GC
# -------------------------------------------------
p_trna_numeric_scatter <- delta_gc %>%
  filter(!is.na(Delta_GC), !is.na(tRNAs)) %>%
  ggplot(aes(x = tRNAs, y = Delta_GC)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  facet_grid(GC_Type ~ host_genus, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Number of phage-encoded tRNAs",
    y = "ΔGC: phage − matched host",
    title = "∆GC by numeric tRNA count"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave(
  file.path(plot_dir, "deltaGC_scatter_numeric_tRNA_count.png"),
  p_trna_numeric_scatter,
  width = 18,
  height = 9,
  dpi = 300
)

# ======
# DONE!
# ======
msg("All done! :) Outputs saved in: ", out_dir)