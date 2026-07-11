# =========================
# RSCU and ∆RSCU ANALYSIS
# =========================

# ----------------
# Settings/inputs
# ----------------
metadata_xlsx <- "/your/file/path/sample_phages.xlsx"

out_dir <- "/your/directory/path/Results/RSCU_deltaRSCU"

# Folder where downloaded CDS FASTA files will be cached
cds_cache_dir <- file.path(out_dir, "downloaded_CDS_fastas")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cds_cache_dir, recursive = TRUE, showWarnings = FALSE)

plot_dir <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

summary_dir <- file.path(out_dir, "summaries")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

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
  "readr",
  "seqinr",
  "stringr",
  "forcats",
  "ggplot2",
  "pheatmap",
  "broom",
  "vegan",
  "coin",
  "effsize",
  "FSA",
  "rentrez",
  "openxlsx"
)

not_installed <- req_pkgs[!req_pkgs %in% installed.packages()[, "Package"]]

if (length(not_installed) > 0) {
  install.packages(not_installed, dependencies = TRUE)
}

invisible(lapply(req_pkgs, library, character.only = TRUE))

theme_set(theme_minimal(base_size = 14))

# ==============
# RSCU Settings
# ==============
lifestyle_order <- c("virulent", "temperate")

trna_presence_order <- c("Without tRNAs", "With tRNAs")

include_stop_codons <- FALSE

rscu_method <- "by_counts"

sense_codons <- c(
  "ttt","ttc","tta","ttg","ctt","ctc","cta","ctg",
  "att","atc","ata","atg",
  "gtt","gtc","gta","gtg",
  "tct","tcc","tca","tcg","agt","agc",
  "cct","ccc","cca","ccg",
  "act","acc","aca","acg",
  "gct","gcc","gca","gcg",
  "tat","tac","cat","cac","caa","cag","aat","aac",
  "aaa","aag","gat","gac","gaa","gag",
  "tgt","tgc","tgg",
  "cgt","cgc","cga","cgg","aga","agg",
  "ggt","ggc","gga","ggg"
)

stop_codons <- c("taa","tag","tga")

# Optional drop Met/Trp for downstream ∆RSCU summaries
excluded_codons <- c("atg", "tgg")
filtered_codons <- setdiff(sense_codons, excluded_codons)

# ======================================================
# Helper functions for accessions and statistical tests
# ======================================================
msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste(...), "\n"))
}

standardize_text <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_squish()
}

clean_accession <- function(x) {
  x <- x %>%
    as.character() %>%
    stringr::str_squish() %>%
    dplyr::na_if("") %>%
    dplyr::na_if("NA")
  x <- stringr::str_replace(
    x,
    "^(NC|NZ|NM|NR|XM|XR|YP|XP)\\s+([0-9]+)",
    "\\1_\\2"
  )
  x <- stringr::str_remove(x, "\\.[0-9]+$")
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

safe_wilcox_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }
  wt <- try(stats::wilcox.test(x, y, exact = FALSE), silent = TRUE)
  if (inherits(wt, "try-error")) {
    return(NA_real_)
  }
  wt$p.value
}

safe_wilcox_tbl <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(tibble::tibble(
      W = NA_real_,
      p_value = NA_real_
    ))
  }
  wt <- try(stats::wilcox.test(x, y, exact = FALSE), silent = TRUE)
  if (inherits(wt, "try-error")) {
    return(tibble::tibble(
      W = NA_real_,
      p_value = NA_real_
    ))
  }
  tibble::tibble(
    W = unname(wt$statistic),
    p_value = wt$p.value
  )
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble::tibble(
      n = length(x),
      spearman_rho = NA_real_,
      p_value = NA_real_
    ))
  }
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman"))
  tibble::tibble(
    n = length(x),
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

cliffs_delta <- function(x, g) {
  g <- droplevels(factor(g))
  if (nlevels(g) != 2) {
    return(NA_real_)
  }
  x1 <- x[g == levels(g)[1]]
  x2 <- x[g == levels(g)[2]]
  x1 <- x1[is.finite(x1)]
  x2 <- x2[is.finite(x2)]
  if (length(x1) == 0 || length(x2) == 0) {
    return(NA_real_)
  }
  m <- outer(x1, x2, "-")
  (sum(m > 0) - sum(m < 0)) / (length(x1) * length(x2))
}

cap_title <- function(x) {
  tools::toTitleCase(gsub("_", " ", x))
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
  # ------------------------
  # Validate returned FASTA
  # ------------------------
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

# ===========================
# RSCU calculation functions
# ===========================
safe_read_fasta <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(list())
  }
  first_line <- tryCatch(readLines(path, n = 1, warn = FALSE), error = function(e) "")
  if (length(first_line) == 0 || !stringr::str_starts(first_line[1], ">")) {
    warning("Skipping non-FASTA file: ", path)
    return(list())
  }
  tryCatch(
    seqinr::read.fasta(file = path, as.string = FALSE, seqtype = "DNA"),
    error = function(e) {
      warning("Could not read FASTA file: ", path)
      return(list())
    }
  )
}

is_valid_cds <- function(x, min_length_nt = 200) {
  len <- length(x)
  len >= min_length_nt &&
    len %% 3 == 0 &&
    !any(!toupper(x) %in% c("A", "C", "G", "T"))
}

uco_to_tidy <- function(uco_named_vec, include_stop = FALSE) {
  v <- uco_named_vec
  nm <- tolower(names(v))
  names(v) <- nm
  full_codons <- if (include_stop) c(sense_codons, stop_codons) else sense_codons
  tibble::tibble(codon = full_codons) %>%
    dplyr::left_join(
      tibble::tibble(codon = nm, rscu = as.numeric(v)),
      by = "codon"
    ) %>%
    dplyr::mutate(rscu = dplyr::coalesce(rscu, 0))
}

compute_rscu_by_counts <- function(fasta_list) {
  if (length(fasta_list) == 0) {
    return(uco_to_tidy(rep(NA_real_, 64)))
  }
  nuc <- unlist(fasta_list, use.names = FALSE)
  nuc <- nuc[nuc %in% c("a", "t", "g", "c", "A", "T", "G", "C")]
  if (length(nuc) < 3) {
    return(uco_to_tidy(rep(NA_real_, 64)))
  }
  r <- seqinr::uco(nuc, index = "rscu")
  uco_to_tidy(r, include_stop = include_stop_codons)
}

compute_rscu_by_gene_mean <- function(fasta_list) {
  if (length(fasta_list) == 0) {
    return(uco_to_tidy(rep(NA_real_, 64)))
  }
  per_gene <- purrr::map(fasta_list, function(seq) {
    seq <- seq[seq %in% c("a", "t", "g", "c", "A", "T", "G", "C")]
    if (length(seq) < 3) {
      return(uco_to_tidy(rep(NA_real_, 64)))
    }
    r <- seqinr::uco(seq, index = "rscu")
    uco_to_tidy(r, include_stop = include_stop_codons)
  })
  dplyr::bind_rows(per_gene, .id = "gene") %>%
    dplyr::group_by(codon) %>%
    dplyr::summarise(
      rscu = mean(rscu, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_rscu_for_one_organism <- function(
    fasta_file,
    sample_id,
    accession,
    sample_type,
    host_key,
    chromosome,
    min_length_nt = 200
) {
  fasta <- safe_read_fasta(fasta_file)
  if (length(fasta) == 0) {
    warning("No sequences found/readable for accession: ", accession)
    return(tibble::tibble())
  }
  fasta <- fasta[vapply(fasta, is_valid_cds, logical(1), min_length_nt = min_length_nt)]
  if (length(fasta) == 0) {
    warning("All CDS filtered out for accession: ", accession)
    return(tibble::tibble())
  }
  rscu_tbl <- if (identical(rscu_method, "by_gene_mean")) {
    compute_rscu_by_gene_mean(fasta)
  } else {
    compute_rscu_by_counts(fasta)
  }
  rscu_tbl %>%
    dplyr::mutate(
      sample_id = sample_id,
      accession = accession,
      sample_type = sample_type,
      host_key = host_key,
      chromosome = chromosome
    ) %>%
    dplyr::select(
      sample_id,
      accession,
      sample_type,
      host_key,
      chromosome,
      codon,
      rscu
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

# ------------------------------------------
# Build organism table for RSCU calculation
# ------------------------------------------
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

# -----------------------------------------
# Fetch FASTAs and calculate per-CDS RSCU
# -----------------------------------------
msg("Fetching CDS FASTAs for ", nrow(organisms), " unique organism records...")

organisms_fetched <- organisms %>%
  dplyr::mutate(
    cds_fasta = purrr::map_chr(accession, fetch_cds_fasta)
  )

# phage RSCU
phage_rscu_raw <- organisms_fetched %>%
  dplyr::filter(sample_type == "phage") %>%
  dplyr::mutate(
    rscu_data = purrr::pmap(
      list(cds_fasta, sample_id, accession, sample_type, host_key, chromosome),
      compute_rscu_for_one_organism
    )
  ) %>%
  dplyr::select(rscu_data) %>%
  tidyr::unnest(rscu_data)

# host RSCU computed by host key - if two chromosomes exist, they are combined BEFORE RSCU calculation
compute_rscu_for_host_key <- function(
    host_key,
    host_records,
    min_length_nt = 200
) {
  fasta_files <- host_records$cds_fasta
  fasta_files <- fasta_files[!is.na(fasta_files)]
  if (length(fasta_files) == 0) {
    warning("No CDS FASTA files available for host_key: ", host_key)
    return(tibble::tibble())
  }
  fasta_list <- purrr::map(fasta_files, safe_read_fasta)
  combined_fasta <- unlist(fasta_list, recursive = FALSE)
  if (length(combined_fasta) == 0) {
    warning("No readable CDS sequences for host_key: ", host_key)
    return(tibble::tibble())
  }
  combined_fasta <- combined_fasta[
    vapply(combined_fasta, is_valid_cds, logical(1), min_length_nt = min_length_nt)
  ]
  if (length(combined_fasta) == 0) {
    warning("All CDS filtered out for host_key: ", host_key)
    return(tibble::tibble())
  }
  rscu_tbl <- if (identical(rscu_method, "by_gene_mean")) {
    compute_rscu_by_gene_mean(combined_fasta)
  } else {
    compute_rscu_by_counts(combined_fasta)
  }
  tibble::tibble(
    sample_id = host_key,
    accession = paste(host_records$accession, collapse = " + "),
    sample_type = "host",
    host_key = host_key,
    chromosome = paste(host_records$chromosome, collapse = " + ")
  ) %>%
    dplyr::bind_cols(rscu_tbl)
}

host_rscu_raw <- organisms_fetched %>%
  dplyr::filter(sample_type == "host") %>%
  dplyr::group_by(host_key) %>%
  dplyr::group_split() %>%
  purrr::map_dfr(function(df) {
    compute_rscu_for_host_key(
      host_key = df$host_key[[1]],
      host_records = df
    )
  })

# -----------------------------
# Combine host and phage RSCU
# -----------------------------
rscu_tidy <- dplyr::bind_rows(
  phage_rscu_raw,
  host_rscu_raw
) %>%
  dplyr::rename(RSCU = rscu)

readr::write_csv(
  rscu_tidy,
  file.path(out_dir, "rscu_tidy_per_organism.csv")
)

rscu_wide <- rscu_tidy %>%
  dplyr::select(sample_id, accession, sample_type, host_key, chromosome, codon, RSCU) %>%
  tidyr::pivot_wider(
    names_from = codon,
    values_from = RSCU,
    values_fill = list(RSCU = 0)
  ) %>%
  dplyr::arrange(sample_type, sample_id)

readr::write_csv(
  rscu_wide,
  file.path(out_dir, "rscu_wide_per_organism.csv")
)

rscu_filt <- rscu_tidy %>%
  dplyr::filter(codon %in% filtered_codons)

# --------------------------------------------------------
# Attach the RSCU data for phage and host to the metadata
# --------------------------------------------------------
phage_rscu <- rscu_tidy %>%
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

host_rscu <- rscu_tidy %>%
  dplyr::filter(sample_type == "host")

readr::write_csv(
  phage_rscu,
  file.path(out_dir, "phage_RSCU_with_metadata.csv")
)

readr::write_csv(
  host_rscu,
  file.path(out_dir, "host_RSCU.csv")
)

# ---------------------------------------------------------
# Calculate ∆RSCU (RSCU phage - RSCU host) for all phages
# ---------------------------------------------------------
msg("Computing matched-host ΔRSCU...")

host_rscu_by_key <- host_rscu %>%
  dplyr::select(
    host_key,
    codon,
    host_RSCU = RSCU
  )

phage_delta <- phage_rscu %>%
  dplyr::filter(!is.na(host_key)) %>%
  dplyr::left_join(
    host_rscu_by_key,
    by = c("host_key", "codon")
  ) %>%
  dplyr::mutate(
    Delta_RSCU = RSCU - host_RSCU
  )

readr::write_csv(
  phage_delta,
  file.path(out_dir, "delta_rscu_per_phage_per_codon.csv")
)

phage_delta_filt <- phage_delta %>%
  dplyr::filter(codon %in% filtered_codons)

phage_delta_filt <- phage_delta_filt %>%
  dplyr::mutate(
    tRNA_count_bin = dplyr::case_when(
      is.na(tRNAs) ~ NA_character_,
      tRNAs == 0 ~ "0",
      tRNAs > 0 & tRNAs <= 5 ~ "1-5",
      tRNAs > 5 & tRNAs <= 10 ~ "6-10",
      tRNAs > 10 ~ ">10"
    ),
    tRNA_count_bin = factor(
      tRNA_count_bin,
      levels = c("0", "1-5", "6-10", ">10")
    )
  )

# ========================================
# RSCU and ∆RSCU summaries and statistics
# ========================================

# ------------------------
# ∆RSCU summary per-phage
# ------------------------
phage_delta_scores <- phage_delta_filt %>%
  dplyr::group_by(
    sample_id,
    accession,
    host_key,
    host_genus,
    host_species,
    host_accession,
    host_accession2,
    lifestyle,
    tRNAs,
    has_tRNAs
  ) %>%
  dplyr::summarise(
    delta_MAE = mean(abs(Delta_RSCU), na.rm = TRUE),
    delta_EUC = sqrt(sum((Delta_RSCU)^2, na.rm = TRUE)),
    mean_delta_RSCU = mean(Delta_RSCU, na.rm = TRUE),
    median_delta_RSCU = median(Delta_RSCU, na.rm = TRUE),
    n_codons = sum(is.finite(Delta_RSCU)),
    .groups = "drop"
  ) %>%
  dplyr::rename(
    phage_id = sample_id,
    phage_accession = accession
  ) %>%
  dplyr::mutate(
    host_genus = factor(as.character(host_genus), levels = genus_order_final),
    lifestyle = factor(lifestyle, levels = lifestyle_order),
    has_tRNAs = factor(has_tRNAs, levels = trna_presence_order)
  )

readr::write_csv(
  phage_delta_scores,
  file.path(out_dir, "delta_rscu_scores_per_phage.csv")
)

# --------------------------------
# Genus-level per-codon summaries
# --------------------------------
delta_rscu_genus_summary <- phage_delta_scores %>%
  dplyr::filter(!is.na(delta_MAE), !is.na(host_genus)) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::summarise(
    n_phages = dplyr::n_distinct(phage_id),
    mean_delta_MAE = mean(delta_MAE, na.rm = TRUE),
    median_delta_MAE = median(delta_MAE, na.rm = TRUE),
    sd_delta_MAE = sd(delta_MAE, na.rm = TRUE),
    mean_delta_EUC = mean(delta_EUC, na.rm = TRUE),
    median_delta_EUC = median(delta_EUC, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  delta_rscu_genus_summary,
  file.path(summary_dir, "deltaRSCU_genus_summary.csv")
)

# -----------------------------
# Per-phage similarity to host
# -----------------------------
per_phage_rscu_similarity <- phage_delta_filt %>%
  dplyr::filter(
    is.finite(RSCU),
    is.finite(host_RSCU)
  ) %>%
  dplyr::group_by(
    sample_id,
    accession,
    host_key,
    host_genus,
    host_species,
    host_accession,
    host_accession2,
    lifestyle,
    tRNAs,
    has_tRNAs
  ) %>%
  dplyr::group_modify(~{
    df <- .x %>%
      dplyr::filter(
        is.finite(RSCU),
        is.finite(host_RSCU)
      )
    n_codons <- nrow(df)
    if (
      n_codons < 3 ||
      dplyr::n_distinct(df$RSCU) < 2 ||
      dplyr::n_distinct(df$host_RSCU) < 2
    ) {
      return(
        tibble::tibble(
          n_codons = n_codons,
          spearman_rho = NA_real_,
          p_value = NA_real_,
          delta_MAE = mean(abs(df$Delta_RSCU), na.rm = TRUE),
          delta_RMSE = sqrt(mean(df$Delta_RSCU^2, na.rm = TRUE)),
          delta_EUC = sqrt(sum(df$Delta_RSCU^2, na.rm = TRUE))
        )
      )
    }
    ct <- suppressWarnings(
      stats::cor.test(
        df$RSCU,
        df$host_RSCU,
        method = "spearman",
        exact = FALSE
      )
    )
    tibble::tibble(
      n_codons = n_codons,
      spearman_rho = unname(ct$estimate),
      p_value = ct$p.value,
      delta_MAE = mean(abs(df$Delta_RSCU), na.rm = TRUE),
      delta_RMSE = sqrt(mean(df$Delta_RSCU^2, na.rm = TRUE)),
      delta_EUC = sqrt(sum(df$Delta_RSCU^2, na.rm = TRUE))
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    p_adj_BH = stats::p.adjust(p_value, method = "BH")
  ) %>%
  dplyr::arrange(host_genus, sample_id) %>%
  dplyr::rename(
    phage_id = sample_id,
    phage_accession = accession
  )

readr::write_csv(
  per_phage_rscu_similarity,
  file.path(
    out_dir,
    "per_phage_RSCU_profile_similarity_to_host.csv"
  )
)

# ---------------------
# lifestyle comparison
# ---------------------
lifestyle <- phage_delta_scores %>%
  dplyr::filter(!is.na(delta_MAE), !is.na(lifestyle)) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::group_modify(~{
    df <- .x
    wt <- safe_wilcox_tbl(
      df$delta_MAE[df$lifestyle == "virulent"],
      df$delta_MAE[df$lifestyle == "temperate"]
    )
    tibble::tibble(
      test = "Wilcoxon",
      W = wt$W,
      p_value = wt$p_value,
      n_virulent = sum(df$lifestyle == "virulent", na.rm = TRUE),
      n_temperate = sum(df$lifestyle == "temperate", na.rm = TRUE),
      median_virulent = median(df$delta_MAE[df$lifestyle == "virulent"], na.rm = TRUE),
      median_temperate = median(df$delta_MAE[df$lifestyle == "temperate"], na.rm = TRUE),
      cliffs_delta = cliffs_delta(df$delta_MAE, df$lifestyle)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  lifestyle,
  file.path(summary_dir, "stats_lifestyle_deltaMAE.csv")
)

# -------------------------
# tRNA presence comparison
# -------------------------
trna_presence <- phage_delta_scores %>%
  dplyr::filter(!is.na(delta_MAE), !is.na(has_tRNAs)) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::group_modify(~{
    df <- .x
    wt <- safe_wilcox_tbl(
      df$delta_MAE[df$has_tRNAs == "Without tRNAs"],
      df$delta_MAE[df$has_tRNAs == "With tRNAs"]
    )
    tibble::tibble(
      test = "Wilcoxon",
      W = wt$W,
      p_value = wt$p_value,
      n_without_tRNAs = sum(df$has_tRNAs == "Without tRNAs", na.rm = TRUE),
      n_with_tRNAs = sum(df$has_tRNAs == "With tRNAs", na.rm = TRUE),
      median_without_tRNAs = median(df$delta_MAE[df$has_tRNAs == "Without tRNAs"], na.rm = TRUE),
      median_with_tRNAs = median(df$delta_MAE[df$has_tRNAs == "With tRNAs"], na.rm = TRUE),
      cliffs_delta = cliffs_delta(df$delta_MAE, df$has_tRNAs)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  trna_presence,
  file.path(summary_dir, "stats_tRNA_presence_deltaMAE.csv")
)

# ------------------------
# tRNA count correlations
# ------------------------
trna_numeric <- phage_delta_scores %>%
  dplyr::filter(!is.na(delta_MAE), !is.na(tRNAs)) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::group_modify(~{
    safe_spearman(.x$tRNAs, .x$delta_MAE)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  trna_numeric,
  file.path(summary_dir, "stats_numeric_tRNA_deltaMAE_spearman.csv")
)

# ======================
# Per-codon statistics
# ======================
safe_sign_p <- function(x, nresample = 10000) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  if (all(x == 0)) return(1)
  df <- data.frame(v = x)
  set.seed(123)
  p <- tryCatch(
    {
      wt <- coin::wilcoxsign_test(
        v ~ 1,
        data = df,
        distribution = coin::approximate(nresample = nresample)
      )
      as.numeric(coin::pvalue(wt))
    },
    error = function(e) {
      suppressWarnings(
        stats::wilcox.test(x, mu = 0, exact = FALSE, correct = TRUE)$p.value
      )
    }
  )
  p
}

percodon_raw <- phage_delta_filt %>%
  dplyr::select(
    sample_id,
    accession,
    host_key,
    host_genus,
    host_species,
    lifestyle,
    tRNAs,
    has_tRNAs,
    codon,
    RSCU,
    host_RSCU,
    Delta_RSCU
  ) %>%
  dplyr::arrange(host_genus, host_species, sample_id, codon)

readr::write_csv(
  percodon_raw,
  file.path(out_dir, "percodon_phage_vs_matched_host_RAW.csv")
)

# Per-codon phage vs host
percodon_phage_vs_host_within_genus <- phage_delta_filt %>%
  dplyr::group_by(host_genus, codon) %>%
  dplyr::summarise(
    n = sum(is.finite(Delta_RSCU)),
    median_delta = stats::median(Delta_RSCU, na.rm = TRUE),
    mean_delta = mean(Delta_RSCU, na.rm = TRUE),
    mad_delta = stats::mad(Delta_RSCU, na.rm = TRUE),
    p_value = safe_sign_p(Delta_RSCU, nresample = 10000),
    .groups = "drop"
  ) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(host_genus, codon)

readr::write_csv(
  percodon_phage_vs_host_within_genus,
  file.path(out_dir, "percodon_phage_vs_host_significance.csv")
)

# Per-codon lifestyle
percodon_lifestyle <- phage_delta_filt %>%
  dplyr::filter(!is.na(lifestyle)) %>%
  dplyr::group_by(host_genus, codon) %>%
  dplyr::group_modify(~{
    df <- .x
    wt <- safe_wilcox_tbl(
      df$Delta_RSCU[df$lifestyle == "virulent"],
      df$Delta_RSCU[df$lifestyle == "temperate"]
    )
    tibble::tibble(
      n_virulent = sum(df$lifestyle == "virulent" & is.finite(df$Delta_RSCU), na.rm = TRUE),
      n_temperate = sum(df$lifestyle == "temperate" & is.finite(df$Delta_RSCU), na.rm = TRUE),
      med_virulent = median(df$Delta_RSCU[df$lifestyle == "virulent"], na.rm = TRUE),
      med_temperate = median(df$Delta_RSCU[df$lifestyle == "temperate"], na.rm = TRUE),
      cliffs_delta = cliffs_delta(df$Delta_RSCU, df$lifestyle),
      W = wt$W,
      p_value = wt$p_value
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(host_genus) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(host_genus, codon)

readr::write_csv(
  percodon_lifestyle,
  file.path(out_dir, "percodon_lifestyle_wilcox.csv")
)

# Per-codon tRNA presence
percodon_trna_presence <- phage_delta_filt %>%
  dplyr::filter(!is.na(has_tRNAs)) %>%
  dplyr::group_by(host_genus, codon) %>%
  dplyr::group_modify(~{
    df <- .x
    wt <- safe_wilcox_tbl(
      df$Delta_RSCU[df$has_tRNAs == "Without tRNAs"],
      df$Delta_RSCU[df$has_tRNAs == "With tRNAs"]
    )
    tibble::tibble(
      n_without_tRNAs = sum(df$has_tRNAs == "Without tRNAs" & is.finite(df$Delta_RSCU), na.rm = TRUE),
      n_with_tRNAs = sum(df$has_tRNAs == "With tRNAs" & is.finite(df$Delta_RSCU), na.rm = TRUE),
      med_without_tRNAs = median(df$Delta_RSCU[df$has_tRNAs == "Without tRNAs"], na.rm = TRUE),
      med_with_tRNAs = median(df$Delta_RSCU[df$has_tRNAs == "With tRNAs"], na.rm = TRUE),
      cliffs_delta = cliffs_delta(df$Delta_RSCU, df$has_tRNAs),
      W = wt$W,
      p_value = wt$p_value
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(host_genus) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(host_genus, codon)

readr::write_csv(
  percodon_trna_presence,
  file.path(out_dir, "percodon_tRNA_presence_wilcox.csv")
)

# Per-codon numeric tRNA count
percodon_trna_numeric <- phage_delta_filt %>%
  dplyr::filter(!is.na(tRNAs)) %>%
  dplyr::group_by(host_genus, codon) %>%
  dplyr::group_modify(~{
    safe_spearman(.x$tRNAs, .x$Delta_RSCU)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(host_genus) %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(host_genus, codon)

readr::write_csv(
  percodon_trna_numeric,
  file.path(out_dir, "percodon_numeric_tRNA_spearman.csv")
)

# Per-codon numeric tRNA count vs ABSOLUTE ΔRSCU
percodon_trna_numeric_abs <- phage_delta_filt %>%
  dplyr::filter(
    !is.na(tRNAs),
    is.finite(Delta_RSCU)
  ) %>%
  dplyr::mutate(
    abs_Delta_RSCU = abs(Delta_RSCU)
  ) %>%
  dplyr::group_by(
    host_genus,
    codon
  ) %>%
  dplyr::group_modify(~{
    safe_spearman(
      .x$tRNAs,
      .x$abs_Delta_RSCU
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(host_genus) %>%
  dplyr::mutate(
    p_adj_BH = stats::p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    host_genus,
    codon
  )

readr::write_csv(
  percodon_trna_numeric_abs,
  file.path(
    out_dir,
    "percodon_numeric_tRNA_ABS_deltaRSCU_WITHIN_GENUS_spearman.csv"
  )
)

# ----------------------
# ∆ RSCU Summary tables
# ----------------------
codon_levels <- filtered_codons

# lifestyle
delta_life_long <- phage_delta_filt %>%
  dplyr::filter(!is.na(lifestyle)) %>%
  dplyr::group_by(host_genus, codon, lifestyle) %>%
  dplyr::summarise(
    n_phages = sum(is.finite(Delta_RSCU)),
    mean_delta = mean(Delta_RSCU, na.rm = TRUE),
    median_delta = stats::median(Delta_RSCU, na.rm = TRUE),
    mad_delta = stats::mad(Delta_RSCU, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(codon = factor(codon, levels = codon_levels)) %>%
  dplyr::arrange(host_genus, codon, lifestyle)

readr::write_csv(
  delta_life_long,
  file.path(out_dir, "deltaRSCU_lifestyle_LONG.csv")
)

delta_life_mean_wide <- delta_life_long %>%
  dplyr::select(host_genus, codon, lifestyle, mean_delta) %>%
  tidyr::pivot_wider(names_from = lifestyle, values_from = mean_delta)

readr::write_csv(
  delta_life_mean_wide,
  file.path(out_dir, "deltaRSCU_lifestyle_MEAN_WIDE.csv")
)

# tRNA presence
delta_trna_presence_long <- phage_delta_filt %>%
  dplyr::filter(!is.na(has_tRNAs)) %>%
  dplyr::group_by(host_genus, codon, has_tRNAs) %>%
  dplyr::summarise(
    n_phages = sum(is.finite(Delta_RSCU)),
    mean_delta = mean(Delta_RSCU, na.rm = TRUE),
    median_delta = stats::median(Delta_RSCU, na.rm = TRUE),
    mad_delta = stats::mad(Delta_RSCU, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(codon = factor(codon, levels = codon_levels)) %>%
  dplyr::arrange(host_genus, codon, has_tRNAs)

readr::write_csv(
  delta_trna_presence_long,
  file.path(out_dir, "deltaRSCU_tRNA_presence_LONG.csv")
)

delta_trna_presence_mean_wide <- delta_trna_presence_long %>%
  dplyr::select(host_genus, codon, has_tRNAs, mean_delta) %>%
  tidyr::pivot_wider(names_from = has_tRNAs, values_from = mean_delta)

readr::write_csv(
  delta_trna_presence_mean_wide,
  file.path(out_dir, "deltaRSCU_tRNA_presence_MEAN_WIDE.csv")
)

# ============================
# PERMANOVA on ∆RSCU vectors
# ============================
phage_delta_wide <- phage_delta_filt %>%
  dplyr::select(
    sample_id,
    codon,
    Delta_RSCU,
    host_genus,
    host_species,
    lifestyle,
    tRNAs,
    has_tRNAs
  ) %>%
  tidyr::pivot_wider(
    names_from = codon,
    values_from = Delta_RSCU,
    values_fill = 0
  )

phage_delta_mat <- phage_delta_wide %>%
  dplyr::select(
    -sample_id,
    -host_genus,
    -host_species,
    -lifestyle,
    -tRNAs,
    -has_tRNAs
  ) %>%
  as.matrix()

rownames(phage_delta_mat) <- phage_delta_wide$sample_id

permanova_df <- phage_delta_wide %>%
  dplyr::select(
    sample_id,
    host_genus,
    host_species,
    lifestyle,
    tRNAs,
    has_tRNAs
  ) %>%
  dplyr::mutate(
    lifestyle = factor(lifestyle, levels = lifestyle_order),
    host_genus = factor(host_genus, levels = genus_order_final),
    tRNAs = as.numeric(tRNAs)
  )

dist_delta <- stats::dist(
  scale(phage_delta_mat, center = TRUE, scale = TRUE),
  method = "euclidean"
)

set.seed(123)

perm_across <- vegan::adonis2(
  dist_delta ~ lifestyle + tRNAs + host_genus,
  data = permanova_df,
  permutations = 4999,
  by = "margin",
  strata = permanova_df$host_genus
)

sink(file.path(out_dir, "permanova_across_marginal_numeric_tRNAs.txt"))
print(perm_across)
sink()

# ------------------
# Dispersion checks
# ------------------
bd_life <- vegan::betadisper(dist_delta, group = permanova_df$lifestyle)
bd_host <- vegan::betadisper(dist_delta, group = permanova_df$host_genus)

permanova_df <- permanova_df %>%
  dplyr::mutate(
    tRNA_count_bin = dplyr::case_when(
      is.na(tRNAs) ~ NA_character_,
      tRNAs == 0 ~ "0",
      tRNAs > 0 & tRNAs <= 5 ~ "1-5",
      tRNAs > 5 & tRNAs <= 10 ~ "6-10",
      tRNAs > 10 ~ ">10"
    ),
    tRNA_count_bin = factor(tRNA_count_bin, levels = c("0", "1-5", "6-10", ">10"))
  )

bd_trna_bin <- vegan::betadisper(
  dist_delta,
  group = permanova_df$tRNA_count_bin
)

sink(file.path(out_dir, "permanova_dispersion_numeric_tRNAs.txt"))

cat("\n-- Dispersion lifestyle --\n")
print(anova(bd_life))
print(vegan::permutest(bd_life, permutations = 4999))

cat("\n-- Dispersion host genus --\n")
print(anova(bd_host))
print(vegan::permutest(bd_host, permutations = 4999))

cat("\n-- Dispersion tRNA count bins, diagnostic only --\n")
print(anova(bd_trna_bin))
print(vegan::permutest(bd_trna_bin, permutations = 4999))

sink()

# ======
# PLOTS
# ======

# -----------------------------
# ∆RSCU MAE by lifestyle
# -----------------------------
p_lifestyle <- phage_delta_scores %>%
  dplyr::filter(!is.na(lifestyle)) %>%
  ggplot(aes(x = lifestyle, y = delta_MAE, fill = lifestyle)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~ host_genus, scales = "free_y", drop = TRUE) +
  scale_fill_manual(
    values = palette_lifestyle,
    labels = c(virulent = "Virulent", temperate = "Temperate"),
    name = "Lifestyle"
  ) +
  scale_x_discrete(labels = c(virulent = "Virulent", temperate = "Temperate")) +
  labs(
    title = "ΔRSCU MAE by Lifestyle within Host Genus",
    x = NULL,
    y = "Mean |ΔRSCU| across codons"
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank()
  )

ggsave(
  file.path(plot_dir, "plot_deltaMAE_by_lifestyle.png"),
  p_lifestyle,
  width = 12,
  height = 8,
  dpi = 300
)

# -----------------------------
# ∆RSCU MAE by tRNA presence
# -----------------------------
p_trna_presence <- phage_delta_scores %>%
  dplyr::filter(!is.na(has_tRNAs)) %>%
  ggplot(aes(x = has_tRNAs, y = delta_MAE, fill = has_tRNAs)) +
  geom_violin(trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~ host_genus, scales = "free_y", drop = TRUE) +
  scale_fill_manual(
    values = palette_trna_presence,
    name = "tRNA presence",
    drop = FALSE
  ) +
  labs(
    title = "ΔRSCU MAE by tRNA Presence within Host Genus",
    x = NULL,
    y = "Mean |ΔRSCU| across codons"
  ) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )


ggsave(
  file.path(plot_dir, "plot_deltaMAE_by_tRNA_presence.png"),
  p_trna_presence,
  width = 12,
  height = 8,
  dpi = 300
)

# -----------------------------
# ∆RSCU MAE by numeric tRNA count
# -----------------------------
p_trna_numeric <- phage_delta_scores %>%
  dplyr::filter(!is.na(tRNAs)) %>%
  ggplot(aes(x = tRNAs, y = delta_MAE)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  facet_wrap(~ host_genus, scales = "free", drop = TRUE) +
  labs(
    title = "ΔRSCU MAE by Numeric tRNA Count",
    x = "Number of phage-encoded tRNAs",
    y = "Mean |ΔRSCU| across codons"
  ) +
  theme(
    panel.grid = element_blank()
  )

ggsave(
  file.path(plot_dir, "plot_deltaMAE_by_numeric_tRNA_count.png"),
  p_trna_numeric,
  width = 12,
  height = 8,
  dpi = 300
)

# ====================================
# Per-genus raw RSCU heatmaps:
# matched host genome(s) + all phages
# ====================================
# Host rows
raw_host_heatmap <- phage_delta_filt %>%
  dplyr::select(
    host_genus,
    host_species,
    host_key,
    codon,
    host_RSCU
  ) %>%
  dplyr::distinct() %>%
  dplyr::transmute(
    host_genus = as.character(host_genus),
    sample_id = paste0("HOST_", host_key),
    sample_label = paste0(
      "HOST | ",
      dplyr::coalesce(host_species, host_key),
      " | ",
      host_key
    ),
    sample_type = "Host",
    codon = codon,
    RSCU = host_RSCU
  )

# Phage rows
raw_phage_heatmap <- phage_delta_filt %>%
  dplyr::select(
    host_genus,
    sample_id,
    codon,
    RSCU
  ) %>%
  dplyr::distinct() %>%
  dplyr::transmute(
    host_genus = as.character(host_genus),
    sample_id = paste0("PHAGE_", sample_id),
    sample_label = paste0(
      "PHAGE | ",
      stringr::str_remove(sample_id, "^PHAGE_")
    ),
    sample_type = "Phage",
    codon = codon,
    RSCU = RSCU
  )

# --------------------------------
# Combine host and phage profiles
# --------------------------------
raw_rscu_heatmap_long <- dplyr::bind_rows(
  raw_host_heatmap,
  raw_phage_heatmap
) %>%
  dplyr::filter(
    !is.na(host_genus),
    codon %in% filtered_codons,
    is.finite(RSCU)
  ) %>%
  dplyr::mutate(
    codon = factor(codon, levels = filtered_codons),
    sample_type = factor(
      sample_type,
      levels = c("Host", "Phage")
    )
  )

readr::write_csv(
  raw_rscu_heatmap_long,
  file.path(
    out_dir,
    "raw_RSCU_host_and_phages_heatmap_LONG.csv"
  )
)

# --------------------------------
# Plot one heatmap per host genus
# --------------------------------
for (gn in unique(na.omit(raw_rscu_heatmap_long$host_genus))) {
  sub <- raw_rscu_heatmap_long %>%
    dplyr::filter(host_genus == gn)
  row_info <- sub %>%
    dplyr::distinct(
      sample_id,
      sample_label,
      sample_type
    ) %>%
    dplyr::arrange(
      sample_type,
      sample_label
    )
  sub_wide <- sub %>%
    dplyr::select(
      sample_id,
      sample_label,
      sample_type,
      codon,
      RSCU
    ) %>%
    tidyr::pivot_wider(
      names_from = codon,
      values_from = RSCU,
      values_fill = list(RSCU = 0)
    ) %>%
    dplyr::left_join(
      row_info,
      by = c(
        "sample_id",
        "sample_label",
        "sample_type"
      )
    ) %>%
    dplyr::arrange(
      sample_type,
      sample_label
    )
  codons_present <- filtered_codons[
    filtered_codons %in% colnames(sub_wide)
  ]
  mat <- sub_wide %>%
    dplyr::select(dplyr::all_of(codons_present)) %>%
    as.matrix()
  rownames(mat) <- sub_wide$sample_label
  row_annotation <- data.frame(
    Genome_Type = as.character(sub_wide$sample_type)
  )
  rownames(row_annotation) <- sub_wide$sample_label
  safe_genus <- clean_name_for_file(gn)
  heatmap_height <- max(
    6,
    2.5 + 0.22 * nrow(mat)
  )
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    main = paste0(
      "Raw RSCU Profiles — ",
      gn,
      "\nMatched Host Genome(s) and Phages"
    ),
    border_color = NA,
    fontsize = 9,
    fontsize_row = 7,
    fontsize_col = 7,
    angle_col = 90,
    filename = file.path(
      plot_dir,
      paste0(
        "heatmap_raw_RSCU_host_and_phages_",
        safe_genus,
        ".png"
      )
    ),
    width = 14,
    height = heatmap_height
  )
}

# ---------------
# Other Heatmaps
# ---------------
# ∆RSCU by lifestyle
heat_df_life <- phage_delta_filt %>%
  dplyr::filter(!is.na(lifestyle)) %>%
  dplyr::group_by(host_genus, lifestyle, codon) %>%
  dplyr::summarise(
    mean_delta = mean(Delta_RSCU, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(codon = factor(codon, levels = filtered_codons))

for (gn in unique(na.omit(heat_df_life$host_genus))) {
  sub <- heat_df_life %>%
    dplyr::filter(host_genus == gn) %>%
    tidyr::pivot_wider(
      names_from = codon,
      values_from = mean_delta,
      values_fill = 0
    ) %>%
    dplyr::arrange(lifestyle)
  mat <- as.matrix(sub[, setdiff(colnames(sub), c("host_genus", "lifestyle"))])
  rownames(mat) <- tools::toTitleCase(as.character(sub$lifestyle))
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste0("Mean ΔRSCU by Codon — ", gn, " (Lifestyle)"),
    filename = file.path(plot_dir, paste0("heatmap_mean_delta_by_lifestyle_", gn, ".png")),
    width = 14,
    height = 4
  )
}

# ∆RSCU by tRNA presence
heat_df_trna <- phage_delta_filt %>%
  dplyr::filter(!is.na(has_tRNAs)) %>%
  dplyr::group_by(host_genus, has_tRNAs, codon) %>%
  dplyr::summarise(
    mean_delta = mean(Delta_RSCU, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(codon = factor(codon, levels = filtered_codons))

for (gn in unique(na.omit(heat_df_trna$host_genus))) {
  sub <- heat_df_trna %>%
    dplyr::filter(host_genus == gn) %>%
    tidyr::pivot_wider(
      names_from = codon,
      values_from = mean_delta,
      values_fill = 0
    ) %>%
    dplyr::arrange(has_tRNAs)
  mat <- as.matrix(sub[, setdiff(colnames(sub), c("host_genus", "has_tRNAs"))])
  rownames(mat) <- as.character(sub$has_tRNAs)
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste0("Mean ΔRSCU by Codon — ", gn, " (tRNA Presence)"),
    filename = file.path(plot_dir, paste0("heatmap_mean_delta_by_tRNA_presence_", gn, ".png")),
    width = 14,
    height = 4
  )
}

# ∆RSCU by tRNA number, binned
heat_df_trna_bins <- phage_delta_filt %>%
  dplyr::filter(!is.na(tRNA_count_bin), is.finite(Delta_RSCU)) %>%
  dplyr::group_by(host_genus, tRNA_count_bin, codon) %>%
  dplyr::summarise(
    mean_delta = mean(Delta_RSCU, na.rm = TRUE),
    median_delta = median(Delta_RSCU, na.rm = TRUE),
    n_phages = dplyr::n_distinct(sample_id),
    .groups = "drop"
  ) %>%
  dplyr::mutate(codon = factor(codon, levels = filtered_codons))

readr::write_csv(
  heat_df_trna_bins,
  file.path(out_dir, "deltaRSCU_numeric_tRNA_bins_LONG.csv")
)

for (gn in unique(na.omit(heat_df_trna_bins$host_genus))) {
  sub <- heat_df_trna_bins %>%
    dplyr::filter(host_genus == gn) %>%
    tidyr::pivot_wider(
      names_from = codon,
      values_from = mean_delta,
      values_fill = 0
    ) %>%
    dplyr::arrange(tRNA_count_bin)
  mat <- as.matrix(sub[, setdiff(colnames(sub), c("host_genus", "tRNA_count_bin"))])
  rownames(mat) <- as.character(sub$tRNA_count_bin)
  pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    main = paste0("Mean ΔRSCU by Codon — ", gn, " (tRNA Count Bin)"),
    filename = file.path(plot_dir, paste0("heatmap_mean_delta_by_numeric_tRNA_bins_", gn, ".png")),
    width = 14,
    height = 4
  )
}

# ======
# DONE!
# ======
msg("All done! :) Outputs saved in: ", out_dir)
