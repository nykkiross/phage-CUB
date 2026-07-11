# ================================
# EnC, ∆EnC, and EnC/GC3 ANALYSIS
# ================================

# ----------------
# Settings/inputs
# ----------------
metadata_xlsx <- "/your/file/path/sample_phages.xlsx"

out_dir <- "/your/directory/path/Results/EnC_deltaEnC"

# Folder where downloaded CDS FASTA files will be cached
cds_cache_dir <- file.path(out_dir, "downloaded_CDS_fastas")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cds_cache_dir, recursive = TRUE, showWarnings = FALSE)

out_sum <- file.path(out_dir, "summaries")
dir.create(out_sum, recursive = TRUE, showWarnings = FALSE)

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
  "broom",
  "openxlsx"
)

not_installed <- req_pkgs[!req_pkgs %in% installed.packages()[, "Package"]]

if (length(not_installed) > 0) {
  install.packages(not_installed, dependencies = TRUE)
}

invisible(lapply(req_pkgs, library, character.only = TRUE))

# ======================================================
# Helper functions for accessions and statistical tests
# ======================================================
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
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

safe_wilcox <- function(x, y) {
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

# ==========================
# EnC calculation functions
# ==========================
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
    seqinr::read.fasta(path, seqtype = "DNA"),
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

compute_sample_Ks <- function(genome, pseudocount = 1.0) {
  first_counts <- seqinr::uco(genome[[1]], index = "eff", frame = 0)
  total_counts <- rep(0, length(first_counts))
  names(total_counts) <- toupper(names(first_counts))
  
  for (seq in genome) {
    cc <- seqinr::uco(seq, index = "eff", frame = 0)
    cc <- cc[match(names(total_counts), toupper(names(cc)))]
    cc[is.na(cc)] <- 0
    total_counts <- total_counts + cc
  }
  codon_table <- list(
    K = c("AAA", "AAG"),
    N = c("AAC", "AAT"),
    Q = c("CAA", "CAG"),
    H = c("CAC", "CAT"),
    D = c("GAC", "GAT"),
    E = c("GAA", "GAG"),
    Y = c("TAC", "TAT"),
    C = c("TGC", "TGT"),
    F = c("TTC", "TTT"),
    L2 = c("TTA", "TTG"),
    
    I = c("ATA", "ATC", "ATT"),
    
    T = c("ACA", "ACC", "ACG", "ACT"),
    P = c("CCA", "CCC", "CCG", "CCT"),
    A = c("GCA", "GCC", "GCG", "GCT"),
    G = c("GGA", "GGC", "GGG", "GGT"),
    V = c("GTA", "GTC", "GTG", "GTT"),
    R = c("CGA", "CGC", "CGG", "CGT"),
    L4 = c("CTA", "CTC", "CTG", "CTT"),
    S4 = c("AGC", "AGT", "TCA", "TCC", "TCG", "TCT")
  )
  eps <- 1e-8
  
  F_values <- sapply(codon_table, function(codons) {
    counts <- total_counts[codons]
    counts[is.na(counts)] <- 0
    counts <- counts + pseudocount
    n <- sum(counts)
    
    if (n <= 1) {
      return(NA_real_)
    }
    p <- counts / n
    sum_f2 <- sum(p * p)
    Fcf <- (n * sum_f2 - 1) / (n - 1)
    Fcf <- max(min(Fcf, 1 - eps), eps)
    Fcf
  })
  list(
    K2 = mean(F_values[c("K", "N", "Q", "H", "D", "E", "Y", "C", "F", "L2")], na.rm = TRUE),
    K3 = mean(F_values[c("I")], na.rm = TRUE),
    K4 = mean(F_values[c("T", "P", "A", "G", "V", "R", "L4", "S4")], na.rm = TRUE)
  )
}

calculate_enc_with_fallback <- function(seq, sample_Ks, pseudocount = 1.0) {
  codon_counts <- seqinr::uco(seq, index = "eff", frame = 0)
  names(codon_counts) <- toupper(names(codon_counts))
  
  codon_table <- list(
    K = c("AAA", "AAG"),
    N = c("AAC", "AAT"),
    Q = c("CAA", "CAG"),
    H = c("CAC", "CAT"),
    D = c("GAC", "GAT"),
    E = c("GAA", "GAG"),
    Y = c("TAC", "TAT"),
    C = c("TGC", "TGT"),
    F = c("TTC", "TTT"),
    L2 = c("TTA", "TTG"),
    
    I = c("ATA", "ATC", "ATT"),
    
    T = c("ACA", "ACC", "ACG", "ACT"),
    P = c("CCA", "CCC", "CCG", "CCT"),
    A = c("GCA", "GCC", "GCG", "GCT"),
    G = c("GGA", "GGC", "GGG", "GGT"),
    V = c("GTA", "GTC", "GTG", "GTT"),
    R = c("CGA", "CGC", "CGG", "CGT"),
    L4 = c("CTA", "CTC", "CTG", "CTT"),
    S4 = c("AGC", "AGT", "TCA", "TCC", "TCG", "TCT")
  )
  eps <- 1e-8
  F_values <- sapply(codon_table, function(codons) {
    raw <- codon_counts[codons]
    if (all(is.na(raw)) || sum(raw, na.rm = TRUE) == 0) {
      return(NA_real_)
    }
    counts <- raw
    counts[is.na(counts)] <- 0
    counts <- counts + pseudocount
    n <- sum(counts)
    if (n <= 1) {
      return(NA_real_)
    }
    p <- counts / n
    sum_f2 <- sum(p * p)
    Fcf <- (n * sum_f2 - 1) / (n - 1)
    max(min(Fcf, 1 - eps), eps)
  })
  
  K2 <- mean(F_values[c("K", "N", "Q", "H", "D", "E", "Y", "C", "F", "L2")], na.rm = TRUE)
  K3 <- mean(F_values[c("I")], na.rm = TRUE)
  K4 <- mean(F_values[c("T", "P", "A", "G", "V", "R", "L4", "S4")], na.rm = TRUE)
  
  if (is.na(K2)) K2 <- sample_Ks$K2
  if (is.na(K3)) K3 <- sample_Ks$K3
  if (is.na(K4)) K4 <- sample_Ks$K4
  
  K2 <- ifelse(is.na(K2), NA_real_, pmax(K2, eps))
  K3 <- ifelse(is.na(K3), NA_real_, pmax(K3, eps))
  K4 <- ifelse(is.na(K4), NA_real_, pmax(K4, eps))
  
  if (any(is.na(c(K2, K3, K4)))) {
    return(NA_real_)
  }
  enc <- 2 + 9 / K2 + 1 / K3 + 5 / K4
  enc <- min(max(enc, 20), 61)
  enc
}

compute_enc_gc3_for_one_organism <- function(
    fasta_file,
    sample_id,
    accession,
    sample_type,
    host_key,
    chromosome,
    min_length_nt = 200,
    pseudocount = 0.5
) {
  genome <- safe_read_fasta(fasta_file)
  if (length(genome) == 0) {
    warning("No sequences found/readable for accession: ", accession)
    return(tibble::tibble())
  }
  genome <- genome[vapply(genome, is_valid_cds, logical(1), min_length_nt = min_length_nt)]
  if (length(genome) == 0) {
    warning("All CDS filtered out for accession: ", accession)
    return(tibble::tibble())
  }
  sample_Ks <- compute_sample_Ks(genome, pseudocount = pseudocount)
  enc_gc3_list <- lapply(genome, function(seq) {
    enc <- calculate_enc_with_fallback(seq, sample_Ks, pseudocount = pseudocount)
    gc3 <- seqinr::GC3(seq)
    c(EnC = enc, GC3 = gc3)
  })
  df <- as.data.frame(do.call(rbind, enc_gc3_list)) %>%
    tibble::rownames_to_column("gene_id") %>%
    mutate(
      gene_id = names(genome),
      sample_id = sample_id,
      accession = accession,
      sample_type = sample_type,
      host_key = host_key,
      chromosome = chromosome
    ) %>%
    select(
      gene_id,
      sample_id,
      accession,
      sample_type,
      host_key,
      chromosome,
      EnC,
      GC3
    )
  df
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

# -----------------------------------------
# Build organism table for EnC calculation
# -----------------------------------------
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

# ---------------------------------------------------
# Fetch FASTAs and calculate per-CDS EnC and EnC/GC3
# ---------------------------------------------------
all_gene <- organisms %>%
  mutate(
    cds_fasta = purrr::map_chr(accession, fetch_cds_fasta)
  ) %>%
  mutate(
    enc_data = purrr::pmap(
      list(cds_fasta, sample_id, accession, sample_type, host_key, chromosome),
      compute_enc_gc3_for_one_organism
    )
  ) %>%
  select(enc_data) %>%
  tidyr::unnest(enc_data)

readr::write_csv(
  all_gene,
  file.path(out_dir, "gene_level_EnC_GC3_all_organisms.csv")
)

all_gene %>%
  split(.$sample_type) %>%
  purrr::iwalk(function(df, st) {
    readr::write_csv(
      df,
      file.path(out_dir, paste0("gene_level_EnC_GC3_", st, ".csv"))
    )
  })

# -------------------------------------------------------------------
# Attach the EnC and EnC/GC3 data for phage and host to the metadata
# -------------------------------------------------------------------
phage_genes <- all_gene %>%
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

host_genes <- all_gene %>%
  filter(sample_type == "host")

readr::write_csv(
  phage_genes,
  file.path(out_dir, "phage_gene_level_EnC_GC3_with_metadata.csv")
)

readr::write_csv(
  host_genes,
  file.path(out_dir, "host_gene_level_EnC_GC3.csv")
)

# -------------------
# Sample-level means
# -------------------
host_sample_means <- host_genes %>%
  group_by(sample_id, accession, sample_type, host_key) %>%
  summarise(
    mean_EnC = mean(EnC, na.rm = TRUE),
    median_EnC = median(EnC, na.rm = TRUE),
    mean_GC3 = mean(GC3, na.rm = TRUE),
    median_GC3 = median(GC3, na.rm = TRUE),
    n_genes = sum(is.finite(EnC)),
    .groups = "drop"
  )

phage_sample_means <- phage_genes %>%
  group_by(
    sample_id,
    accession,
    sample_type,
    host_key,
    host_genus,
    host_species,
    host_accession,
    host_accession2,
    lifestyle,
    tRNAs,
    has_tRNAs
  ) %>%
  summarise(
    mean_EnC = mean(EnC, na.rm = TRUE),
    median_EnC = median(EnC, na.rm = TRUE),
    mean_GC3 = mean(GC3, na.rm = TRUE),
    median_GC3 = median(GC3, na.rm = TRUE),
    n_genes = sum(is.finite(EnC)),
    .groups = "drop"
  )

sample_means <- bind_rows(
  host_sample_means,
  phage_sample_means
)

readr::write_csv(
  sample_means,
  file.path(out_dir, "sample_means_EnC_GC3.csv")
)

# -----------------------------------------------------
# Calculate ∆EnC (EnC phage - EnC host) for all phages
# -----------------------------------------------------
per_phage_results <- phage_sample_means %>%
  filter(!is.na(host_key)) %>%
  group_by(
    sample_id,
    accession,
    host_key,
    host_accession,
    host_accession2,
    host_genus,
    host_species,
    lifestyle,
    tRNAs,
    has_tRNAs
  ) %>%
  group_modify(~{
    phage_row <- .x
    keys <- .y
    phage_id <- keys$sample_id[[1]]
    hk <- keys$host_key[[1]]
    ph_g <- phage_genes %>%
      filter(sample_id == phage_id)
    host_g <- host_genes %>%
      filter(host_key == hk)
    if (nrow(ph_g) < 2 || nrow(host_g) < 2) {
      return(tibble::tibble(
        mean_EnC_phage = phage_row$mean_EnC[[1]],
        median_EnC_phage = phage_row$median_EnC[[1]],
        mean_GC3_phage = phage_row$mean_GC3[[1]],
        mean_EnC_host = NA_real_,
        median_EnC_host = NA_real_,
        mean_GC3_host = NA_real_,
        delta_EnC = NA_real_,
        abs_delta_EnC = NA_real_,
        delta_GC3 = NA_real_,
        W = NA_real_,
        p_value = NA_real_,
        n_phage_genes = nrow(ph_g),
        n_host_genes = nrow(host_g)
      ))
    }
    host_mean_enc <- mean(host_g$EnC, na.rm = TRUE)
    host_median_enc <- median(host_g$EnC, na.rm = TRUE)
    host_mean_gc3 <- mean(host_g$GC3, na.rm = TRUE)
    phage_mean_enc <- phage_row$mean_EnC[[1]]
    phage_median_enc <- phage_row$median_EnC[[1]]
    phage_mean_gc3 <- phage_row$mean_GC3[[1]]
    
    wt <- safe_wilcox(ph_g$EnC, host_g$EnC)
    tibble::tibble(
      mean_EnC_phage = phage_mean_enc,
      median_EnC_phage = phage_median_enc,
      mean_GC3_phage = phage_mean_gc3,
      
      mean_EnC_host = host_mean_enc,
      median_EnC_host = host_median_enc,
      mean_GC3_host = host_mean_gc3,
      
      delta_EnC = phage_mean_enc - host_mean_enc,
      abs_delta_EnC = abs(phage_mean_enc - host_mean_enc),
      delta_GC3 = phage_mean_gc3 - host_mean_gc3,
      
      W = wt$W,
      p_value = wt$p_value,
      n_phage_genes = nrow(ph_g),
      n_host_genes = nrow(host_g)
    )
  }) %>%
  ungroup() %>%
  rename(
    phage_id = sample_id,
    phage_accession = accession
  ) %>%
  mutate(
    host_genus = factor(as.character(host_genus), levels = genus_order_final),
    lifestyle = factor(lifestyle, levels = c("virulent", "temperate")),
    has_tRNAs = factor(has_tRNAs, levels = c("Without tRNAs", "With tRNAs")),
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(host_genus, host_species, phage_id)

readr::write_csv(
  per_phage_results,
  file.path(out_sum, "per_phage_vs_matched_host_deltaEnC.csv")
)

# ================================
# Statistical summaries and tests
# ================================

# ----------
# Summaries
# ----------
enc_genus_summary <- per_phage_results %>%
  filter(!is.na(delta_EnC), !is.na(host_genus)) %>%
  group_by(host_genus) %>%
  summarise(
    n_phages = n_distinct(phage_id),
    mean_phage_EnC = mean(mean_EnC_phage, na.rm = TRUE),
    median_phage_EnC = median(mean_EnC_phage, na.rm = TRUE),
    mean_host_EnC = mean(mean_EnC_host, na.rm = TRUE),
    median_host_EnC = median(mean_EnC_host, na.rm = TRUE),
    mean_delta_EnC = mean(delta_EnC, na.rm = TRUE),
    median_delta_EnC = median(delta_EnC, na.rm = TRUE),
    sd_delta_EnC = sd(delta_EnC, na.rm = TRUE),
    mean_abs_delta_EnC = mean(abs_delta_EnC, na.rm = TRUE),
    median_abs_delta_EnC = median(abs_delta_EnC, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  enc_genus_summary,
  file.path(out_sum, "EnC_deltaEnC_genus_summary.csv")
)

# lifestyle summary
lifestyle_tests <- per_phage_results %>%
  filter(!is.na(abs_delta_EnC), !is.na(lifestyle)) %>%
  group_by(host_genus) %>%
  group_modify(~{
    df <- .x
    wt <- safe_wilcox(
      df$abs_delta_EnC[df$lifestyle == "virulent"],
      df$abs_delta_EnC[df$lifestyle == "temperate"]
    )
    tibble::tibble(
      test = "Wilcoxon",
      W = wt$W,
      p_value = wt$p_value,
      n_virulent = sum(df$lifestyle == "virulent", na.rm = TRUE),
      n_temperate = sum(df$lifestyle == "temperate", na.rm = TRUE),
      median_virulent = median(df$abs_delta_EnC[df$lifestyle == "virulent"], na.rm = TRUE),
      median_temperate = median(df$abs_delta_EnC[df$lifestyle == "temperate"], na.rm = TRUE),
      cliffs_delta = cliffs_delta(df$abs_delta_EnC, df$lifestyle)
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  lifestyle_tests,
  file.path(out_sum, "lifestyle_Wilcoxon_absDeltaEnC.csv")
)

# tRNA presence/absence summary
trna_presence_tests <- per_phage_results %>%
  filter(!is.na(abs_delta_EnC), !is.na(has_tRNAs)) %>%
  group_by(host_genus) %>%
  group_modify(~{
    df <- .x
    wt <- safe_wilcox(
      df$abs_delta_EnC[df$has_tRNAs == "Without tRNAs"],
      df$abs_delta_EnC[df$has_tRNAs == "With tRNAs"]
    )
    tibble::tibble(
      test = "Wilcoxon",
      W = wt$W,
      p_value = wt$p_value,
      n_without_tRNAs = sum(df$has_tRNAs == "Without tRNAs", na.rm = TRUE),
      n_with_tRNAs = sum(df$has_tRNAs == "With tRNAs", na.rm = TRUE),
      median_without_tRNAs = median(df$abs_delta_EnC[df$has_tRNAs == "Without tRNAs"], na.rm = TRUE),
      median_with_tRNAs = median(df$abs_delta_EnC[df$has_tRNAs == "With tRNAs"], na.rm = TRUE),
      cliffs_delta = cliffs_delta(df$abs_delta_EnC, df$has_tRNAs)
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  trna_presence_tests,
  file.path(out_sum, "tRNA_presence_Wilcoxon_absDeltaEnC.csv")
)

# tRNA count summary
trna_numeric_tests <- per_phage_results %>%
  filter(!is.na(abs_delta_EnC), !is.na(tRNAs)) %>%
  group_by(host_genus) %>%
  group_modify(~{
    safe_spearman(.x$tRNAs, .x$abs_delta_EnC)
  }) %>%
  ungroup() %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  trna_numeric_tests,
  file.path(out_sum, "numeric_tRNA_spearman_absDeltaEnC.csv")
)

# ===================================
# Spearman correlations: EnC vs GC3
# ===================================
spearman_grouped <- function(df, group_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    group_modify(~{
      dd <- .x %>%
        filter(is.finite(GC3), is.finite(EnC))
      
      if (nrow(dd) < 3 || length(unique(dd$GC3)) < 2 || length(unique(dd$EnC)) < 2) {
        tibble(n = nrow(dd), rho = NA_real_, p_value = NA_real_)
      } else {
        ct <- suppressWarnings(cor.test(dd$GC3, dd$EnC, method = "spearman"))
        tibble(n = nrow(dd), rho = unname(ct$estimate), p_value = ct$p.value)
      }
    }) %>%
    ungroup() %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
}

corr_host_vs_phage <- gene_clean %>%
  filter(!is.na(host_genus)) %>%
  spearman_grouped(group_cols = c("host_genus", "sample_type"))
readr::write_csv(
  corr_host_vs_phage,
  file.path(out_sum, "CORR_ALL_host_vs_phage.csv")
)

corr_phage_only <- gene_clean %>%
  filter(sample_type == "phage", !is.na(host_genus)) %>%
  spearman_grouped(group_cols = c("host_genus"))
readr::write_csv(
  corr_phage_only,
  file.path(out_sum, "CORR_PHAGE_only.csv")
)

corr_lifestyle <- gene_clean %>%
  filter(sample_type == "phage", !is.na(host_genus), !is.na(lifestyle)) %>%
  spearman_grouped(group_cols = c("host_genus", "lifestyle"))
readr::write_csv(
  corr_lifestyle,
  file.path(out_sum, "CORR_PHAGE_by_lifestyle.csv")
)

corr_trna_presence <- gene_clean %>%
  filter(sample_type == "phage", !is.na(host_genus), !is.na(has_tRNAs)) %>%
  spearman_grouped(group_cols = c("host_genus", "has_tRNAs"))
readr::write_csv(
  corr_trna_presence,
  file.path(out_sum, "CORR_PHAGE_by_tRNA_presence.csv")
)

corr_trna_numeric <- gene_clean %>%
  filter(
    sample_type == "phage",
    !is.na(host_genus),
    !is.na(tRNAs)
  ) %>%
  mutate(
    tRNAs = as.numeric(tRNAs)
  ) %>%
  spearman_grouped(group_cols = c("host_genus", "tRNAs"))

readr::write_csv(
  corr_trna_numeric,
  file.path(out_sum, "CORR_PHAGE_by_numeric_tRNA_count.csv")
)

# =======
# PLOTS
# =======

# -----------------------
# EnC/GC3 plotting table
# -----------------------
gene_clean <- bind_rows(
  host_genes %>%
    mutate(
      host_genus = NA,
      host_species = NA,
      lifestyle = NA,
      tRNAs = NA_real_,
      has_tRNAs = NA
    ),
  phage_genes
) %>%
  filter(
    is.finite(EnC),
    is.finite(GC3),
    GC3 >= 0,
    GC3 <= 1
  ) %>%
  mutate(
    sample_type = factor(as.character(sample_type), levels = c("host", "phage")),
    host_genus = factor(as.character(host_genus), levels = genus_order_final),
    lifestyle = factor(lifestyle, levels = c("virulent", "temperate")),
    has_tRNAs = factor(has_tRNAs, levels = c("Without tRNAs", "With tRNAs"))
  )

host_key_to_genus <- meta %>%
  distinct(host_key, host_genus)

gene_clean <- gene_clean %>%
  left_join(
    host_key_to_genus,
    by = "host_key",
    suffix = c("", "_from_hostkey")
  ) %>%
  mutate(
    host_genus = dplyr::coalesce(
      as.character(host_genus),
      as.character(host_genus_from_hostkey)
    ),
    host_genus = factor(host_genus, levels = genus_order_final)
  ) %>%
  select(-host_genus_from_hostkey)

readr::write_csv(
  gene_clean,
  file.path(out_dir, "gene_level_EnC_GC3_clean_for_plots.csv")
)

# Expected Nc curve, Wright 1990
nc_expected <- function(gc3) {
  2 + gc3 + 29 / (gc3^2 + (1 - gc3)^2)
}

curve_df <- data.frame(GC3 = seq(0, 1, length.out = 500)) %>%
  mutate(Nc_expected = nc_expected(GC3))

# --------------------------------
# Host vs. phage by genus Nc plot
# --------------------------------
p_hosts_phages_by_genus <- ggplot(gene_clean, aes(x = GC3, y = EnC)) +
  geom_point(
    aes(color = sample_type, shape = sample_type),
    alpha = 0.35,
    size = 1.4
  ) +
  geom_line(
    data = curve_df,
    aes(x = GC3, y = Nc_expected),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(  host = "#F8766D",
                 phage = "#9F79EE"),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(host = 16, phage = 17),
    drop = FALSE
  ) +
  labs(
    x = "GC3",
    y = "EnC / Nc",
    color = "",
    shape = "",
    title = "EnC/GC3 plot: hosts vs phages by host genus"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(8, "pt"),
    panel.grid = element_blank()
  ) +
  facet_wrap(~ host_genus, scales = "free", drop = TRUE)

ggsave(
  file.path(plot_dir, "Ncplot_hosts_vs_phages.png"),
  p_hosts_phages_by_genus,
  width = 12,
  height = 9,
  dpi = 300
)

# -----------------------------
# lifestyle comparison Nc plot
# -----------------------------
phage_only <- gene_clean %>%
  filter(sample_type == "phage")

p_lifestyle <- phage_only %>%
  filter(!is.na(lifestyle)) %>%
  ggplot(aes(x = GC3, y = EnC, color = lifestyle)) +
  geom_point(alpha = 0.45, size = 1.3) +
  geom_line(
    data = curve_df,
    aes(x = GC3, y = Nc_expected),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(
      virulent = "#EE30A7",
      temperate = "#00B0F6"
    ),
    drop = FALSE
  ) +
  labs(
    x = "GC3",
    y = "EnC / Nc",
    color = "Lifestyle",
    title = "Phage EnC/GC3 by lifestyle"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(8, "pt"),
    panel.grid = element_blank()
  ) +
  facet_wrap(~ host_genus, scales = "free", drop = TRUE)

ggsave(
  file.path(plot_dir, "Ncplot_PHAGE_by_lifestyle.png"),
  p_lifestyle,
  width = 12,
  height = 9,
  dpi = 300
)

# ---------------------------------
# tRNA presence comparison Nc plot
# ---------------------------------
p_trna_presence <- phage_only %>%
  filter(!is.na(has_tRNAs)) %>%
  ggplot(aes(x = GC3, y = EnC, color = has_tRNAs)) +
  geom_point(alpha = 0.45, size = 1.3) +
  geom_line(
    data = curve_df,
    aes(x = GC3, y = Nc_expected),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Without tRNAs" = "#999999",
      "With tRNAs" = "#66CDAA"
    ),
    drop = FALSE
  ) +
  labs(
    x = "GC3",
    y = "EnC / Nc",
    color = "tRNA presence",
    title = "Phage EnC/GC3 by tRNA presence"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(8, "pt"),
    panel.grid = element_blank()
  ) +
  facet_wrap(~ host_genus, scales = "free", drop = TRUE)

ggsave(
  file.path(plot_dir, "Ncplot_PHAGE_by_tRNA_presence.png"),
  p_trna_presence,
  width = 12,
  height = 9,
  dpi = 300
)

# ------------------------------
# tRNA count comparison Nc plot
# ------------------------------
p_trna_numeric <- phage_only %>%
  filter(!is.na(tRNAs)) %>%
  ggplot(aes(x = GC3, y = EnC, color = tRNAs)) +
  geom_point(alpha = 0.45, size = 1.3) +
  geom_line(
    data = curve_df,
    aes(x = GC3, y = Nc_expected),
    inherit.aes = FALSE,
    linewidth = 0.8
  ) +
  labs(
    x = "GC3",
    y = "EnC / Nc",
    color = "Number of tRNAs",
    title = "Phage EnC/GC3 by numeric tRNA count"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(8, "pt"),
    panel.grid = element_blank()
  ) +
  facet_wrap(~ host_genus, scales = "free", drop = TRUE)

ggsave(
  file.path(plot_dir, "Ncplot_PHAGE_by_numeric_tRNA_count.png"),
  p_trna_numeric,
  width = 12,
  height = 9,
  dpi = 300
)

# -----------
# ∆EnC plots
# -----------
# lifestyle
p_delta_lifestyle <- per_phage_results %>%
  filter(!is.na(abs_delta_EnC), !is.na(lifestyle)) %>%
  ggplot(aes(x = lifestyle, y = abs_delta_EnC, fill = lifestyle)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, color = "black") +
  facet_wrap(~ host_genus, scales = "free_y") +
  labs(
    x = "",
    y = "|ΔEnC|",
    title = "Absolute ∆EnC by lifestyle"
  ) +
  scale_fill_manual(
    values = c(
      virulent = "#EE30A7",
      temperate = "#00B0F6"
    ),
    drop = FALSE
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  file.path(plot_dir, "absDeltaEnC_violin_lifestyle_by_genus.png"),
  p_delta_lifestyle,
  width = 12,
  height = 8,
  dpi = 300
)

# tRNA presence/absence
p_delta_trna_presence <- per_phage_results %>%
  filter(!is.na(abs_delta_EnC), !is.na(has_tRNAs)) %>%
  ggplot(aes(x = has_tRNAs, y = abs_delta_EnC, fill = has_tRNAs)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.4, color = "black") +
  facet_wrap(~ host_genus, scales = "free_y") +
  labs(
    x = "",
    y = "|ΔEnC|",
    title = "Absolute ∆EnC by tRNA presence"
  ) +
  scale_fill_manual(
    values = c(
      "Without tRNAs" = "#999999",
      "With tRNAs" = "#66CDAA"
    ),
    drop = FALSE
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  file.path(plot_dir, "absDeltaEnC_violin_tRNA_presence_by_genus.png"),
  p_delta_trna_presence,
  width = 12,
  height = 8,
  dpi = 300
)

# tRNA count
p_delta_trna_numeric <- per_phage_results %>%
  filter(!is.na(abs_delta_EnC), !is.na(tRNAs)) %>%
  ggplot(aes(x = tRNAs, y = abs_delta_EnC)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  facet_wrap(~ host_genus, scales = "free") +
  labs(
    x = "Number of phage-encoded tRNAs",
    y = "|ΔEnC|",
    title = "Absolute ∆EnC by numeric tRNA count"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggsave(
  file.path(plot_dir, "absDeltaEnC_scatter_numeric_tRNA_count_by_genus.png"),
  p_delta_trna_numeric,
  width = 12,
  height = 8,
  dpi = 300
)

# ======
# DONE!
# ======
msg("All done! :) Outputs saved in: ", out_dir)
