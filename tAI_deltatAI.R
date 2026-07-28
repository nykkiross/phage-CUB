# ===========================================
# tAI AND ∆tAI ANALYSIS
# WITH LOCAL tRNAscan-SE HOST PREPROCESSING
# ===========================================

# ----------------
# Settings/inputs
# ----------------
metadata_xlsx <- "/your/file/path/sample_phages.xlsx"

out_dir <- "/your/directory/path/Results/RSCU_deltaRSCU"

# Folder where downloaded CDS FASTA files will be cached
cds_cache_dir <- file.path(out_dir, "downloaded_CDS_fastas")

# Set this to TRUE if you want the script to run tRNAscan-SE
run_host_trnascan <- TRUE
# If R cannot find tRNAscan-SE automatically, paste the full path here
trnascan_exe <- Sys.which("tRNAscan-SE")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cds_cache_dir, recursive = TRUE, showWarnings = FALSE)

trnascan_work_dir <- file.path(path.expand("~"), "trnascan_work")

host_fasta_dir <- file.path(trnascan_work_dir, "host_genome_fastas")
host_trnascan_dir <- file.path(trnascan_work_dir, "host_tRNAscan_outputs")

dir.create(trnascan_work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(host_fasta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(host_trnascan_dir, recursive = TRUE, showWarnings = FALSE)

plot_dir <- file.path(out_dir, "plots")
summary_dir <- file.path(out_dir, "summaries")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
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

lifestyle_order <- c("virulent", "temperate")
trna_presence_order <- c("Without tRNAs", "With tRNAs")

eps <- 1e-6
min_weight <- 0.01

# =========
# PACKAGES
# =========
req_pkgs <- c(
  "readxl",
  "dplyr",
  "readr",
  "tidyr",
  "purrr",
  "stringr",
  "forcats",
  "seqinr",
  "rentrez",
  "ggplot2",
  "broom.mixed",
  "emmeans",
  "rstatix",
  "tibble",
  "openxlsx"
)

not_installed <- req_pkgs[!req_pkgs %in% installed.packages()[, "Package"]]

if (length(not_installed) > 0) {
  install.packages(not_installed, dependencies = TRUE)
}

invisible(lapply(req_pkgs, library, character.only = TRUE))

theme_set(theme_minimal(base_size = 14))

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

clean_path <- function(x) {
  x <- as.character(x)
  x <- stringr::str_squish(x)
  x <- dplyr::na_if(x, "")
  x <- dplyr::na_if(x, "NA")
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

safe_wilcox_tbl <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) {
    return(tibble::tibble(W = NA_real_, p_value = NA_real_))
  }
  wt <- try(stats::wilcox.test(x, y, exact = FALSE), silent = TRUE)
  if (inherits(wt, "try-error")) {
    return(tibble::tibble(W = NA_real_, p_value = NA_real_))
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
fetch_ncbi_fasta <- function(
    accession,
    rettype = "fasta",
    out_dir,
    suffix = ".fasta"
) {
  accession <- clean_accession(accession)
  if (is.na(accession)) {
    return(NA_character_)
  }
  safe_acc <- clean_name_for_file(accession)
  out_file <- file.path(out_dir, paste0(safe_acc, suffix))
  if (file.exists(out_file) && file.info(out_file)$size > 0) {
    first_line <- tryCatch(readLines(out_file, n = 1, warn = FALSE), error = function(e) "")
    if (length(first_line) > 0 && stringr::str_starts(first_line[1], ">")) {
      return(out_file)
    } else {
      warning("Cached file is not valid FASTA, deleting: ", out_file)
      file.remove(out_file)
    }
  }
  msg("Fetching ", rettype, " from NCBI: ", accession)
  fasta_txt <- tryCatch(
    rentrez::entrez_fetch(
      db = "nuccore",
      id = accession,
      rettype = rettype,
      retmode = "text"
    ),
    error = function(e) NA_character_
  )
  if (
    is.na(fasta_txt) ||
    !nzchar(fasta_txt) ||
    !stringr::str_starts(strsplit(fasta_txt, "\n")[[1]][1], ">")
  ) {
    msg("Direct fetch failed; trying Entrez search for: ", accession)
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
      msg("Resolved ", accession, " to NCBI UID: ", uid)
      fasta_txt <- tryCatch(
        rentrez::entrez_fetch(
          db = "nuccore",
          id = uid,
          rettype = rettype,
          retmode = "text"
        ),
        error = function(e) NA_character_
      )
    }
  }
  if (is.na(fasta_txt) || !nzchar(fasta_txt)) {
    warning("Empty FASTA returned for accession: ", accession)
    return(NA_character_)
  }
  first_line <- strsplit(fasta_txt, "\n")[[1]][1]
  if (!stringr::str_starts(first_line, ">")) {
    warning("NCBI did not return valid FASTA for accession: ", accession)
    return(NA_character_)
  }
  write_ok <- tryCatch(
    {
      writeLines(fasta_txt, out_file)
      TRUE
    },
    error = function(e) {
      warning(
        "Could not write FASTA for accession: ",
        accession,
        " to file: ",
        out_file,
        "\nError: ",
        conditionMessage(e)
      )
      FALSE
    }
  )
  if (!write_ok) {
    return(NA_character_)
  }
  out_file
}

fetch_cds_fasta <- function(accession, cache_dir = cds_cache_dir) {
  fetch_ncbi_fasta(
    accession = accession,
    rettype = "fasta_cds_na",
    out_dir = cache_dir,
    suffix = "_CDS.fasta"
  )
}

fetch_genome_fasta <- function(accession, fasta_dir = host_fasta_dir) {
  fetch_ncbi_fasta(
    accession = accession,
    rettype = "fasta",
    out_dir = fasta_dir,
    suffix = ".fna"
  )
}

# =========================
# FASTA and RSCU functions
# =========================
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

calculate_gene_rscu_from_fasta <- function(
    fasta_file,
    sample_id,
    accession,
    sample_type,
    host_key,
    chromosome,
    min_length_nt = 200
) {
  seqs <- safe_read_fasta(fasta_file)
  if (length(seqs) == 0) {
    warning("No sequences found/readable for accession: ", accession)
    return(tibble::tibble())
  }
  seqs <- seqs[vapply(seqs, is_valid_cds, logical(1), min_length_nt = min_length_nt)]
  if (length(seqs) == 0) {
    warning("All CDS filtered out for accession: ", accession)
    return(tibble::tibble())
  }
  rscu_list <- lapply(seqs, function(seq) {
    gene_name <- attr(seq, "name")
    codon_counts <- seqinr::uco(seq, index = "rscu")
    tibble::tibble(
      sample_id = sample_id,
      accession = accession,
      sample_type = sample_type,
      host_key = host_key,
      chromosome = chromosome,
      gene_id = gene_name,
      codon = tolower(names(codon_counts)),
      RSCU = as.numeric(codon_counts)
    )
  })
  dplyr::bind_rows(rscu_list)
}

# ==============================================
# Local tRNAscan-SE functions and tAI functions
# ==============================================
check_trnascan_available <- function(trnascan_exe) {
  if (is.na(trnascan_exe) || trnascan_exe == "") {
    stop(
      "Could not find tRNAscan-SE on PATH. ",
      "Activate your conda environment first or set trnascan_exe to the full path from `which tRNAscan-SE`."
    )
  }
  test <- tryCatch(
    system2(trnascan_exe, args = "-h", stdout = TRUE, stderr = TRUE),
    error = function(e) NA_character_
  )
  if (all(is.na(test))) {
    stop("tRNAscan-SE was found but could not be executed: ", trnascan_exe)
  }
  invisible(TRUE)
}

run_trnascan_one <- function(
    fasta_file,
    accession,
    output_dir = host_trnascan_dir,
    trnascan_exe_path = trnascan_exe
) {
  accession <- clean_accession(accession)
  if (is.na(fasta_file) || !file.exists(fasta_file)) {
    warning("Missing FASTA for tRNAscan-SE accession: ", accession)
    return(NA_character_)
  }
  if (is.na(trnascan_exe_path) || trnascan_exe_path == "" || !file.exists(trnascan_exe_path)) {
    warning("tRNAscan-SE executable not found for accession: ", accession)
    return(NA_character_)
  }
  safe_acc <- clean_name_for_file(accession)
  out_txt <- file.path(output_dir, paste0(safe_acc, "_tRNAscan.txt"))
  out_struct <- file.path(output_dir, paste0(safe_acc, "_tRNAscan.struct"))
  if (file.exists(out_txt) && file.info(out_txt)$size > 0) {
    return(out_txt)
  }
  unlink(c(out_txt, out_struct), force = TRUE)
  
  msg("Running tRNAscan-SE for: ", accession)
  
  fasta_file <- normalizePath(fasta_file, mustWork = TRUE)
  out_txt <- normalizePath(out_txt, mustWork = FALSE)
  out_struct <- normalizePath(out_struct, mustWork = FALSE)
  trnascan_exe_path <- normalizePath(trnascan_exe_path, mustWork = TRUE)
  args <- c(
    "-B",
    "-o", out_txt,
    "-f", out_struct,
    fasta_file
  )
  res_status <- tryCatch(
    system2(
      command = trnascan_exe_path,
      args = args,
      stdout = FALSE,
      stderr = FALSE
    ),
    error = function(e) {
      warning(
        "system2 failed for tRNAscan-SE accession ",
        accession,
        ": ",
        conditionMessage(e)
      )
      return(1)
    }
  )
  if (is.null(res_status)) {
    res_status <- 0
  }
  if (res_status != 0) {
    warning("tRNAscan-SE returned non-zero exit status for accession: ", accession)
    return(NA_character_)
  }
  if (!file.exists(out_txt) || file.info(out_txt)$size == 0) {
    warning("tRNAscan-SE produced empty/missing output for accession: ", accession)
    return(NA_character_)
  }
  out_txt
}

empty_trna_tbl <- function() {
  tibble::tibble(
    Name = character(),
    tRNA_Type = character(),
    Anti_Codon = character(),
    Score = numeric(),
    host_key = character(),
    source_label = character()
  )
}

process_trnascan_txt <- function(
    file_path,
    host_key,
    source_label = NA_character_,
    optional = FALSE
) {
  file_path <- clean_path(file_path)
  if (is.na(file_path) || !file.exists(file_path)) {
    if (!optional) {
      warning("Missing REQUIRED tRNAscan file for host_key ", host_key, ": ", file_path)
    }
    return(empty_trna_tbl())
  }
  df <- tryCatch(
    readr::read_tsv(
      file_path,
      skip = 2,
      col_names = c(
        "Name",
        "tRNA_number",
        "Begin",
        "End",
        "tRNA_Type",
        "Anti_Codon",
        "Intron_Begin",
        "Intron_End",
        "Score",
        "Isotype_CM",
        "Isotype_Score",
        "Note"
      ),
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    ),
    error = function(e) {
      warning("Could not read tRNAscan file: ", file_path, "\n", conditionMessage(e))
      return(empty_trna_tbl())
    }
  )
  if (nrow(df) == 0) {
    return(empty_trna_tbl())
  }
  df %>%
    dplyr::select(Name, tRNA_Type, Anti_Codon, Score) %>%
    dplyr::filter(
      !is.na(Name),
      !is.na(Anti_Codon),
      !is.na(Score),
      !stringr::str_detect(Name, "^-+$"),
      !stringr::str_detect(Anti_Codon, "^-+$"),
      !stringr::str_detect(Score, "^-+$"),
      !stringr::str_detect(Name, "^Name$"),
      !stringr::str_detect(Name, "^Sequence$")
    ) %>%
    dplyr::mutate(
      host_key = host_key,
      source_label = source_label,
      Anti_Codon = tolower(Anti_Codon),
      Score = suppressWarnings(as.numeric(Score))
    ) %>%
    dplyr::filter(
      !is.na(Score),
      !is.na(Anti_Codon),
      Anti_Codon != "",
      Anti_Codon != "-"
    )
}

generate_anticodon_dictionary <- function(trna_data) {
  tapply(as.numeric(trna_data$Score), trna_data$Anti_Codon, sum, na.rm = TRUE)
}

reverse_complement <- function(codon) {
  comp <- c("a" = "t", "t" = "a", "g" = "c", "c" = "g")
  paste0(rev(comp[strsplit(codon, NULL)[[1]]]), collapse = "")
}

calculate_tai_with_rscu <- function(rscu_data, anticodon_dict, min_weight = 0.01) {
  rscu_data <- rscu_data %>%
    dplyr::mutate(
      tAI_weight = mapply(
        function(codon, rscu_value) {
          anticodon <- reverse_complement(codon)
          matching_scores <- anticodon_dict[anticodon]
          
          if (!is.na(matching_scores)) {
            (rscu_value + 1e-6) * pmax(matching_scores, min_weight)
          } else {
            min_weight
          }
        },
        codon,
        RSCU,
        SIMPLIFY = TRUE
      )
    )
  gene_tai <- rscu_data %>%
    dplyr::group_by(
      sample_id,
      accession,
      sample_type,
      host_key,
      chromosome,
      gene_id
    ) %>%
    dplyr::summarise(
      raw_tAI = exp(mean(log(pmax(tAI_weight, min_weight)), na.rm = TRUE)),
      .groups = "drop"
    )
  max_tai <- max(gene_tai$raw_tAI, na.rm = TRUE)
  gene_tai %>%
    dplyr::mutate(tAI = raw_tAI / max_tai)
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

if (!("Host_tRNAscan" %in% names(meta_raw))) {
  meta_raw$Host_tRNAscan <- NA_character_
}

if (!("Host_tRNAscan2" %in% names(meta_raw))) {
  meta_raw$Host_tRNAscan2 <- NA_character_
}

meta <- meta_raw %>%
  dplyr::mutate(
    dplyr::across(
      .cols = where(is.character) &
        !any_of(c(
          "Accession",
          "Host_Accession",
          "Host_Accession2",
          "Host_tRNAscan",
          "Host_tRNAscan2"
        )),
      .fns = standardize_text
    ),
    
    phage_accession = clean_accession(Accession),
    host_accession = clean_accession(Host_Accession),
    host_accession2 = clean_accession(Host_Accession2),
    
    host_trnascan = clean_path(Host_tRNAscan),
    host_trnascan2 = clean_path(Host_tRNAscan2),
    
    phage_id = as.character(Phage_ID),
    sample_name = phage_id,
    
    host_genus = as.character(Host_Genus),
    host_species = as.character(Host_Species),
    
    host_key = dplyr::case_when(
      !is.na(host_accession2) ~ paste(host_accession, host_accession2, sep = " + "),
      TRUE ~ host_accession
    ),
    
    lifestyle = tolower(as.character(Lifestyle)),
    lifestyle = dplyr::case_when(
      lifestyle %in% c("temperate", "temp") ~ "temperate",
      lifestyle %in% c("virulent", "vir") ~ "virulent",
      TRUE ~ lifestyle
    ),
    
    tRNAs = as.numeric(tRNAs),
    
    has_tRNAs = dplyr::case_when(
      is.na(tRNAs) ~ NA_character_,
      tRNAs > 0 ~ "With tRNAs",
      tRNAs == 0 ~ "Without tRNAs"
    ),
    
    has_tRNAs = factor(has_tRNAs, levels = trna_presence_order),
    
    host_genus = factor(host_genus, levels = genus_order),
    lifestyle = factor(lifestyle, levels = lifestyle_order)
  ) %>%
  dplyr::filter(
    !is.na(phage_accession),
    !is.na(host_accession),
    !is.na(host_genus)
  )

meta <- meta %>%
  dplyr::mutate(
    host_genus = factor(as.character(host_genus), levels = genus_order_final)
  )

readr::write_csv(
  meta,
  file.path(out_dir, "metadata_cleaned_accession_based.csv")
)

# ------------------------------------------
# Build organism table for tAI calculation
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

# ============================================================
# Run host tRNAscan-SE preprocessing
# ============================================================
trnascan_exe <- "/your/path/miniforge3/envs/trnascan/bin/tRNAscan-SE"
check_trnascan_available(trnascan_exe)

file.exists(trnascan_exe)

host_accessions_chr1 <- meta %>%
  dplyr::distinct(host_key, host_accession) %>%
  dplyr::transmute(
    host_key = host_key,
    accession = host_accession,
    chromosome = "chromosome_1"
  )

host_accessions_chr2 <- meta %>%
  dplyr::filter(!is.na(host_accession2)) %>%
  dplyr::distinct(host_key, host_accession2) %>%
  dplyr::transmute(
    host_key = host_key,
    accession = host_accession2,
    chromosome = "chromosome_2"
  )

host_accessions_for_trnascan <- dplyr::bind_rows(
  host_accessions_chr1,
  host_accessions_chr2
) %>%
  dplyr::distinct(host_key, accession, chromosome) %>%
  dplyr::filter(!is.na(accession))

readr::write_csv(
  host_accessions_for_trnascan,
  file.path(out_dir, "host_accessions_for_tRNAscan.csv")
)

one_test <- host_accessions_for_trnascan %>%
  dplyr::slice(1) %>%
  dplyr::mutate(
    genome_fasta = purrr::map_chr(accession, fetch_genome_fasta),
    trnascan_txt = purrr::map2_chr(
      genome_fasta,
      accession,
      ~ run_trnascan_one(
        fasta_file = .x,
        accession = .y,
        output_dir = host_trnascan_dir,
        trnascan_exe_path = trnascan_exe
      )
    )
  )

print(one_test)

if (run_host_trnascan) {
  check_trnascan_available(trnascan_exe)
  msg("Fetching host genome FASTAs and running tRNAscan-SE...")
  host_trnascan_results <- host_accessions_for_trnascan %>%
    dplyr::mutate(
      genome_fasta = purrr::map_chr(accession, fetch_genome_fasta),
      trnascan_txt = purrr::map2_chr(
        genome_fasta,
        accession,
        ~ run_trnascan_one(
          fasta_file = .x,
          accession = .y,
          output_dir = host_trnascan_dir,
          trnascan_exe_path = trnascan_exe
        )
      )
    )
  readr::write_csv(
    host_trnascan_results,
    file.path(out_dir, "host_tRNAscan_run_results_long.csv")
  )
  host_trnascan_paths <- host_trnascan_results %>%
    dplyr::select(host_key, chromosome, trnascan_txt) %>%
    tidyr::pivot_wider(
      names_from = chromosome,
      values_from = trnascan_txt
    ) %>%
    dplyr::rename(
      host_trnascan = chromosome_1,
      host_trnascan2 = chromosome_2
    )
  readr::write_csv(
    host_trnascan_paths,
    file.path(out_dir, "host_tRNAscan_paths_by_host_key.csv")
  )
  meta <- meta %>%
    dplyr::select(-any_of(c("host_trnascan", "host_trnascan2"))) %>%
    dplyr::left_join(host_trnascan_paths, by = "host_key")
} else {
  msg("Skipping local tRNAscan-SE run. Using Host_tRNAscan columns from metadata.")
}

# ===============================================
# Fetch CDS FASTAs and calculate gene-level RSCU
# ===============================================
msg("Fetching CDS FASTAs for ", nrow(organisms), " unique organism records...")

organisms_fetched <- organisms %>%
  dplyr::mutate(
    cds_fasta = purrr::map_chr(accession, fetch_cds_fasta)
  )

readr::write_csv(
  organisms_fetched,
  file.path(out_dir, "unique_organisms_with_CDS_paths.csv")
)

msg("Calculating gene-level RSCU...")

all_gene_rscu <- organisms_fetched %>%
  dplyr::mutate(
    rscu_data = purrr::pmap(
      list(cds_fasta, sample_id, accession, sample_type, host_key, chromosome),
      calculate_gene_rscu_from_fasta
    )
  ) %>%
  dplyr::select(rscu_data) %>%
  tidyr::unnest(rscu_data)

readr::write_csv(
  all_gene_rscu,
  file.path(out_dir, "rscu_tidy_per_gene.csv")
)

# --------------------------------------------------------
# Attach the RSCU data for phage and host to the metadata
# --------------------------------------------------------
phage_gene_rscu <- all_gene_rscu %>%
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

host_gene_rscu <- all_gene_rscu %>%
  dplyr::filter(sample_type == "host")

readr::write_csv(
  phage_gene_rscu,
  file.path(out_dir, "phage_gene_RSCU_with_metadata.csv")
)

readr::write_csv(
  host_gene_rscu,
  file.path(out_dir, "host_gene_RSCU.csv")
)

# ===================================
# Generate tRNAscan-SE dictionaries
# ===================================
host_trnascan_files <- meta %>%
  dplyr::distinct(
    host_key,
    host_accession,
    host_accession2,
    host_trnascan,
    host_trnascan2
  )

missing_trnascan <- host_trnascan_files %>%
  dplyr::filter(is.na(host_trnascan))

readr::write_csv(
  missing_trnascan,
  file.path(out_dir, "missing_host_tRNAscan_files.csv")
)

if (nrow(missing_trnascan) > 0) {
  warning(
    "Some host_key values are missing Host_tRNAscan files. ",
    "Those host backgrounds cannot be used for tAI. ",
    "See missing_host_tRNAscan_files.csv"
  )
}

host_trna_data <- host_trnascan_files %>%
  dplyr::mutate(
    trna_chr1 = purrr::map2(
      host_trnascan,
      host_key,
      ~ process_trnascan_txt(
        file_path = .x,
        host_key = .y,
        source_label = "chromosome_1",
        optional = FALSE
      )
    ),
    # chromosome 2 is optional
    trna_chr2 = purrr::map2(
      host_trnascan2,
      host_key,
      ~ process_trnascan_txt(
        file_path = .x,
        host_key = .y,
        source_label = "chromosome_2",
        optional = TRUE
      )
    )
  ) %>%
  dplyr::transmute(
    host_key,
    trna_data = purrr::map2(trna_chr1, trna_chr2, dplyr::bind_rows)
  )

host_dict_tbl <- host_trna_data %>%
  dplyr::mutate(
    anticodon_dict = purrr::map(trna_data, generate_anticodon_dictionary),
    n_tRNAscan_rows = purrr::map_int(trna_data, nrow)
  )

host_dict_tbl %>%
  dplyr::select(host_key, n_tRNAscan_rows) %>%
  print(n = 100)

# ===================================
# Calculate tAI by matched host_key
# ===================================
msg("Calculating tAI by matched host_key...")

all_gene_tai_list <- list()

for (hk in unique(meta$host_key)) {
  msg("Processing host_key: ", hk)
  dict_row <- host_dict_tbl %>%
    dplyr::filter(host_key == hk)
  if (nrow(dict_row) == 0 || dict_row$n_tRNAscan_rows[[1]] == 0) {
    warning("Skipping host_key with no tRNAscan dictionary: ", hk)
    next
  }
  anticodon_dict <- dict_row$anticodon_dict[[1]]
  subset_rscu <- dplyr::bind_rows(
    host_gene_rscu %>% dplyr::filter(host_key == hk),
    phage_gene_rscu %>% dplyr::filter(host_key == hk)
  )
  if (nrow(subset_rscu) == 0) {
    warning("No RSCU rows for host_key: ", hk)
    next
  }
  gene_tai <- calculate_tai_with_rscu(
    subset_rscu,
    anticodon_dict,
    min_weight = min_weight
  )
  all_gene_tai_list[[hk]] <- gene_tai
}

gene_tai_results_raw <- dplyr::bind_rows(all_gene_tai_list)

gene_tai_results <- gene_tai_results_raw %>%
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
  ) %>%
  dplyr::mutate(
    host_genus = dplyr::coalesce(
      as.character(host_genus),
      as.character(meta$host_genus[match(host_key, meta$host_key)])
    ),
    host_species = dplyr::coalesce(
      as.character(host_species),
      as.character(meta$host_species[match(host_key, meta$host_key)])
    ),
    host_genus = factor(host_genus, levels = genus_order_final),
    lifestyle = factor(lifestyle, levels = lifestyle_order),
    has_tRNAs = factor(has_tRNAs, levels = trna_presence_order)
  )

readr::write_csv(
  gene_tai_results,
  file.path(out_dir, "tai_results_per_gene_all_hosts.csv")
)

# ---------------------------
# Genome-level tAI and ∆tAI
# ---------------------------
genome_tai <- gene_tai_results %>%
  dplyr::group_by(
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
  dplyr::summarise(
    mean_raw_tAI = mean(raw_tAI, na.rm = TRUE),
    median_raw_tAI = median(raw_tAI, na.rm = TRUE),
    mean_tAI = mean(tAI, na.rm = TRUE),
    median_tAI = median(tAI, na.rm = TRUE),
    n_genes = dplyr::n_distinct(gene_id),
    .groups = "drop"
  )

host_tai_by_key <- genome_tai %>%
  dplyr::filter(sample_type == "host") %>%
  dplyr::group_by(host_key) %>%
  dplyr::summarise(
    host_mean_tAI = mean(mean_tAI, na.rm = TRUE),
    host_median_tAI = median(median_tAI, na.rm = TRUE),
    host_mean_raw_tAI = mean(mean_raw_tAI, na.rm = TRUE),
    .groups = "drop"
  )

tai_results <- genome_tai %>%
  dplyr::left_join(host_tai_by_key, by = "host_key") %>%
  dplyr::mutate(
    delta_tAI = dplyr::if_else(
      sample_type == "phage",
      mean_tAI - host_mean_tAI,
      NA_real_
    ),
    abs_delta_tAI = abs(delta_tAI)
  )

readr::write_csv(
  tai_results,
  file.path(out_dir, "tai_results_genome_level_all_hosts.csv")
)

host_gene_tai_by_key <- gene_tai_results %>%
  dplyr::filter(sample_type == "host") %>%
  dplyr::group_by(host_key) %>%
  dplyr::summarise(
    host_gene_mean_tAI = mean(tAI, na.rm = TRUE),
    .groups = "drop"
  )

gene_tai_results <- gene_tai_results %>%
  dplyr::left_join(host_gene_tai_by_key, by = "host_key") %>%
  dplyr::mutate(
    delta_tAI = dplyr::if_else(
      sample_type == "phage",
      tAI - host_gene_mean_tAI,
      NA_real_
    ),
    abs_delta_tAI = abs(delta_tAI)
  )

readr::write_csv(
  gene_tai_results,
  file.path(out_dir, "tai_results_per_gene_all_hosts_with_delta.csv")
)

# ================================
# Summaries and statistical tests
# ================================

# ---------------------------
# Summaries of tAI and ∆tAI
# ---------------------------
host_phage_summary <- tai_results %>%
  dplyr::group_by(host_genus, sample_type) %>%
  dplyr::summarise(
    n = dplyr::n_distinct(sample_id),
    median_tAI = median(mean_tAI, na.rm = TRUE),
    mean_tAI = mean(mean_tAI, na.rm = TRUE),
    sd_tAI = sd(mean_tAI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    id_cols = host_genus,
    names_from = sample_type,
    values_from = c(n, median_tAI, mean_tAI, sd_tAI),
    names_sep = "_"
  )

readr::write_csv(
  host_phage_summary,
  file.path(summary_dir, "summary_tAI_host_phage.csv")
)

delta_summary <- tai_results %>%
  dplyr::filter(sample_type == "phage") %>%
  dplyr::group_by(host_genus) %>%
  dplyr::summarise(
    n_phages = dplyr::n_distinct(sample_id),
    mean_tAI = mean(mean_tAI, na.rm = TRUE),
    median_tAI = median(mean_tAI, na.rm = TRUE),
    mean_delta_tAI = mean(delta_tAI, na.rm = TRUE),
    median_delta_tAI = median(delta_tAI, na.rm = TRUE),
    sd_delta_tAI = sd(delta_tAI, na.rm = TRUE),
    mean_abs_delta_tAI = mean(abs_delta_tAI, na.rm = TRUE),
    median_abs_delta_tAI = median(abs_delta_tAI, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  delta_summary,
  file.path(summary_dir, "summary_delta_tAI.csv")
)

# ------------------------------------
# Genome-level ∆tAI statistical tests
# ------------------------------------
phage_tai <- tai_results %>%
  dplyr::filter(sample_type == "phage", !is.na(delta_tAI), !is.na(host_genus))

stats_delta <- phage_tai %>%
  dplyr::group_by(host_genus) %>%
  rstatix::wilcox_test(delta_tAI ~ 1, mu = 0, alternative = "two.sided") %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::left_join(
    phage_tai %>%
      dplyr::group_by(host_genus) %>%
      dplyr::summarise(
        n = dplyr::n(),
        median_delta = median(delta_tAI, na.rm = TRUE),
        iqr_delta = IQR(delta_tAI, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "host_genus"
  )

readr::write_csv(
  stats_delta,
  file.path(summary_dir, "stats_delta_tAI_phage_vs_zero.csv")
)

# lifestyle comparison
lifestyle_tests <- phage_tai %>%
  dplyr::filter(!is.na(lifestyle), !is.na(abs_delta_tAI)) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::group_modify(~{
    df <- .x
    wt <- safe_wilcox_tbl(
      df$abs_delta_tAI[df$lifestyle == "virulent"],
      df$abs_delta_tAI[df$lifestyle == "temperate"]
    )
    tibble::tibble(
      test = "Wilcoxon",
      W = wt$W,
      p_value = wt$p_value,
      n_virulent = sum(df$lifestyle == "virulent", na.rm = TRUE),
      n_temperate = sum(df$lifestyle == "temperate", na.rm = TRUE),
      median_virulent = median(df$abs_delta_tAI[df$lifestyle == "virulent"], na.rm = TRUE),
      median_temperate = median(df$abs_delta_tAI[df$lifestyle == "temperate"], na.rm = TRUE),
      cliffs_delta = cliffs_delta(df$abs_delta_tAI, df$lifestyle)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  lifestyle_tests,
  file.path(summary_dir, "stats_lifestyle_absDeltaTAI.csv")
)

# tRNA presence comparison
trna_presence_tests <- phage_tai %>%
  dplyr::filter(!is.na(has_tRNAs), !is.na(abs_delta_tAI)) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::group_modify(~{
    df <- .x
    wt <- safe_wilcox_tbl(
      df$abs_delta_tAI[df$has_tRNAs == "Without tRNAs"],
      df$abs_delta_tAI[df$has_tRNAs == "With tRNAs"]
    )
    tibble::tibble(
      test = "Wilcoxon",
      W = wt$W,
      p_value = wt$p_value,
      n_without_tRNAs = sum(df$has_tRNAs == "Without tRNAs", na.rm = TRUE),
      n_with_tRNAs = sum(df$has_tRNAs == "With tRNAs", na.rm = TRUE),
      median_without_tRNAs = median(df$abs_delta_tAI[df$has_tRNAs == "Without tRNAs"], na.rm = TRUE),
      median_with_tRNAs = median(df$abs_delta_tAI[df$has_tRNAs == "With tRNAs"], na.rm = TRUE),
      cliffs_delta = cliffs_delta(df$abs_delta_tAI, df$has_tRNAs)
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  trna_presence_tests,
  file.path(summary_dir, "stats_tRNA_presence_absDeltaTAI.csv")
)

# tRNA count correlations
trna_numeric_tests <- phage_tai %>%
  dplyr::filter(!is.na(tRNAs), !is.na(abs_delta_tAI)) %>%
  dplyr::group_by(host_genus) %>%
  dplyr::group_modify(~{
    safe_spearman(.x$tRNAs, .x$abs_delta_tAI)
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

readr::write_csv(
  trna_numeric_tests,
  file.path(summary_dir, "stats_numeric_tRNA_absDeltaTAI.csv")
)

# ======
# PLOTS
# ======
theme_large <- theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.grid = element_blank()
  )

clean_labels <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(x)
}

# phage vs. host
p_host_vs_phage <- ggplot(tai_results, aes(x = host_genus, y = mean_tAI)) +
  geom_violin(
    data = dplyr::filter(tai_results, sample_type == "phage"),
    aes(fill = "Phages"),
    trim = FALSE
  ) +
  geom_boxplot(
    data = dplyr::filter(tai_results, sample_type == "phage"),
    width = 0.15,
    alpha = 0.5,
    outlier.shape = NA
  ) +
  geom_point(
    data = dplyr::filter(tai_results, sample_type == "host"),
    aes(color = "Host"),
    size = 3,
    shape = 18
  ) +
  scale_fill_manual(values = c("Phages" = "steelblue")) +
  scale_color_manual(values = c("Host" = "red")) +
  labs(
    title = "Host vs Phage tAI by Host Genus",
    y = "tAI",
    x = "Host Genus",
    fill = NULL,
    color = NULL
  ) +
  theme_large

ggsave(
  file.path(plot_dir, "violin_tAI_host_vs_phage.png"),
  p_host_vs_phage,
  width = 11,
  height = 7,
  dpi = 300
)

# lifestyle plot ∆tAI
p_lifestyle_delta <- phage_tai %>%
  dplyr::filter(!is.na(lifestyle)) %>%
  ggplot(aes(x = lifestyle, y = delta_tAI, fill = lifestyle)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55) +
  facet_wrap(~ host_genus, scales = "free_y") +
  scale_fill_manual(values = palette_lifestyle, name = "Lifestyle") +
  scale_x_discrete(labels = clean_labels) +
  labs(
    title = "ΔtAI by Lifestyle",
    x = "Lifestyle",
    y = "ΔtAI"
  ) +
  theme_large

ggsave(
  file.path(plot_dir, "violin_delta_tAI_by_lifestyle.png"),
  p_lifestyle_delta,
  width = 12,
  height = 8,
  dpi = 300
)

# tRNA presence plot ∆tAI
p_trna_presence_delta <- phage_tai %>%
  dplyr::filter(!is.na(has_tRNAs)) %>%
  ggplot(aes(x = has_tRNAs, y = delta_tAI, fill = has_tRNAs)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.55) +
  facet_wrap(~ host_genus, scales = "free_y") +
  scale_fill_manual(values = palette_trna_presence, name = "tRNA presence") +
  labs(
    title = "ΔtAI by tRNA Presence",
    x = "",
    y = "ΔtAI"
  ) +
  theme_large

ggsave(
  file.path(plot_dir, "violin_delta_tAI_by_tRNA_presence.png"),
  p_trna_presence_delta,
  width = 12,
  height = 8,
  dpi = 300
)

# tRNA numeric plot ∆tAI
p_trna_numeric_delta <- phage_tai %>%
  dplyr::filter(!is.na(tRNAs)) %>%
  ggplot(aes(x = tRNAs, y = delta_tAI)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  facet_wrap(~ host_genus, scales = "free") +
  labs(
    title = "ΔtAI by Numeric tRNA Count",
    x = "Number of phage-encoded tRNAs",
    y = "ΔtAI"
  ) +
  theme_large

ggsave(
  file.path(plot_dir, "scatter_delta_tAI_by_numeric_tRNA_count.png"),
  p_trna_numeric_delta,
  width = 12,
  height = 8,
  dpi = 300
)

# tRNA numeric plot tAI
p_trna_numeric_tai <- phage_tai %>%
  dplyr::filter(!is.na(tRNAs)) %>%
  ggplot(aes(x = tRNAs, y = mean_tAI)) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  facet_wrap(~ host_genus, scales = "free") +
  labs(
    title = "Phage tAI by Numeric tRNA Count",
    x = "Number of phage-encoded tRNAs",
    y = "Mean tAI"
  ) +
  theme_large

ggsave(
  file.path(plot_dir, "scatter_tAI_by_numeric_tRNA_count.png"),
  p_trna_numeric_tai,
  width = 12,
  height = 8,
  dpi = 300
)

# ======
# DONE!
# ======
msg("All done! :) Outputs saved in: ", out_dir)
