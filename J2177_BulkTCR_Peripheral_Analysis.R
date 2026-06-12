# This script loads Adaptive Bulk TCRSeq Data to identify mKRAS-reactive T-cell clonotypes
# Clonotypes are tracked by their unique TCRBeta CDR3 amino acid sequence
# Last edited by Henry Wang on 09.12.25

# Loading Libraries
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(writexl)
library(readxl)

##### 0. Variable Selection ####
# Sets the patient ID for analysis
Patient_ID <- '08'

# Sets FDR for mKRAS TCRs
FDR <- 0.01

# Sets min. fold-exp for mKRAS TCRs
Ratio <- 5

# Sets min. post-exp freq for mKRAS TCRs
exp_freq <- 2.5e-4

##### A. Loading Data #####
# Set the working directory. Change as needed
setwd()

# Setting input and output directories
Dirwd <- getwd()
DirMeta <- paste0(Dirwd, "/Data_Bulk/") #Location of metadata file
DirDataBulk <- paste0(Dirwd, "/Data_Bulk/", Patient_ID, "/") #Location of sequencing data
DirOutput <- paste0(Dirwd, "/Output/", Patient_ID, "/") #Output folder

# Retrieves file paths for all .tsv files in /Data_Bulk
Dirtsv_file_names <- list.files(path = DirDataBulk, pattern = ".tsv$")
Dirtsv_files <- paste0(DirDataBulk, Dirtsv_file_names)

# Loads metadata
meta_path <- file.path(DirMeta, "metadata.xlsx")
meta <- read_excel(meta_path)
rm(meta_path)

# Read the .tsv files into a list of data frames - BulkTCR
BulkTCR <- lapply(Dirtsv_files, function(file) {
  read.delim(file, header = TRUE, sep = "\t")
})
names(BulkTCR) <- gsub("\\.tsv$", "", Dirtsv_file_names)
rm(Dirtsv_file_names, Dirtsv_files)

# Renames each file based on metadata
for (i in seq_along(BulkTCR)) {
  current_name <- names(BulkTCR)[i] 
  
  if (current_name %in% meta$`File Name`) {
    matching_row <- meta[meta$`File Name` == current_name, ]
    new_name <- paste0(matching_row$Timepoint, "_", matching_row$Condition)
    names(BulkTCR)[i] <- new_name
  }
  rm(matching_row, new_name, current_name)
}

#### B. Data Clean-up ####
# Only keep columns that have an in frame CDR3aa (productive TCRs)
columns_to_keep <- c("aminoAcid",
                     "count..templates.reads.",
                     "vMaxResolved",
                     "jMaxResolved")
new_column_names <- c("CDR3aa",
                      "Counts",
                      "vResolved",
                      "jResolved")

# Only keep columns that have an in frame CDR3aa (productive TCRs)
# CDR3aa must also begin with C and end in F or W
BulkTCR <- lapply(BulkTCR, function(df) {
  df <- df[df$sequenceStatus == "In" &
             grepl("^C.*[FW]$", df$aminoAcid), ]
  return(df)
})

# Adjust columns in BulkTCR
BulkTCR <- lapply(BulkTCR, function(df) {
  df <- df %>% select(all_of(columns_to_keep))
  colnames(df) <- new_column_names
  return(df)
})
rm(columns_to_keep, new_column_names)

# Create a unique identifier for each clonotype based on the CDR3aa, TCRV, TCRJ
add_unique_column <- function(df) {
  df$Unique <- paste(df$CDR3aa, df$vResolved, df$jResolved, sep = "_")
  return(df)
}
BulkTCR <- lapply(BulkTCR, add_unique_column)

# Merge rows with same unique identifier, but sum the counts
BulkTCR_unique <- lapply(BulkTCR, function(df) {
  df %>%
    group_by(Unique) %>%
    summarise(
      Counts = sum(Counts, na.rm = TRUE),
      Clonotypes = n(),
      across(-Counts, first),
      .groups = "drop"
    )
})

# Create a new summary data frame 'meta_lim'
meta_lim <- data.frame(names(BulkTCR))
colnames(meta_lim) <- 'Sample'
rm(BulkTCR)

# Renames the data frames in BulkTCR_unique
names(BulkTCR_unique) <- meta_lim$Sample

# Save # productive counts in meta_lim file
for(i in seq_along(BulkTCR_unique)){
  Total <- sum(BulkTCR_unique[[i]]$Counts)
  meta_lim$Total_Counts[[i]] <- Total
  meta_lim$Unique_Clonotypes[[i]] <- nrow(BulkTCR_unique[[i]])
  rm(Total)
}
meta_lim$Total_Counts <- as.numeric(as.character(meta_lim$Total_Counts))
meta_lim$Unique_Clonotypes <- as.numeric(as.character(meta_lim$Unique_Clonotypes))

# For log-log plotting purposes, create a pseudo-frequency count
meta_lim$Unobs_Freq <- (3* meta_lim$Total_Counts)^(-1)

# Calculates clonotype frequencies based on Total_Counts
BulkTCR_unique <- lapply(names(BulkTCR_unique), function(sample_name) {
  df <- BulkTCR_unique[[sample_name]]
  total_counts <- meta_lim$Total_Counts[meta_lim$Sample == sample_name]
  df$Frequency <- df$Counts / total_counts
  df
})
names(BulkTCR_unique) <- meta_lim$Sample

#### C. Identifying mKRAS TCRs ####
#### C.1 Merging Datasets ####
# Merging the datasets for each peptide against FLC control
# For each row, replace 0 with NA under counts
# Note that Counts.x represents the mKRAS peptide exp whereas Count.y represents the FLC exp

# First, determine the peptides used for expansion
peptide_names <- meta %>%
  filter(Condition != "Unstim") %>%
  filter(Patient == Patient_ID) %>%
  pull(Condition) %>%
  unique()

peptide_names_no_FLC <- peptide_names[peptide_names != "FLC"]
Baseline <- BulkTCR_unique[["Baseline_Unstim"]]
BulkTCR_unique$Baseline_Unstim <- NULL

# Creating new merged data frames
BulkTCR_merged <- list()

# Determine the timepoints
timepoints <- meta %>%
  filter(Timepoint != "Baseline") %>%
  filter(Patient == Patient_ID) %>%
  pull(Timepoint) %>%
  unique()

# Create merged data frames for each timepoint & peptide pair
for (tp in timepoints) {
  flc_name <- paste0(tp, "_FLC")
  
  flc_df <- BulkTCR_unique[[flc_name]]
  
  for (pep in peptide_names_no_FLC) {
    peptide_name <- paste0(tp, "_", pep)
    peptide_df <- BulkTCR_unique[[peptide_name]]
    
    merged_df <- full_join(
      peptide_df, flc_df,
      by = c("CDR3aa", "vResolved", "jResolved", "Unique")
    ) %>%
      mutate(
        Counts.x = ifelse(is.na(Counts.x), 0, Counts.x),
        Counts.y = ifelse(is.na(Counts.y), 0, Counts.y)
      )
    
    merged_name <- paste0(tp, "_", pep)
    BulkTCR_merged[[merged_name]] <- merged_df
  }
  rm(flc_name,
     pep,
     peptide_name,
     peptide_df,
     flc_df,
     merged_name,
     merged_df,
     tp)
}

# Setting unobserved pseudo-frequencies
for (df_name in names(BulkTCR_merged)) {
  
  timepoint <- strsplit(df_name, "_")[[1]][1]
  flc_name <- paste0(timepoint, "_FLC")
  
  freq_x <- meta_lim$Unobs_Freq[meta_lim$Sample == df_name]
  freq_y <- meta_lim$Unobs_Freq[meta_lim$Sample == flc_name]
  
  df <- BulkTCR_merged[[df_name]]
  
  df$Frequency.x <- ifelse(is.na(df$Frequency.x), freq_x, df$Frequency.x)
  df$Frequency.y <- ifelse(is.na(df$Frequency.y), freq_y, df$Frequency.y)
  
  BulkTCR_merged[[df_name]] <- df
  rm(df,
     df_name,
     timepoint,
     flc_name,
     freq_x,
     freq_y)
}

#### C.2 Statistically Enriched ####
# For each mKRAS peptide exp, conduct Fisher's exact test against the FLC control expansion
# Then, calculate p adjusted using Benjamini-Yekutieli's correction for FDR

for (j in seq_along(BulkTCR_merged)) {
  df <- BulkTCR_merged[[j]]
  sample_name <- names(BulkTCR_merged)[j]
  timepoint <- strsplit(sample_name, "_")[[1]][1]
  flc_name <- paste0(timepoint, "_FLC")
  total_peptide <- meta_lim$Total_Counts[meta_lim$Sample == sample_name]
  total_FLC <- meta_lim$Total_Counts[meta_lim$Sample == flc_name]
  
  contig_tables <- mapply(
    function(count_x, count_y) {
      matrix(c(
        count_x, total_peptide - count_x,
        count_y, total_FLC - count_y
      ), nrow = 2)
    },
    df$Counts.x, df$Counts.y, SIMPLIFY = FALSE
  )
  
  p_values <- sapply(contig_tables, function(tbl) fisher.test(tbl, alternative = "greater")$p.value)
  
  df$P_value <- p_values
  df$P_value_BY <- p.adjust(p_values, method = "BY")
  BulkTCR_merged[[j]] <- df
  
  rm(contig_tables,
     df,
     sample_name,
     timepoint,
     flc_name,
     total_peptide,
     total_FLC,
     p_values)
}

#### C.3 FDR & Ratio Filter ####
# Filter mKRAS TCRs based on previous set FDR cutoffs
BulkTCR_merged_sig <- list()
for( i in seq_along(BulkTCR_merged)){
  df <- BulkTCR_merged[[i]] %>% filter(P_value_BY < FDR)
  df <- df %>% filter(Frequency.x / Frequency.y > Ratio)
  BulkTCR_merged_sig[[names(BulkTCR_merged)[i]]] <- df
  rm(df)
}

#### C.4 Min Freq Filter ####
# Filter mKRAS TCRs based on post-exp freq
BulkTCR_merged_sig_final <- list()
for( i in seq_along(BulkTCR_merged_sig)){
  BulkTCR_merged_sig_final[[i]] <- BulkTCR_merged_sig[[i]] %>% filter(Frequency.x >= exp_freq)
}
names(BulkTCR_merged_sig_final) <- names(BulkTCR_merged_sig)

# For the expanded clonotypes, only keep ones not present at baseline
BulkTCR_merged_sig_final <- lapply(BulkTCR_merged_sig_final,function(df) {
    df %>% filter(!Unique %in% Baseline$Unique)
  }
)

# Now identify in the original merged list data frames whether a clonotype is expanded
add_significance <- function(df_list, input_list, col_name) {
  Map(function(name, df) {
    input_df <- input_list[[name]]
    unique_ids <- if (!is.null(input_df)) input_df$Unique else character(0)
    df[[col_name]] <- df$Unique %in% unique_ids
    df
  }, names(df_list), df_list)
}

BulkTCR_merged <- add_significance(BulkTCR_merged, BulkTCR_merged_sig_final, "Sig")

#### D. Pairwise Plots####
# Plot expanded clonotypes
for (sample_name in names(BulkTCR_merged)) {
  plot_data <- BulkTCR_merged[[sample_name]] %>%
    mutate(group = ifelse(Sig, "highlight", "background")) %>%
    group_by(Frequency.x, Frequency.y, group) %>%
    summarise(n = n(), .groups = "drop")
  
  accent_col <- if (grepl("G12V", sample_name)) "#FB8808"
  else if (grepl("G12D", sample_name)) "#5770FF"
  else "orange"
  
  x_cutoff <- plot_data %>%
    distinct(Frequency.y) %>%
    arrange(Frequency.y) %>%
    slice_head(n = 2) %>%
    summarise(geom_mean = exp(mean(log(Frequency.y)))) %>%
    pull(geom_mean)
  
  y_cutoff <- plot_data %>%
    distinct(Frequency.x) %>%
    arrange(Frequency.x) %>%
    slice_head(n = 2) %>%
    summarise(geom_mean = exp(mean(log(Frequency.x)))) %>%
    pull(geom_mean)   
  
  p <- ggplot(plot_data, aes(Frequency.y, Frequency.x, size = n, color = group)) +
    geom_point(alpha = 0.5, show.legend = FALSE) +
    scale_color_manual(values = c(background = "black", highlight = accent_col)) +
    geom_vline(xintercept = x_cutoff, color = "black", linetype = "dashed", linewidth = 1) +
    geom_hline(yintercept = y_cutoff, color = "black", linetype = "dashed", linewidth = 1) +
    
    scale_size_continuous(
      trans = "log2",
      name = "# Clonotypes",
      breaks = c(1, 10, 100, 1000),
      labels = c("1", "10", "100", "1000"),
      range = c(2,6)) +
    
    labs(x = "FLC Frequency",
      y = "Peptide Frequency",
      title = paste("Expanded clonotypes: ", sample_name)) +
    
    theme(legend.position = "right", 
          legend.key = element_rect(fill = "white"),
          aspect.ratio = 1) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    coord_equal()
  
  ggsave(
    filename = file.path(DirOutput, paste0(sample_name, "_clonotype_plot.tiff")),
    plot = p,
    device = "tiff", width = 8, height = 6
  )
  rm(p, plot_data, sample_name)
}

#### E. Additional Stats ####
# All additional statistics are stored in meta_lim table for export

#### E.1 Total Frequencies ####
# Tracks total mKRAS clonotype frequencies across conditions
Frequency <- matrix(nrow = length(names(BulkTCR_merged_sig_final)), ncol = 2)
rownames(Frequency) <- names(BulkTCR_merged_sig_final)
colnames(Frequency) <- c('FLC_Exp_Freq', 'mKRAS_Exp_Freq')

# Including the baseline PBMC frequencies in BulkTCR_merged_sig_final
for (df_name in names(BulkTCR_merged_sig_final)) {
  df <- BulkTCR_merged_sig_final[[df_name]]
  matched_indices <- match(df$Unique, Baseline$Unique)
  BulkTCR_merged_sig_final[[df_name]] <- df
}
rm(df, df_name, matched_indices)

# Now calculating cumulative clonotype frequencies
for (df_name in names(BulkTCR_merged_sig_final)) {
  df <- BulkTCR_merged_sig_final[[df_name]]
  
  mKRAS_Exp_Freq  <- sum(df$Frequency.x, na.rm = TRUE)
  FLC_Exp_Freq    <- sum(df$Frequency.y[df$Counts.y != 0], na.rm = TRUE)
  
  Frequency[df_name, "mKRAS_Exp_Freq"]   <- mKRAS_Exp_Freq
  Frequency[df_name, "FLC_Exp_Freq"]     <- FLC_Exp_Freq
  
  rm(df_name, df)
}

#### E.2 Number of Clonotypes ####
Frequency <- cbind(Frequency, Clonotypes = NA)

for (df_name in names(BulkTCR_merged_sig_final)) {
  Frequency[df_name, "Clonotypes"] <- nrow(BulkTCR_merged_sig_final[[df_name]])
  rm(df_name)
}

#### E.3 Merge into meta_lim ####
for (col in colnames(Frequency)) {
  matched_idx <- match(meta_lim$Sample, rownames(Frequency))
  meta_lim[[col]] <- Frequency[matched_idx, col]
  rm(col)
}

#### E.4 Clonality ####
# Create a function for calculating clonality from BulkTCR_unique
calculate_clonality <- function(df) {
  freq <- df$Counts / sum(df$Counts)
  entropy <- -sum(freq * log(freq))
  clonality <- 1 - (entropy / log(sum(df$Counts)))
  return(clonality)
}

# Calculate clonality for all conditions and store in meta_lim
clonality_calc <- tibble(
  Sample    = names(BulkTCR_unique),
  Clonality = map_dbl(BulkTCR_unique, calculate_clonality)
)

meta_lim <- meta_lim %>%
  left_join(clonality_calc, by = "Sample") %>%
  relocate(Clonality, .after = Unobs_Freq)
  
meta_lim$Clonality[ meta_lim$Sample == "Baseline_Unstim" ] <- calculate_clonality(Baseline)

rm(clonality_calc, Frequency)

#### E.5 Merge Timepoints ####
# Merge TCRs based on each peptide
# First separate into different data frames
for (peptide in peptide_names_no_FLC) {
  assign(paste0("TCR_", peptide), list())
  assign(paste0("TCR_", peptide), list())
}

for (name in names(BulkTCR_merged_sig_final)) {
  parts <- strsplit(name, "_")[[1]]
  timepoint <- parts[1]
  peptide <- parts[2]
  df_name <- paste0("TCR_", peptide)
  df_TCR <- get(df_name)
  df_TCR[[timepoint]] <- BulkTCR_merged_sig_final[[name]]
  assign(df_name, df_TCR)
  
  rm(parts, timepoint, peptide, df_name, df_TCR)
}

# Clean up columns
clean_tcr_list <- function(tcr_list) {
  lapply(names(tcr_list), function(df_name) {
    df <- tcr_list[[df_name]]
    df <- df[, c("Unique", "CDR3aa", "vResolved", "jResolved", "Frequency.x", "Counts.x", "Counts.y")]
    
    freq_col_name <- paste0(df_name, "_freq")
    names(df)[names(df) == "Frequency.x"] <- freq_col_name
    foldexp_col_name <- paste0(df_name, "_foldexp")
    df[[foldexp_col_name]] <- with(df, ifelse(Counts.y == 0, Inf, Counts.x / Counts.y))
    df <- df[, !names(df) %in% c("Counts.x", "Counts.y")]
    return(df)
  }) |> setNames(names(tcr_list))
}

for (peptide in peptide_names_no_FLC){
  df_name <- paste0("TCR_", peptide)
  list_TCR <- get(df_name)
  list_TCR <- clean_tcr_list(list_TCR)
  assign(df_name, list_TCR)
  rm(df_name, list_TCR)
}

# Now merging peptides into single data frame
merge_tcr_list <- function(tcr_list) {
  list_name <- deparse(substitute(tcr_list))
  merge_cols <- c("Unique", "CDR3aa", "vResolved", "jResolved")
  
  merged_df <- Reduce(function(x, y) {
    full_join(x, y, by = merge_cols)
  }, tcr_list)
  
  merged_name <- paste0(list_name, "_merged")
  assign(merged_name, merged_df, envir = .GlobalEnv)
  
  rm(merged_df, merge_cols, list_name)
}

if ("G12V" %in% peptide_names_no_FLC) {
  merge_tcr_list(TCR_G12V)
  TCR_G12V_merged <- TCR_G12V_merged %>%
    { df <- .
    foldexp_cols <- grep("_foldexp$", names(df), value = TRUE)
    
    for (fc in foldexp_cols) {
      sig_col <- sub("_foldexp$", "_sig", fc)
      df[[sig_col]] <- !is.na(df[[fc]])
      df <- df %>% relocate(all_of(sig_col), .after = all_of(fc))
    }
    df
    }
  rm(TCR_G12V)
}

if ("G12D" %in% peptide_names_no_FLC) {
  merge_tcr_list(TCR_G12D)
  TCR_G12D_merged <- TCR_G12D_merged %>%
    { df <- .
    foldexp_cols <- grep("_foldexp$", names(df), value = TRUE)
    
    for (fc in foldexp_cols) {
      sig_col <- sub("_foldexp$", "_sig", fc)
      df[[sig_col]] <- !is.na(df[[fc]])
      df <- df %>% relocate(all_of(sig_col), .after = all_of(fc))
    }
    df
    }
  rm(TCR_G12D)
}

# Now finds missing expansion data for those not called
fill_data <- function(peptide, timepoints) {
  df_name <- paste0("TCR_", peptide, "_merged")
  tcr_df <- get(df_name, envir = .GlobalEnv)
  
  for (i in seq_len(nrow(tcr_df))) {
    unique_val <- tcr_df$Unique[i]
    
    for (tp in timepoints) {
      freq_col <- paste0(tp, "_freq")
      foldexp_col <- paste0(tp, "_foldexp")
      
      if (is.na(tcr_df[[freq_col]][i])) {
        bulk_df_name <- paste0(tp, "_", peptide)
        bulk_df <- BulkTCR_merged[[bulk_df_name]]
        
        match_row <- bulk_df[bulk_df$Unique == unique_val & bulk_df$Counts.x > 0, ]
        
        if (nrow(match_row) == 0) {
          tcr_df[[freq_col]][i] <- 0
          tcr_df[[foldexp_col]][i] <- 0
        } else {
          tcr_df[[freq_col]][i] <- match_row$Frequency.x[1]
          tcr_df[[foldexp_col]][i] <- ifelse(match_row$Counts.y[1] == 0, Inf, match_row$Counts.x[1] / match_row$Counts.y[1])
        }
      }
    }
  }
  
  assign(df_name, tcr_df, envir = .GlobalEnv)
  invisible(tcr_df)
}

if("G12V" %in% peptide_names_no_FLC){
  fill_data("G12V", timepoints)
  
  if ("2y_sig" %in% names(TCR_G12V_merged)) {
    TCR_G12V_merged <- TCR_G12V_merged %>%
      mutate(
        Prime = `5w_sig` | `13w_sig`,
        Persistent = Prime & (`1y_sig` | `2y_sig`)
      )
  } else {
    TCR_G12V_merged <- TCR_G12V_merged %>%
      mutate(
        Prime = `5w_sig` | `13w_sig`,
        Persistent = Prime & `1y_sig`
      )
  }
}

if("G12D" %in% peptide_names_no_FLC){
  fill_data("G12D", timepoints)
  
  if ("2y_sig" %in% names(TCR_G12D_merged)) {
    TCR_G12D_merged <- TCR_G12D_merged %>%
      mutate(
        Prime = `5w_sig` | `13w_sig`,
        Persistent = Prime & (`1y_sig` | `2y_sig`)
      )
  } else {
    TCR_G12D_merged <- TCR_G12D_merged %>%
      mutate(
        Prime = `5w_sig` | `13w_sig`,
        Persistent = Prime & `1y_sig`
      )
  }
}

#### F. Exporting Data ####
#### F.1 Reformatting ####
if("G12V" %in% peptide_names_no_FLC){
  TCR_G12V_merged <- TCR_G12V_merged %>%
    select(-Unique) %>%
    rename(TRBV = vResolved, TRBJ = jResolved)
}

if("G12D" %in% peptide_names_no_FLC){
  TCR_G12D_merged <- TCR_G12D_merged %>%
    select(-Unique) %>%
    rename(TRBV = vResolved, TRBJ = jResolved)
}

meta_lim <- meta_lim %>%
  separate(Sample, into = c("Timepoint", "Exp"), sep = "_") %>%
  relocate(Timepoint, Exp, .before = 1)

meta_lim <- meta_lim %>%
  mutate(Timepoint = factor(Timepoint, levels = c(timepoints, "Baseline"))) %>%
  arrange(Timepoint)

#### F.2 Final Export ####
write.csv(meta_lim, paste0(DirOutput, Patient_ID, "_TCR_Numbers.csv"), row.names=FALSE)

if("G12V" %in% peptide_names_no_FLC){
  write.csv(TCR_G12V_merged, paste0(DirOutput, Patient_ID, "_G12V_Summary.csv"), row.names=FALSE)
}

if("G12D" %in% peptide_names_no_FLC){
  write.csv(TCR_G12D_merged, paste0(DirOutput, Patient_ID, "_G12D_Summary.csv"), row.names=FALSE)
}
