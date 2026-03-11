# This script loads mKRAS TCR data and constructs Sankey diagrams
# Last edited by Henry Wang on 09.08.25

# Loading Libraries
library(writexl)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)

##### 0. Variable Selection ####
# Sets the patient ID for analysis
Patient_ID <- '01'

##### A. Loading Data #####
# Set the working directory. Change as needed
setwd("C:/Users/henry/Dropbox/Research/Jaffee Lab/J2177/Sequencing")

# Setting input and output directories
Dirwd <- getwd()
DirMeta <- paste0(Dirwd, "/Data_Bulk/") #Location of metadata file
DirDataBulk <- paste0(Dirwd, "/Data_Bulk/", Patient_ID, "/") #Location of input sequencing
DirOutput <- paste0(Dirwd, "/Output/", Patient_ID, "/") #Output folder

# Loads metadata
meta_path <- file.path(DirMeta, "metadata.xlsx")
meta <- read_excel(meta_path)
rm(meta_path)

# Identify the mKRAS peptides for the patient
peptides <- meta %>%
  filter(Patient == Patient_ID) %>%
  filter(!Condition %in% c("Unstim", "FLC")) %>%
  pull(Condition) %>%
  unique()
rm(meta)

# Reads the merged data frame files
TCR_merged <- setNames(
  lapply(peptides, function(peptide) {
    filename <- file.path(
      DirOutput,
      paste0(Patient_ID, "_", peptide, "_Summary.csv")
    )
    read.csv(filename, header = TRUE, check.names = FALSE)
  }),
  peptides
)

##### B. Reformatting #####
# Each clonotype is defined as Persistent if it expands at 5w or 13w AND
# either is expanded at 1y or 2y otherwise, it is a 'transient' clonotype
BulkTCR_sankey <- lapply(
  TCR_merged,
  function(df) df[ , setdiff(names(df), c("TRBV","TRBJ")), drop = FALSE]
)

sig_cols <- grep("_sig$", names(BulkTCR_sankey[[1]]), value = TRUE)
timepoints <- sub("_sig$", "", sig_cols)

for (peptide in names(BulkTCR_sankey)) {
  TCR_df <- BulkTCR_sankey[[peptide]]
  
  TCR_df$Signature <- if_else(
    TCR_df$Persistent,
    "Persistent",
    "Transient"
  )
  
  for (col in sig_cols) {
    TCR_df[[col]] <- ifelse(
      TCR_df[[col]],
      TCR_df$Signature,
      "Not expanded"
    )
  }
  
  TCR_df <- TCR_df %>%
    mutate(
      `13w_sig` = if_else(
        (`13w_sig` == "Not expanded") & (`5w_sig` == "Persistent"),
        "Not expanded_2",
        `13w_sig`
      )
    )
  
  if ("2y" %in% timepoints) {
    TCR_df <- TCR_df %>%
      mutate(
        `1y_sig` = if_else(
          (`1y_sig` == "Not expanded") & (`5w_sig` == "Persistent"),
          "Not expanded_2",
          `1y_sig`
        )
      )
  }
  
  TCR_df$row_id <- seq_len(nrow(TCR_df))
  
  BulkTCR_sankey[[peptide]] <- TCR_df
  rm(TCR_df, col)
}

##### C. Plotting #####
plots <- list()
plot_width <- ifelse(("2y" %in% timepoints), 9, 6) 

# Plots with weights based on expansion frequency
for (peptide in names(BulkTCR_sankey)) {
  
  tcr_df <- BulkTCR_sankey[[peptide]]
  total_5w <- sum(tcr_df$"5w_freq")/150
  
  # Prepare long format
  long_df <- tcr_df %>%
    pivot_longer(cols = all_of(sig_cols), names_to = "Timepoint", values_to = "State") %>%
    mutate(
      State = factor(State, levels = rev(c("Not expanded", "Transient", "Persistent", "Not expanded_2"))),
      Timepoint = sub("_.*$", "", Timepoint),
      Timepoint = factor(Timepoint, levels = timepoints)
    )
  
  long_df <- long_df %>%
    rowwise() %>%
    mutate(
      Freq = case_when(
        State == "Not expanded_2" ~ 0,
        State == "Not expanded"   ~ -total_5w,
        State %in% c("Persistent", "Transient") ~ 
          100 * cur_data()[[paste0(Timepoint, "_freq")]]
      ),
      Alpha = case_when(
        State == "Not expanded"                     ~ 0,
        State %in% c("Persistent", "Not expanded_2") ~ 1,
        State == "Transient"                  ~ 1,
        TRUE                                  ~ NA_real_
      )
    ) %>%
    ungroup()
  
  # Plot
  p <- ggplot(long_df,
              aes(x = Timepoint, stratum = State, alluvium = row_id, y = Freq)) +
    geom_flow(aes(fill = State, alpha = Alpha), width = 1/8) +
    geom_stratum(aes(fill = State), width = 1/8, color = "black") +
    coord_cartesian(clip = "off") +
    scale_x_discrete(expand = c(0.05, 0)) +
    scale_y_continuous(labels = abs, expand = c(0, 0)) +
    scale_alpha_identity() +
    scale_fill_manual(
      values = c(
        "Persistent" = "#1f77b4",
        "Transient" = "lightgray",
        "Not expanded" = "#ff7f0e",
        "Not expanded_2" = "#1f77b4"
      ),
      breaks = c("Persistent", "Transient", "Not expanded"),
      name = NULL
    ) +
    theme_minimal() +
    theme(
      axis.title.y = element_text(size = 14),
      axis.title.x = element_text(size = 14), 
      axis.text.y = element_text(size = 12),
      axis.ticks.y = element_line(),
      axis.line.y = element_line(),
      axis.ticks.x = element_line(),
      axis.line.x = element_line(),
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    labs(x = "Timepoint",
      y = "Post-expansion Frequency (%)",
      title = paste("Freq. Clonotypes Plot - Patient", Patient_ID, " - ", peptide))
  
  # Save the plot
  ggsave(
    filename = file.path(DirOutput, paste0(Patient_ID, "_", peptide, "_freq_clonotypes.tiff")),
    plot = p,
    width = plot_width,
    height = 6
  )
  
  rm(long_df, p, tcr_df, total_5w)
}