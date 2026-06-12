# This script integrates called TCRs
# For J2177 sequencing
# Clonotypes are tracked by their unique TCR Beta CDR3 amino acid sequence

##### Libraries to load
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(writexl)
library(readxl)
library(stringr)

##### A. Data Import #####
# This sets the working directory. Please change as needed
setwd()

# Setting the data input directory
Dirwd <- getwd()
DirInput <- paste0(Dirwd, "/Input/")
DirOutput <- paste0(Dirwd, "/Output/")

# Listing of the .csv files in the folder and getting their full path
Dirtsv_file_names <- list.files(path = DirInput, pattern = ".csv$")
Dirtsv_files <- paste0(DirInput, Dirtsv_file_names)

#Gets the metadata file in the Data_Bulk folder
meta = read_excel(paste0(DirInput, "metadata.xlsx"))

# Now inputs all of the .csv files into a list of data frames called BulkTCR
BulkTCR <- lapply(Dirtsv_files, function(file) {
  read.csv(file, header = TRUE, check.names = FALSE)
})
names(BulkTCR) <- gsub("\\.csv$", "", Dirtsv_file_names)

# Renames each sample based on the metadata file
for (i in seq_along(BulkTCR)) {
  current_name <- names(BulkTCR)[i] 
  if (current_name %in% meta$`File Name`) {
    matching_row <- meta[meta$`File Name` == current_name, ]
    new_name <- matching_row$`File Name` 
    names(BulkTCR)[i] <- new_name
  }
  rm(current_name)
}
rm(matching_row)

#### B. Data Cleanup ####
# Creating columns for Patient and Antigen
for (df_name in names(BulkTCR)) {
  df <- BulkTCR[[df_name]]
  
  meta_row <- meta[meta$`File Name` == df_name, ]
  
  df$Patient <- meta_row$Patient
  
  if (meta_row$Patient != "Published") {
    df$Antigen <- meta_row$Mutation
  }
  
  if (meta_row$Patient == "Published") {
    df$De_novo <- TRUE
    df$Prime <- FALSE
    df$Persistent <- FALSE
  }
  
  BulkTCR[[df_name]] <- df
  
  BulkTCR[[df_name]] <- BulkTCR[[df_name]] %>%
    mutate(Unique = paste(CDR3aa, TRBV, TRBJ, sep = "_"))  %>%
    select(CDR3aa, TRBV, TRBJ, Patient, Antigen, De_novo, Unique, Prime, Persistent)
  
  rm(df, meta_row)
}

#### C. Data Merging ####
# Creates a single table with all of the TCRs from each patient
TCR_Merged <- do.call(rbind, BulkTCR)
rm(BulkTCR)

# Merge TCR clonotypes from the same patient
TCR_Merged <- TCR_Merged %>%
  group_by(Unique, Patient, CDR3aa, TRBV, TRBJ, Prime, Persistent) %>%
  summarise(
    De_novo = any(De_novo, na.rm = TRUE),
    Antigen = paste(unique(na.omit(Antigen)), collapse = ", "),
    .groups = "drop"
  )

# Split into rows with and without allele information
TCR_missing_V <- TCR_Merged[
  is.na(TCR_Merged$TRBV) | TCR_Merged$TRBV == "" |
    is.na(TCR_Merged$TRBJ) | TCR_Merged$TRBJ == "", 
]

TCR_with_V <- TCR_Merged[!rownames(TCR_Merged) %in% rownames(TCR_missing_V), ]

#### D. Public Clonotypes ####
# 1. Identify public clonotypes via identical Unique
Public <- TCR_with_V %>%
  group_by(Unique) %>%
  filter(n_distinct(Patient) > 1) %>%
  ungroup()

# 2. Now, for clonotypes without a TCRV, match by CDR3aa only
Public_noV <- map_dfr(seq_len(nrow(TCR_missing_V)), function(i) {
  row_i <- TCR_missing_V[i, ]
  
  matches <- TCR_with_V %>%
    filter(
      CDR3aa == row_i$CDR3aa,
      Patient != row_i$Patient
    )
  
  if (nrow(matches) > 0) {
    bind_rows(row_i, matches)
  } else {
    NULL
  }
})

# Merge both types of public clonotypes and assign each Public clonotype as a pair
TCR_Public <- bind_rows(Public, Public_noV)
rm(Public)
rm(Public_noV)
rm(TCR_missing_V)
rm(TCR_with_V)

TCR_Public <- TCR_Public %>%
  group_by(CDR3aa) %>%
  mutate(Public_ID = cur_group_id()) %>%
  ungroup() %>%
  relocate('Public_ID', .after = 'TRBJ')

# Identify the public clonotypes in the merged clonotype dataframe
TCR_Merged_annotated <- TCR_Merged %>%
  left_join(
    TCR_Public %>% select(Unique, Patient, Public_ID),
    by = c("Unique", "Patient")
  )

# Filter for Public TCR rows and collapse them into one row
# Identifies all patients that have the TCR
# Identifies all antigens that it is reactive to across patients
# For frequencies, identifies the maximum frequency & count across each row
TCR_Public_collapsed <- TCR_Merged_annotated %>%
  filter(!is.na(Public_ID)) %>%
  group_by(Public_ID) %>%
  mutate(
    Public = TRUE,
    Patient = toString(unique(Patient)),
    Antigen = toString(sort(unique(unlist(strsplit(Antigen, ",\\s*"))))),
    De_novo = any(De_novo, na.rm = TRUE)
  ) %>%
  filter(!is.na(TRBV) & TRBV != "") %>%
  filter(!is.na(TRBJ) & TRBJ != "") %>%
  slice_head(n = 1) %>%
  ungroup()
rm(TCR_Merged_annotated)

# Remove all Public clonotypes (by Unique) from the original TCR_Merged
TCR_Merged <- TCR_Merged %>%
  filter(!Unique %in% TCR_Public$Unique)

# Add in the collapsed Public clonotype rows
TCR_Merged <- bind_rows(TCR_Merged, TCR_Public_collapsed)

# Cleaning up
TCR_Merged <- TCR_Merged %>%
  mutate(
    Public = ifelse(is.na(Public), FALSE, Public)
  ) %>%
  select(-Public_ID)  # Remove the Unique_ID column

rm(TCR_Public_collapsed)

#### F. Exporting Data ####
write.csv(TCR_Public, paste0(DirOutput, "All_TCR_Public.csv"), row.names=FALSE)
write.csv(TCR_Merged, paste0(DirOutput, "All_TCR_Merged.csv"), row.names=FALSE)