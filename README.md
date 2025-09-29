# J2177-NCT05013216
R code for analysis of bulk TCRseq data for cohort A of NCT05013216. Instructions for running code are below.

A. Data structure
1. TCR sequencing files should be placed in /Data_Bulk/ within individual subfolders according to patient ID (two digit IDs). File names should match the provided metadata.xlsx file
2. Code output will be in /Output/ folder with individualized patient output in subfolders according to patient ID.

B. Peripheral TCRseq analysis
1. 'J2177_BulkTCR_Peripheral_Analysis' is a R script with code for annotating mKRAS-specific TCRs and plotting clonotype expansion scatter plots.
2. Code will need to be rerun for each patient based on patient ID.

C. Sankey diagrams
1. 'J2177_Sankey' is a R script with code for generating Sankey diagrams tracking persistent mKRAS-specific TCRs. 'J2177_BulkTCR_Peripheral_Analysis' must have been run first.
2. Code will need to be rerun for each patient based on patient ID.
