### Script to load in the file names for the EDI zipped files
# 11 May 2021, A Hounshell

wd <- getwd()
setwd(wd)

# Load potential libraries
pacman::p_load(tidyverse)

# Load in file names
files <- list.files("./EDI_2021/EEMs_pfiles",full.names=TRUE)
files

strings <- files %>% str_replace(".*/p", "p")

strings <- as.data.frame(strings)

write_csv(strings,"./EDI_2021/EEMs_pfiles_readme_2.csv")

# Load Abs file names
abs_files <- list.files("./EDI_2021/CDOM Correction",full.names=TRUE)

abs_strings <- abs_files %>% str_replace(".*/20","20")

abs_strings <- as.data.frame(abs_strings)

write_csv(abs_strings,"./EDI_2021/CDOM_Correction_readme.csv")
