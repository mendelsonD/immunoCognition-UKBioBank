library(readr)
library(glue)
library(stats)

# Define files ----
path <- "/home/doodlefish/Documents/Research/LepageLab/immunologyAndSz/Analysis/immunoCognition-new/Outputs/Summary/Med_C"
fileNames <- c("UKBB_med_C_results_CRPWeights_05_05_2022","UKBB_excludingDx_med_C_results_05_06_2022","UKBB_excludingSSRI_med_C_results_05_06_2022","UKBB_oldestTert_med_C_results_05_06_2022","UKBB_onlyDx_med_C_results_05_06_2022","UKBB_onlyMedUsers_med_C_results_05_06_2022","UKBB_youngestTert_med_C_results_05_06_2022")
# fileNames <- c("UKBB_med_C_results_CRPWeights_05_05_2022") # for testing purposes

# correlation analysis --------------
## goals: Determine significant CRP - cog var associations from mediation analysis output
## a) Make df with results from direct effect of CRP on cognition vars
## b) FDR correction on these associations

for(file in fileNames){
  fileName <- glue("{path}/{file}.csv")
  outputDF <- read_csv(fileName, show_col_types = F)
  variableNameCols <- c("X", "Y")
  colOfInt_names <- c("dirEff_df2","dirEff_b","dirEff_SE","dirEff_X_t","dirEff_X_p")
  colOfInt_nums <- which(colnames(outputDF) %in% c(colOfInt_names,variableNameCols))
  trimmedDF <- outputDF[colOfInt_nums]
  trimmedDF <- unique(trimmedDF)
  trimmedDF <- trimmedDF[order(trimmedDF$dirEff_X_p),]
  p_cor <- p.adjust(trimmedDF$dirEff_X_p, method = "fdr")
  trimmedDF$p_cor <- p_cor
  
  saveName <- glue("{file}_cogSummary.csv")
  write.csv(trimmedDF, file = saveName)
  cat("\n Summary of associations between crp_log_z and cognition variables with corrected p-values saved as: ", saveName)
}

# Mediation analysis --------------
## Goals: 
## a) compute p-values for indirect effects
## b) FDR correction for indirect effects

for(file in fileNames){
  fileName <- glue("{path}/{file}.csv")
  outputDF <- read_csv(fileName, show_col_types = F)
  
  z <- outputDF$IndEff_est/outputDF$IndEff_SE
  p <- pnorm(z) # compute p-value
  outputDF <- outputDF %>% 
    mutate("IndEff_p" = p, .after = "IndEff_SE")
    
  p_cor <- p.adjust(outputDF$IndEff_p, method = "fdr") # correct p-value
  outputDF <- outputDF %>%   
    mutate("IndEff_p_cor" = p_cor, .after = "IndEff_p")
  
  # order by corrected p-value
  outputDF <- outputDF[order(outputDF$IndEff_p_cor),]
  saveName <- glue("{file}_pCor.csv")
  write.csv(outputDF, file = saveName)
}
