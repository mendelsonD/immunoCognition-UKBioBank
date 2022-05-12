library(readr)
library(glue)
library(stats)
library(dplyr)

# Define files ----
path <- "/home/doodlefish/Documents/Research/LepageLab/immunologyAndSz/Analysis/immunoCognition-new/Outputs/Summary/Med_C/Raw" # path to raw mediation output summary files
fileNames <- c("UKBB_med_C_results_CRPWeights_05_05_2022","UKBB_excludingDx_med_C_results_05_06_2022","UKBB_excludingSSRI_med_C_results_05_06_2022","UKBB_oldestTert_med_C_results_05_06_2022","UKBB_onlyDx_med_C_results_05_06_2022","UKBB_onlyMedUsers_med_C_results_05_06_2022","UKBB_youngestTert_med_C_results_05_06_2022", "UKBB_NoSSRINoDx_med_C_results_05_09_2022", "UKBB_NoMedNoDx_med_C_results_05_08_2022", "UKBB_noMed_med_C_results_05_08_2022") # name of summary mediation outcome files
# fileNames <- c("UKBB_med_C_results_CRPWeights_05_05_2022")

# correlation analysis --------------
## goals: Determine significant CRP - cog var associations from mediation analysis output
## a) Make df with results from direct effect of CRP on cognition vars
## b) FDR correction on these associations

for(file in fileNames){
  fileName <- glue("{path}/{file}.csv")
  outputDF <- read_csv(fileName, show_col_types = F)
  variableNameCols <- c("X", "Y")
  colOfInt_names <- c("dirEff_df2","dirEff_b","dirEff_SE","dirEff_X_t","dirEff_X_p", "dirEff_R^2adj","dirEff_F","dirEff_df1","dirEff_p")
  colOfInt_nums <- which(colnames(outputDF) %in% c(colOfInt_names,variableNameCols))
  trimmedDF <- outputDF[colOfInt_nums]
  trimmedDF <- unique(trimmedDF)
  trimmedDF <- trimmedDF[order(trimmedDF$dirEff_X_p),]
  p_cor <- p.adjust(trimmedDF$dirEff_X_p, method = "fdr")
  trimmedDF <- trimmedDF %>%   
    mutate("X_p_cor" = p_cor, .after = "dirEff_X_p")
  
  saveName <- glue("{file}_cogSummary.csv")
  write.csv(trimmedDF, file = saveName)
  cat("\n Summary of associations between crp_log_z and cognition variables with corrected p-values saved as: ", saveName)
}

# Mediation analysis --------------
## Goals: 
## a) compute p-values for indirect effects
## b) FDR correction for indirect effects oof cognition variables that are retained

keyCognition <- c("cog_fluidIntel_score_t2_z", "cog_matrix_cor_t2_z", "cog_numMem_maxDigitRemem_t2_z") # name of cognition variables of interest for mediation analyses

outputDf_p_cor <- c()

for(file in fileNames){
  fileName <- glue("{path}/{file}.csv")
  outputDF <- read_csv(fileName, show_col_types = F)
  variableNameCols <- c("modelName", "X", "Y", "M")
  colOfInt_names <- c("IndEff_est","IndEff_SE","IndEff_p","IndEff_95%CI-Lo","IndEff_95%CI-Hi","bPath_model_R^2adj","bPath_model_df1","dirEff_df2","bPath_model_p")
  colOfInt_nums <- which(colnames(outputDF) %in% c(colOfInt_names,variableNameCols))
  outputDF <- outputDF[colOfInt_nums]
  
  z <- outputDF$IndEff_est/outputDF$IndEff_SE
  p <- pnorm(z) # compute p-value
  outputDF <- outputDF %>%
    mutate("IndEff_p" = p, .after = "IndEff_SE")
    
  for(cog in keyCognition){
    outputDF_keyCogOnly <- outputDF %>% filter(Y %in% cog)
    
    p_cor <- p.adjust(outputDF_keyCogOnly$IndEff_p, method = "fdr") # correct p-value
    outputDF_keyCogOnly <- outputDF_keyCogOnly %>%   
      mutate("IndEff_p_cor" = p_cor, .after = "IndEff_p")
    outputDF_keyCogOnly <- outputDF_keyCogOnly[order(outputDF_keyCogOnly$IndEff_p_cor),]   # order by corrected p-value
    outputDF_keyCogOnly_sig <- outputDF_keyCogOnly %>% filter(IndEff_p_cor < .05)
    
    cat("\n File: ", file, "\n Variable: ", cog, ". Significant mediators:")
    if(length(outputDF_keyCogOnly_sig$M) == 0){
      smallest_NonSigP <- round(min(outputDF_keyCogOnly$IndEff_p_cor),3)
      cat("\n\t None. All corrected p > ", smallest_NonSigP, "\n")
    } else{
      cat("(num. = ", length(outputDF_keyCogOnly_sig$M), ") ")
      for(sigMediator in outputDF_keyCogOnly_sig$M){
        cat("\n\t",sigMediator)
      }
      outputDF_keyCogOnly_nonSig <- outputDF_keyCogOnly %>% filter(IndEff_p_cor > .05)
      smallest_NonSigP <- round(min(outputDF_keyCogOnly_nonSig$IndEff_p_cor),3)
      cat("\n Remaining indirect effects, corrected p > ", smallest_NonSigP)
    }
    cat("\n---------------------------")
   
    outputDf_p_cor <- rbind(outputDf_p_cor, outputDF_keyCogOnly_sig)
  }
  
  saveName <- glue("{file}_sig_pCor.csv")
  # write.csv(outputDf_p_cor, file = saveName)
  outputDf_p_cor <- c()
}


View(outputDf_p_cor)

for(file in fileNames){
  fileName <- glue("{path}/{file}.csv")
  outputDF <- read_csv(fileName, show_col_types = F)
  variableNameCols <- c("X", "Y", "M")
  colOfInt_names <- c("IndEff_est","IndEff_SE","IndEff_p","IndEff_p_cor","IndEff_95%CI-Lo","IndEff_95%CI-Hi")
  colOfInt_nums <- which(colnames(outputDF) %in% c(colOfInt_names,variableNameCols))
  trimmedDF_filtered <- trimmedDF %>% filter(IndEff_p_cor < .1)
  
  # 
  # trimmedDF_keyCogOnly_filtered <- trimmedDF_keyCogOnly %>% filter(IndEff_p_cor < .1)
  
  # trimmedDF <- trimmedDF[order(trimmedDF$dirEff_X_p),]
  # p_cor <- p.adjust(trimmedDF$dirEff_X_p, method = "fdr")
  # trimmedDF <- trimmedDF %>%   
  #   mutate("X_p_cor" = p_cor, .after = "dirEff_X_p")
  # 
  # saveName <- glue("{file}_cogSummary.csv")
  # write.csv(trimmedDF, file = saveName)
  # cat("\n Summary of associations between crp_log_z and cognition variables with corrected p-values saved as: ", saveName)
}

table(trimmedDF_keyCogOnly_filtered$)
