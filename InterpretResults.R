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

getPMed <- function(path, fileNames, keyCognition, keyMediators, pCor, sigOnly, alpha, outputAffix) {
  outputDf_p_cor <- c()

  for (file in fileNames) {
    fileName <- glue("{path}/{file}.csv")
    outputDF <- read_csv(fileName, show_col_types = F)

    # defines columns of interest
    variableNameCols <- c("modelName", "X", "Y", "M")
    colOfInt_names <- c("IndEff_est", "IndEff_SE", "IndEff_p", "IndEff_95%CI-Lo", "IndEff_95%CI-Hi", "bPath_model_R^2adj", "bPath_model_df1", "dirEff_df2", "bPath_model_p")
    colOfInt_nums <- which(colnames(outputDF) %in% c(colOfInt_names, variableNameCols))
    outputDF <- outputDF[colOfInt_nums]

    z <- outputDF$IndEff_est / outputDF$IndEff_SE
    p <- pnorm(z) # compute p-value
    outputDF <- outputDF %>%
      mutate("IndEff_p" = p, .after = "IndEff_SE")
    
    fileNameAffix <- "_p"

    if (pCor == TRUE) {
      fileNameAffix <- "_pCor"
      for (cog in keyCognition) {
        outputDF_keyYM <- outputDF %>%
          filter(Y %in% cog) %>%
          filter(M %in% keyMediators)

        p_cor <- p.adjust(outputDF_keyYM$IndEff_p, method = "fdr") # correct p-value
        outputDF_keyYM <- outputDF_keyYM %>%
          mutate("IndEff_p_cor" = p_cor, .after = "IndEff_p")
        
        outputDF_keyYM <- outputDF_keyYM[order(outputDF_keyYM$IndEff_p_cor), ] # order by corrected p-value
                
        outputDF_keyYM_sig <- outputDF_keyYM %>% filter(IndEff_p_cor <= alpha)
        cat("\n File: ", file, "\n Variable: ", cog, ". Significant mediators:")
        if (length(outputDF_keyYM_sig$M) == 0) {
          smallest_NonSigP <- round(min(outputDF_keyYM$IndEff_p_cor), 3)
          cat("\n\t None. All corrected p >= ", smallest_NonSigP, "\n")
        } else {
          cat("(num. = ", length(outputDF_keyYM_sig$M), ") ")
          for (sigMediator in outputDF_keyYM_sig$M) {
            cat("\n\t", sigMediator)
          }
          outputDF_keyYM_nonSig <- outputDF_keyYM %>% filter(IndEff_p_cor > alpha)
          smallest_NonSigP <- round(min(outputDF_keyYM_nonSig$IndEff_p_cor), 3)
          cat("\n Remaining indirect effects, corrected p >= ", smallest_NonSigP)
        }
        cat("\n---------------------------")

        if (sigOnly == TRUE) {
          outputDF_keyYM <- outputDF_keyYM_sig
        }

        outputDf_p_cor <- rbind(outputDf_p_cor, outputDF_keyYM)
      }
    }

    saveName <- glue("{file}_{fileNameAffix}_{outputAffix}.csv")
    write.csv(outputDf_p_cor, file = saveName)
    outputDf_p_cor <- c()
  }
} # path: path to raw mediation output summary files; fileNames: name of summary mediation outcome files; pCor: logical specifying if p-vqlues qre to be corrected. If TRUE, then will correct considering all mediators listed in 'keyMediators' ; sigOnly: logicql specifying if output file should return only significant mediation models (uses alpha given by 'alpha'); outputAffix: string to append to end of output file name.

### UPDATE BELOW AS NEEDED ###
path <- "./AnalysisSummary" # Path for file containing anaylisis summary csv files.
##############################

fileNames <- c("UKBB_All_med_C_results_05_11_2022", "UKBB_noDx_med_C_results_05_11_2022", "UKBB_noDxNoSSRI_med_C_results_05_11_2022", "UKBB_noMed_med_C_results_05_11_2022", "UKBB_NoMedNoDx_med_C_results_05_11_2022", "UKBB_noSSRI_med_C_results_05_11_2022", "UKBB_oldestTert_med_C_results_05_11_2022", "UKBB_onlyDx_med_C_results_05_11_2022", "UKBB_onlyMed_med_C_results_05_11_2022", "UKBB_onlySSRI_med_C_results_05_11_2022", "UKBB_youngestTert_med_C_results_05_11_2022") # name of summary mediation outcome files

keyCognition <- c("cog_fluidIntel_score_t2_z", "cog_numMem_maxDigitRemem_t2_z") # name of cognition variables of interest for mediation analyses

# List of brain variables
## all Freesurfer DKT volumes, area, mThickness variables in UKBB
generalMetrics <- c("vol_BrainSeg_WB_t2", "vol_BrainSegNotVent_WB_t2", "vol_BrainSegNotVentSurf_WB_t2", "vol_SubCortGray_WB_t2", "vol_TotalGray_WB_t2", "vol_SupraTentorial_WB_t2", "vol_SupraTentorialNotVent_WB_t2", "vol_EstimatedTotalIntraCranial_WB_t2", "vol_BrainStem_WB_t2", "vol_CSF_WB_t2", "vol_wmHyperintensities_WB_t2", "vol_WMhypointensities_WB_t2", "vol_Hippocampus_L_t2", "vol_Hippocampus_R_t2", "vol_CerebellumCortex_R_t2", "vol_CerebellumCortex_L_t2", "vol_CerebellumWhiteMatter_L_t2", "vol_CerebellumWhiteMatter_R_t2")

aseg_brainVars <- c("vol_VentricleChoroid_WB_t2", "vol_3rdVentricle_WB_t2", "vol_4thVentricle_WB_t2", "vol_5thVentricle_WB_t2", "vol_nonWMhypointensities_WB_t2", "vol_OpticChiasm_WB_t2", "vol_CCPosterior_WB_t2", "vol_CCMidPosterior_WB_t2", "vol_CCCentral_WB_t2", "vol_CCMidAnterior_WB_t2", "vol_CCAnterior_WB_t2", "area_Accumbensarea_L_t2", "vol_Cortex_L_t2", "vol_CerebralWhiteMatter_L_t2", "vol_LateralVentricle_L_t2", "vol_InfLatVent_L_t2", "vol_thalamusProper_L_t2", "vol_Caudate_L_t2", "vol_Putamen_L_t2", "vol_Pallidum_L_t2", "vol_Amygdala_L_t2", "vol_Accumbensarea_L_t2", "vol_VentralDC_L_t2", "vol_choroidplexus_L_t2", "area_Accumbensarea_R_t2", "vol_Cortex_R_t2", "vol_CerebralWhiteMatter_R_t2", "vol_LateralVentricle_R_t2", "vol_InfLatVent_R_t2", "vol_thalamusProper_R_t2", "vol_Caudate_R_t2", "vol_Putamen_R_t2", "vol_Pallidum_R_t2", "vol_Amygdala_R_t2", "vol_Accumbensarea_R_t2", "vol_VentralDC_R_t2", "vol_vessel_L_t2", "vol_vessel_R_t2", "vol_choroidplexus_R_t2") # all Freesurfer ASEG variables in UKBB

dkt_brainVars_area <- c("area_Caudalanteriorcingulate_L_t2", "area_Caudalmiddlefrontal_L_t2", "area_Cuneus_L_t2", "area_Entorhinal_L_t2", "area_Fusiform_L_t2", "area_Inferiorparietal_L_t2", "area_Inferiortemporal_L_t2", "area_Isthmuscingulate_L_t2", "area_Lateraloccipital_L_t2", "area_Lateralorbitofrontal_L_t2", "area_Lingual_L_t2", "area_Medialorbitofrontal_L_t2", "area_Middletemporal_L_t2", "area_Parahippocampal_L_t2", "area_Paracentral_L_t2", "area_Parsopercularis_L_t2", "area_Parsorbitalis_L_t2", "area_Parstriangularis_L_t2", "area_Pericalcarine_L_t2", "area_Postcentral_L_t2", "area_Posteriorcingulate_L_t2", "area_Precentral_L_t2", "area_Precuneus_L_t2", "area_Rostralanteriorcingulate_L_t2", "area_Rostralmiddlefrontal_L_t2", "area_Superiorfrontal_L_t2", "area_Superiorparietal_L_t2", "area_Superiortemporal_L_t2", "area_Supramarginal_L_t2", "area_Transversetemporal_L_t2", "area_Insula_L_t2", "area_Caudalanteriorcingulate_R_t2", "area_Caudalmiddlefrontal_R_t2", "area_Cuneus_R_t2", "area_Entorhinal_R_t2", "area_Fusiform_R_t2", "area_Inferiorparietal_R_t2", "area_Inferiortemporal_R_t2", "area_Isthmuscingulate_R_t2", "area_Lateraloccipital_R_t2", "area_Lateralorbitofrontal_R_t2", "area_Lingual_R_t2", "area_Medialorbitofrontal_R_t2", "area_Middletemporal_R_t2", "area_Parahippocampal_R_t2", "area_Paracentral_R_t2", "area_Parsopercularis_R_t2", "area_Parsorbitalis_R_t2", "area_Parstriangularis_R_t2", "area_Pericalcarine_R_t2", "area_Postcentral_R_t2", "area_Posteriorcingulate_R_t2", "area_Precentral_R_t2", "area_Precuneus_R_t2", "area_Rostralanteriorcingulate_R_t2", "area_Rostralmiddlefrontal_R_t2", "area_Superiorfrontal_R_t2", "area_Superiorparietal_R_t2", "area_Superiortemporal_R_t2", "area_Supramarginal_R_t2", "area_Transversetemporal_R_t2", "area_Insula_R_t2")

dkt_brainVars_vol <- c("vol_Caudalanteriorcingulate_R_t2", "vol_Caudalmiddlefrontal_R_t2", "vol_Cuneus_R_t2", "vol_Entorhinal_R_t2", "vol_Fusiform_R_t2", "vol_Inferiorparietal_R_t2", "vol_Inferiortemporal_R_t2", "vol_Isthmuscingulate_R_t2", "vol_Lateraloccipital_R_t2", "vol_Lateralorbitofrontal_R_t2", "vol_Lingual_R_t2", "vol_Medialorbitofrontal_R_t2", "vol_Middletemporal_R_t2", "vol_Parahippocampal_R_t2", "vol_Paracentral_R_t2", "vol_Parsopercularis_R_t2", "vol_Parsorbitalis_R_t2", "vol_Parstriangularis_R_t2", "vol_Pericalcarine_R_t2", "vol_Postcentral_R_t2", "vol_Posteriorcingulate_R_t2", "vol_Precentral_R_t2", "vol_Precuneus_R_t2", "vol_Rostralanteriorcingulate_R_t2", "vol_Rostralmiddlefrontal_R_t2", "vol_Superiorfrontal_R_t2", "vol_Superiorparietal_R_t2", "vol_Superiortemporal_R_t2", "vol_Supramarginal_R_t2", "vol_Transversetemporal_R_t2", "vol_Insula_R_t2", "vol_Caudalanteriorcingulate_L_t2", "vol_Caudalmiddlefrontal_L_t2", "vol_Cuneus_L_t2", "vol_Entorhinal_L_t2", "vol_Fusiform_L_t2", "vol_Inferiorparietal_L_t2", "vol_Inferiortemporal_L_t2", "vol_Isthmuscingulate_L_t2", "vol_Lateraloccipital_L_t2", "vol_Lateralorbitofrontal_L_t2", "vol_Lingual_L_t2", "vol_Medialorbitofrontal_L_t2", "vol_Middletemporal_L_t2", "vol_Parahippocampal_L_t2", "vol_Paracentral_L_t2", "vol_Parsopercularis_L_t2", "vol_Parsorbitalis_L_t2", "vol_Parstriangularis_L_t2", "vol_Pericalcarine_L_t2", "vol_Postcentral_L_t2", "vol_Posteriorcingulate_L_t2", "vol_Precentral_L_t2", "vol_Precuneus_L_t2", "vol_Rostralanteriorcingulate_L_t2", "vol_Rostralmiddlefrontal_L_t2", "vol_Superiorfrontal_L_t2", "vol_Superiorparietal_L_t2", "vol_Superiortemporal_L_t2", "vol_Supramarginal_L_t2", "vol_Transversetemporal_L_t2", "vol_Insula_L_t2")

dkt_brainVars_mThick <- c("mThick_Caudalanteriorcingulate_L_t2", "mThick_Caudalmiddlefrontal_L_t2", "mThick_Cuneus_L_t2", "mThick_Entorhinal_L_t2", "mThick_Fusiform_L_t2", "mThick_Inferiorparietal_L_t2", "mThick_Inferiortemporal_L_t2", "mThick_Isthmuscingulate_L_t2", "mThick_Lateraloccipital_L_t2", "mThick_Lateralorbitofrontal_L_t2", "mThick_Lingual_L_t2", "mThick_Medialorbitofrontal_L_t2", "mThick_Middletemporal_L_t2", "mThick_Parahippocampal_L_t2", "mThick_Paracentral_L_t2", "mThick_Parsopercularis_L_t2", "mThick_Parsorbitalis_L_t2", "mThick_Parstriangularis_L_t2", "mThick_Pericalcarine_L_t2", "mThick_Postcentral_L_t2", "mThick_Posteriorcingulate_L_t2", "mThick_Precentral_L_t2", "mThick_Precuneus_L_t2", "mThick_Rostralanteriorcingulate_L_t2", "mThick_Rostralmiddlefrontal_L_t2", "mThick_Superiorfrontal_L_t2", "mThick_Superiorparietal_L_t2", "mThick_Superiortemporal_L_t2", "mThick_Supramarginal_L_t2", "mThick_Transversetemporal_L_t2", "mThick_Insula_L_t2", "mThick_Caudalanteriorcingulate_R_t2", "mThick_Caudalmiddlefrontal_R_t2", "mThick_Cuneus_R_t2", "mThick_Entorhinal_R_t2", "mThick_Fusiform_R_t2", "mThick_Inferiorparietal_R_t2", "mThick_Inferiortemporal_R_t2", "mThick_Isthmuscingulate_R_t2", "mThick_Lateraloccipital_R_t2", "mThick_Lateralorbitofrontal_R_t2", "mThick_Lingual_R_t2", "mThick_Medialorbitofrontal_R_t2", "mThick_Middletemporal_R_t2", "mThick_Parahippocampal_R_t2", "mThick_Paracentral_R_t2", "mThick_Parsopercularis_R_t2", "mThick_Parsorbitalis_R_t2", "mThick_Parstriangularis_R_t2", "mThick_Pericalcarine_R_t2", "mThick_Postcentral_R_t2", "mThick_Posteriorcingulate_R_t2", "mThick_Precentral_R_t2", "mThick_Precuneus_R_t2", "mThick_Rostralanteriorcingulate_R_t2", "mThick_Rostralmiddlefrontal_R_t2", "mThick_Superiorfrontal_R_t2", "mThick_Superiorparietal_R_t2", "mThick_Superiortemporal_R_t2", "mThick_Supramarginal_R_t2", "mThick_Transversetemporal_R_t2", "mThick_Insula_R_t2")

## Names of Brain variables computed
lobularBrainMetrics <- c("area_insula_WB_t2", "vol_insula_WB_t2", "mThick_insula_WB_t2")
for (i in c("area", "mThick", "vol")) {
  for (j in c("frontal", "parietal", "occipital", "temporal")) {
    for (k in c("L", "R", "WB")) {
      lobularBrainMetrics <- c(lobularBrainMetrics, glue("{i}_{j}_{k}_t2"))
      # print(glue("{i}_{j}_{k}"))
    }
  }
}

keyMediators_1 <- c(lobularBrainMetrics, generalMetrics)
keyMediators_1_z <- unlist(lapply(keyMediators_1, function(.)paste(., "_z", sep = "")))
getPMed(path = path, fileNames = fileNames, keyMediators = keyMediators_1_z, keyCognition = keyCognition, pCor = T, sigOnly = F, alpha = .1, outputAffix = "round1") # for round 1, general brain metrics

keyMediators_2 <- c(aseg_brainVars, dkt_brainVars_area, dkt_brainVars_vol, dkt_brainVars_mThick) # list of mediators to correct p-values for round 2
keyMediators_2_z <- unlist(lapply(keyMediators_2, function(.)paste(., "_z", sep = "")))
getPMed(path = path, fileNames, keyMediators = keyMediators_2_z, keyCognition = keyCognition, pCor = T, sigOnly = F, alpha = .1, outputAffix = "round2") # for round 2, specific brain regions
