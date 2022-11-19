# Create mediation diagrams from mediation analyses previously conducted.
## * sig at .05, ** at .01, *** at .001

med_diagram <- function(mediation_data, save, outputName){
  # adapted from Omar Wasow's post on https://stackoverflow.com/questions/46465752/drawing-simple-mediation-diagram-in-r
  require(glue)
  require(DiagrammeR)
  require(DiagrammeRsvg)
  require(rsvg)
  
  height = .5
  width = 2.5
  graph_label = mediation_data$main
  node_text_size = 14
  edge_text_size = 14
  color = "black"
  ranksep = .35
  minlen = 3
  engine = "dot"
  
  mediation_data$height  <- height   # node height
  mediation_data$width   <- width    # node width
  mediation_data$color   <- color    # node + edge border color
  mediation_data$ranksep <- ranksep  # separation btwn mediator row and x->y row
  mediation_data$minlen  <- minlen   # minimum edge length
  mediation_data$engine  <- engine
  
  mediation_data$node_text_size  <- node_text_size
  mediation_data$edge_text_size  <- edge_text_size
  
  mediation_data$graph_label <- ifelse(is.na(graph_label), "", paste0("label = '", graph_label, "'"))
  
  diagram_out <- glue::glue_data(mediation_data,
                                 "digraph flowchart {
        fontname = Times
        <<graph_label>>
        graph [ranksep = <<ranksep>>]
  
        # node definitions with substituted label text
        node [fontname = Times, shape = rectangle, fixedsize = TRUE, width = <<width>>, height = <<height>>, fontsize = <<node_text_size>>, color = <<color>>]        
          mm [label = '<<lab_m>>']
          xx [label = '<<lab_x>>']
          yy [label = '<<lab_y>>']
  
        # edge definitions with the node IDs
        edge [minlen = <<minlen>>, fontname = Times, fontsize = <<edge_text_size>>, color = <<color>>]
          mm -> yy [label = '<<coef_my>>'];
          xx -> mm [label = '<<coef_xm>>'];
          xx -> yy [label = '<<coef_xy>>'];
        
        { rank = same; mm }
        { rank = same; xx; yy }
        
        }
        ", .open = "<<", .close = ">>")  
  
  if(save == T){
    temp <- export_svg(DiagrammeR::grViz(diagram_out)) 
    temp <- charToRaw(temp)
    temp <- rsvg::rsvg_png(temp, outputName)
  } else {
    DiagrammeR::grViz(diagram_out)
  }
} # Creates mediation diagram. 'mediation_data' matrix with strings to be added to med diagram; 'save' logical indicating if the diagram should be saved; 'outputName' name for file to be saved.

round1_MedDiagrams <- function(pCorFilePath, pCorFileName, medFilePath, medFileName, rois, hemispheres, metrics, IndEffpCor, outputPath){
  require(glue)
  require(stringr)
  require(readr)

  # Ensure both df files are for same subset
  pCor_subset <- unlist(strsplit(pCorFileName, split = "_"))[2] # extract subset name
  medDet_subset <-  unlist(strsplit(medFileName, split = "_"))[2]
  if(pCor_subset != medDet_subset){
    cat("Error. Subset of the files provided are not the same. Check file names provided as \`pCorFileName\` and \`medFileName\` . \n Aborting mediation figure creation.")
  } else {
    date <- str_c(format(Sys.time(), "%m"), "_", format(Sys.time(), "%d"), "_", format(Sys.time(), "%Y"))
    
    pToStar <- function(p){
      if(p <= .001){
        p <- "***"
      } else if(p <= .01){
        p <- "**"
      } else if(p <= .05){
        p <- "*"
      } else {
        p <- ""
      }
      return(p)
    } # Format p-values as '*', '**', '***' or ''.
    
    subset <- pCor_subset
    cat("Making mediation figures for: ", subset, "\n")
    df_pCor <- read_csv(glue("{pCorFilePath}/{pCorFileName}"), show_col_types = F)
    df_medDetails <- read_csv(glue("{medFilePath}/{medFileName}"), show_col_types = F)
    
    # find colNum of relevant columns in both dfs
    medDet_colNum_aPath_b <- which(colnames(df_medDetails) == "aPath_b")
    medDet_colNum_aPath_p <- which(colnames(df_medDetails) == "aPath_model_p")
    medDet_colNum_bPath_b <- which(colnames(df_medDetails) == "bPath_b")
    medDet_colNum_bPath_p <- which(colnames(df_medDetails) == "bPath_model_p")
    medDet_colNum_dirEff_b <- which(colnames(df_medDetails) == "dirEff_b")
    medDet_colNum_dirEff_p <- which(colnames(df_medDetails) == "dirEff_p")
    medDet_colNum_indEff_b <- which(colnames(df_medDetails) == "IndEff_est_correct")
    
    pCor_colNum_indEff_pCor <- which(colnames(df_pCor) == "IndEff_p_cor")
    pCor_colNum_indEff_p <- which(colnames(df_pCor) == "IndEff_p")
    pCor_colNum_Y <- which(colnames(df_pCor) == "Y")
    pCor_colNum_M <- which(colnames(df_pCor) == "M")
    
    # For each roi, hemisphere, metric combination ----
    for(roi in rois){
      for(hemi in hemispheres){
        for(metric in metrics){
          
          medName <- glue("{metric}_{roi}_{hemi}") # for output file name
          cat("Analysis name: ", medName, "\n")
          
          dfMed <- glue("{medName}_t2_z") # for finding the unique mediator in df
          
          if(dfMed %in% df_medDetails$M && dfMed %in% df_pCor$M) { # Check that this combination exists. If it doesn't skip
            
            pCorRows <- which(df_pCor$M == dfMed) # determine the rows with analyses from this unique mediator in the corrected p-value file
            
            for(pCorRow in pCorRows) { # iterate through rows containing this mediator
              ## Determine labels for mediation diagram. ----
              ### Determine equivalent rows in both dfs. ----
              Yval <- unlist(df_pCor[pCorRow, pCor_colNum_Y]) # Extract value in column Y
              medDet_row <- which(df_medDetails$M == dfMed)[which(df_medDetails$M == dfMed) %in% which(df_medDetails$Y == Yval)] # find row in df_medDetails that corresponds to pCorRow
              
              ### Extract values for diagram ----
              # format of parameter estimates: 0.###
              
              a_eff <- round(unlist(df_medDetails[medDet_row, medDet_colNum_aPath_b]), 3) # 'aPath_b'
              a_p <- pToStar(df_medDetails[medDet_row, medDet_colNum_aPath_p]) # 'aPath_model_p'
              
              b_eff <- round(unlist(df_medDetails[medDet_row, medDet_colNum_bPath_b]), 3) # 'bPath_b'
              b_p <- pToStar(df_medDetails[medDet_row,medDet_colNum_bPath_p]) # 'bPath_model_p'
              
              c_eff <- round(unlist( df_medDetails[medDet_row, medDet_colNum_dirEff_b]), 3) # 'dirEff_b'
              c_p <- pToStar(df_medDetails[medDet_row, medDet_colNum_dirEff_p]) # 'dirEff_p'
              
              ab_indEff <- format(round(unlist(df_medDetails[medDet_row, medDet_colNum_indEff_b]), 3), nsmall = 3) # colName: "IndEff_est_correct"
              
              if(IndEffpCor == T){
                ab_p <- pToStar(df_pCor[pCorRow, pCor_colNum_indEff_pCor]) # colname: "IndEff_p_cor"
              } else if (IndEffpCor == F){
                ab_p <- pToStar(df_pCor[pCorRow, pCor_colNum_indEff_p]) # colname: "IndEff_p"
              } else {
                return(car("Invalid value for 'IndEffpCor'. Possible values are: \`T\` (corrected p is desired) or \`F\` (non-corrected p is desired)"))
              }
              
              #### Make labels. ----
              ##### Line labels -----
              a <- glue("a = {a_eff}{a_p}") # Text for A path in form: "{ind eff estimate (0.###)} {significance ('*', '**', '***' or '')}"
              b <- glue("b = {b_eff}{b_p}") # Text for B path in form: "{ind eff estimate (0.###)} {significance ('*', '**', '***' or '')}"
              c <- glue("c = {c_eff}{c_p} \n ab = {ab_indEff}{ab_p}") # Text for C and indirect effect (ab) in form: "c = {ind eff estimate (0.###)} {significance ('*', '**', '***' or '')} \n ab = {ind eff estimate (0.###)} {significance ('*', '**', '***' or '')}"
              
              ##### Mediator box labels -----
              if(metric == "vol"){
                metricLabel = "vol."
              } else if(metric == "mThick"){
                metricLabel = "thick."
              } else {
                metricLabel = metric
              }
              
              if(roi == "frontal"){
                medLabel = glue("Frontal lobe {metricLabel} ({hemi})")
              } else if(roi == "parietal"){
                medLabel = glue("Parietal lobe {metricLabel} ({hemi})")
              } else if(roi == "occipital"){
                medLabel = glue("Occipital lobe {metricLabel} ({hemi})")
              } else if(roi == "temporal"){
                medLabel = glue("Temporal lobe {metricLabel} ({hemi})")
              } else if(roi == "hippocampus"){
                medLabel = glue("Hippocampal {metricLabel} ({hemi})")
              } else {
                medLabel = medName
              }
              
              ### Cognitive-variable box label. ----
              if(Yval == "cog_fluidIntel_score_t2_z"){
                lab_y <- "Fluid Intelligence Score"
                cog <- "FluidIntelScore"
              } else if(Yval == "cog_numMem_maxDigitRemem_t2_z"){
                lab_y <- "Num. Mem. (Max dig.)"
                cog <- "NumMemMaxDig"
              } else {
                lab_y <- Yval
                cog <- Yval
              }
              
              ## Make and save diagram. ----
              outputName <- glue("medDia_{subset}_{cog}_{medName}_{date}.png") # name of diagram when saved
              med_diagram(mediation_data =  data.frame(
                main = "",  
                lab_x   = "C-Reactive Protein", # it is assumed that the independent variable is CRP
                lab_m   = medLabel, # name of mediator
                lab_y   = lab_y, # name of cognitive test
                coef_xm = a, 
                coef_my = b, 
                coef_xy = c), 
                save = T, outputName = outputName)  
              cat("\t File \`", outputName, "\` saved. \n")
            }
          } else {
            cat("Error. \`", medName, "\` is an invalid combination; it is not a mediator in at least one of the provided files. \n This combination will be skipped without making a mediation diagram. \n")
          }
        }
      }
    }
  }
} # Automatically make and save diagrams for specifies brain regions and metrics. 'pCorFilePath' string of path to folder containing output of mediation analyses !!with corrected p-values!!; 'pCorFileName' string of name of df for specific data subset containing mediation analyses with corrected p-value for indirect effect. Df must contain the column names: 'M', 'Y', 'IndEff_est', 'IndEff_p_cor'; 'medFilePath' string of path to folder containing output of mediation analyses. 'medFileName' string of name of df for specific data subset. Df must contain the column names: 'M', 'Y', 'aPath_model_p', 'aPath_model_p' 'bPath_model_p', 'bPath_model_p'; 'rois' list of brain regions to make mediation models for. Values can be any name found; 'hemispheres' hemispheres of interest. Can be 'L', 'R' or 'WB'; 'metrics' list of metrics of interest. Can be 'area', 'mThick', and/or 'vol'; 'IndEffpCor' boolean specifying if corrected p-value of ind effect is to be used; 'outputPath' string specifying where to save output diagrams

round1_MedDiagrams(pCorFilePath = "./OutputFiles/round1Output", 
                pCorFileName = "UKBB_noDxNoSSRI_med_C_results_05_11_2022_corrected_06_27_2022_pCor_round1.csv", 
                medFilePath = "./OutputFiles/CorrectedSummaryFiles", 
                medFileName = "UKBB_noDxNoSSRI_med_C_results_05_11_2022_corrected_06_27_2022.csv", 
                rois = c("frontal", "parietal", "occipital", "temporal", "hippocampus"),
                hemispheres = c("WB", "L", "R"),
                metrics = c("area", "mThick", "vol"),
                IndEffpCor = T, 
                outputPath = "./OutputFiles/Figures/MedModels")

# pCorFilePath = "./OutputFiles/round1Output"
# pCorFileName = "UKBB_noDxNoSSRI_med_C_results_05_11_2022_corrected_06_27_2022_pCor_round1.csv" 
# medFilePath = "./OutputFiles/CorrectedSummaryFiles"
# medFileName = "UKBB_noDxNoSSRI_med_C_results_05_11_2022_corrected_06_27_2022.csv"
# rois = c("frontal") # c("frontal", "parietal", "occipital", "temporal", "hippocampus")
# hemispheres = c("WB")  # c("WB", "L", "R")
# metrics = c("area")  # c("area", "mThick", "vol")
# IndEffpCor = T
# outputPath = "./OutputFiles/Figures/MedModels"

# Manual diagrams ----
# medDigram_a <- med_diagram(mediation_data = data.frame(
#   main = "",  
#   lab_x   = "C-Reactive Protein",
#   lab_m   = "L Cerebellum cortex volume",
#   lab_y   = "Fluid Intelligence Score",
#   coef_xm = "-0.026 ***",
#   coef_my = "-0.033 ***",
#   coef_xy = "c\` = -0.034 *** \n ab = -0.001"), save = T, outputName = "test.png")
