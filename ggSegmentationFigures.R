# Creates figure with colour scale associated to each area on the DKT atlas

library(glue)

formatData <- function(file){
  require(readr)
  require(dplyr)
  require(glue)
  
  df <- tibble(read_csv(file, show_col_types = FALSE)) # Import data
  df_split <- split(df, df$Y) # split dataset by unique cognitive outcomes
  df_formatted <- list()
  
  for(i in 1:length(df_split)){
    
    df_i <- as.data.frame(df_split[i])
    dfName <- names(df_split)[i]
    colOfInt <- glue("{dfName}.M") # By splitting df, columns get rename to {dfName}.[colName]
    df_i <- df_i %>%
      tidyr::separate(col = colOfInt, sep = "_", into = c("Metric", "Region", "Hand", glue("{dfName}.discard1"), glue("{dfName}.discard2"))) %>% 
      mutate(Metric = case_when(
        Metric == "vol" ~ "Volume",
        Metric == "mThick" ~ "Thickness",
        Metric == "area" ~ "Area",
        TRUE ~ "other"
      )) %>% 
      mutate(Region = case_when(
        Hand == "L" ~ tolower(as.character(glue("lh_{Region}"))),
        Hand == "R" ~ tolower(as.character(glue("rh_{Region}"))),
        Hand == "WB" ~ tolower(as.character(glue("wb_{Region}"))),
        TRUE ~ Region
      )) %>%
      select(-glue("{dfName}.discard1"), -glue("{dfName}.discard2"), -Hand)
    df_formatted[[i]] <- df_i
  } # Iterate through unique cognitive outcomes, formatting appropriately for each
  
  names(df_formatted) <- names(df_split)
  return(df_formatted)
} # Format mediation summary file to appropriate format for 'createBrainFigure' function. I.e., one row per ROI, metric column, proper names for regions. Input: dataframe with results for all cognitive outcomes of interest. Returns list of dataframes, one per cognitive outcome of interest.

createBrainFigure <- function(file, figureValue, savePath, Y, subset){
  require(stringr)
  require(ggplot2)
  require(ggseg)
  require(ggthemes)
  require(tidyr)
  
  date <- str_c(format(Sys.time(), "%m"), "_", format(Sys.time(), "%d"), "_", format(Sys.time(), "%Y")) # set todays date for easier output filenaming
  
  # Extract name of cognitive performance variable
  Yname <- unlist(strsplit(Y, split = "_"))
  Yname <- Yname[c(-which(Yname == "cog"),-(which(Yname == "t2"):length(Yname)))]
  Yname_title <- ""
  Yname_file <- ""
  for(i in 1:length(Yname)){
    if(i == 1){
      Yname_file <- Yname[i]
      Yname_title <- Yname[i]
    } else {
      Yname_file <- glue("{Yname_file}_{Yname[i]}")
      Yname_title <- glue("{Yname_title} {Yname[i]}")
    }
  }

  if(Yname_title == "numMem maxDigitRemem"){
    Yname_title <- "Numeric Memory - Max. digits remembered"
  } else if(Yname_title == "fluidIntel score") {
    Yname_title <- "Fluid Intelligence - Score"
  }
  
  if(figureValue == "p"){
    
    valueName <- "p (uncorrected) of indirect effect"
    value <- glue("{Y}.IndEff_p")
    figureName <- glue("{subset}_{Yname_file}_p")
    title <- "p-value"

    breakValues <- c(.01, .05, .1, .2, .5)
    gradientSettings <- scale_fill_gradientn(
      name = valueName,
      colours = c("#062635", "#0A4C6A", "#12719E", "#1696D2", "#73BFE2","#A2D4EC","#CFE8F3"), # Lower is better; full colour scale: c("#062635", "#0A4C6A", "#12719E", "#1696D2", "#46ABDB","#73BFE2","#A2D4EC","#CFE8F3"),
      limits = c(0,1),
      guide = "colourbar",
      values = c(0, breakValues, 1),
      breaks = breakValues,
      labels = c(".01", ".05", ".1", ".2", ".5") # Make scale for break values defined in if statements above
    )
    
  } else if(figureValue == "pCor") {
    
    valueName <- "corrected p of indirect effect"
    value <- glue("{Y}.IndEff_p_cor")
    figureName <- glue("{subset}_{Yname_file}_pCor")
    title <- "p-value (corrected)"

    breakValues <- c(.01, .05, .1, .2, .5)
    gradientSettings <- scale_fill_gradientn(
      name = valueName,
      colours = c("#062635", "#0A4C6A", "#12719E", "#1696D2", "#73BFE2","#A2D4EC","#CFE8F3"), # Lower is better; full colour scale: c("#062635", "#0A4C6A", "#12719E", "#1696D2", "#46ABDB","#73BFE2","#A2D4EC","#CFE8F3"),
      limits = c(0,1),
      guide = "colourbar",
      values = c(0, breakValues, 1),
      breaks = breakValues,
      labels = c(".01", ".05", ".1", ".2", ".5") # Make scale for break values defined in if statements above
    )
    
  } else if(figureValue == "indEff"){
    
    valueName <- "Indirect effect parameter estimate"
    value <- glue("{Y}.IndEff_est_correct")
    figureName <- glue("{subset}_{Yname_file}_indEff_est_cor")
    title <- "Indirect effect"
    
    lowBound <- 0
    hiBound <- .2 # !!! TBD !!!
    gradientCodes <- c("#CFE8F3", "#A2D4EC", "#73BFE2","#46ABDB", "#1696D2", "#12719E", "#0A4C6A", "#062635") # Higher is better
    breakValues <- c(.0001, .001, .01, .05, .1) # !!! PROPER BOUNDS TO BE DETERMINED
    breakValueLabels <- c("", ".001", ".01", ".05", ".1")
    
    gradientSettings <- scale_fill_gradientn(
      name = valueName,
      colours = c("#062635", "#0A4C6A", "#12719E", "#1696D2", "#73BFE2","#A2D4EC","#CFE8F3"), # Lower is better; full colour scale: c("#062635", "#0A4C6A", "#12719E", "#1696D2", "#46ABDB","#73BFE2","#A2D4EC","#CFE8F3"),
      limits = c(lowBound,hiBound),
      guide = "colourbar",
      values = c(lowBound, breakValues, hiBound ),
      breaks = breakValues,
      labels = breakValueLabels # Make scale for break values defined in if statements above
    )
    
  } else {
    
    return(cat("Error: ", figureValue, " is invalid. \n Options are 'p', 'pCor', or 'indEff'"))
    
  }
  
  cat("\n --------- \n Value:", value, "\n")
    
  regionColNum <- which(colnames(file) == glue("{Y}.Region"))

  #1 Extract DKT brain regions ---------
  ##A makes a list of all ROIs to keep for DKT figures.
  DKTRegions <- brain_labels(dk) # name of all DKT labels
  DKT_matches <- c()
  for(roi in DKTRegions){
    new_matches <- file[[regionColNum]][grepl(roi, file[[regionColNum]], ignore.case = TRUE, fixed = FALSE) == TRUE]
    DKT_matches <- c(DKT_matches, unique(new_matches))
  } 
  DKT_matches <- unique(DKT_matches)
  # str(matches)
  
  ##B Extract only rows with value of column 'Region' specified in matches. Extract all columns.
  reallyShort_DK <- subset(file, file[[regionColNum]] %in% DKT_matches) 
  # str(reallyShort_DK)
  
  ##C Make and format df with only necessary variables. 
  valueColNum <- which(colnames(reallyShort_DK) == glue("{Y}.{value}"))
  regionColNum <- which(colnames(reallyShort_DK) == glue("{Y}.Region"))
  metricColNum <- which(colnames(reallyShort_DK) == glue("{Y}.Metric"))

  # cat("- Debugging - \n \t Structure: \n ", str(file[valueColNum]), "\n \t Table: \n ", table(file[valueColNum]))
  
  reallyShort_DK <- reallyShort_DK[,c(regionColNum, metricColNum, valueColNum)] 
  colnames(reallyShort_DK) <- c("label", "measure", "dfValue")
  reallyShort_DK$label <- as.factor(reallyShort_DK$label)
  # View(reallyShort_DK)
  # str(reallyShort_DK)
  
  #3 Make figures
  cortical_pos <- c("left lateral", "left medial", "right medial", "right lateral")
  ##A DKT figure
  DKTFigure <- reallyShort_DK %>%
    group_by(measure) %>%
    ggplot() +
    ggtitle(Yname_title) +
    theme_tufte() +
    geom_brain(atlas = dk,
               aes(fill = dfValue),
               position =  position_brain(cortical_pos),
               show.legend = TRUE,
    ) +
    gradientSettings + 
    scale_x_continuous(
      breaks = c(140, 520, 900, 1280), 
      labels = str_to_title(cortical_pos)) +
    theme(
      plot.background = element_rect(fill = "white", size = 0),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(vjust = .75),
      legend.text = element_text(angle = 300, vjust = .8, hjust = .3),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.y = element_text(angle = 0, hjust = 0)) +
    guides(fill = guide_colorbar(barwidth = 15, raster = TRUE, direction = "horizontal", label.position = "bottom"))+ # label.vjust = 0, label.hjust = 0
    facet_grid(measure ~ .)
  
  (cat("Saving DKT figure: ", glue("{figureName}_DK_{date}.jpg"), "\n") + ggsave(glue("{savePath}/{figureName}_DK_{date}.jpg"), dpi=300))
} # Creates ggSegmentation for DSK atlas. Takes: 'file' a data file formatted appropriately (see ROIS.csv for reference); 'figureValue' string specifying the values the figure is to represent (either 'p', 'indEff', or 'pCor'); 'savePath' path where figure should be saved; 'Y' dependent variable name
## N.b. BANKSSTS not showing data

# input files of interest
PCorFilePath <- "OutputFiles/round2Output"
prefix <- "UKBB_"
subsets <- c("All", "noMed", "NoMedNoDx", "noDxNoSSRI", "noDx", "noSSRI", "oldestTert",  "onlyDx", "onlyMed","onlySSRI", "youngestTert")
infix <- "_med_C_results_05_11_2022_corrected_06_27_2022_pCor_round2" # round 2 has the specific brain regions
subset <- c("noDxNoSSRI", "noMed") # for testing purposes
figureValueToPlot <- "indEff" # either 'p', 'indEff', or 'pCor'

for(subset in subset){
  cat("\n---- ", subset, "-----------------------\n")
  # load original file for that subset
  summaryFileName <- glue("./{PCorFilePath}/{prefix}{subset}{infix}.csv") 
  formattedFile <- formatData(file = summaryFileName)
  # str(formattedFile)
  
  for(i in 1:length(formattedFile)){
    cogVar <- names(formattedFile[i])
    createBrainFigure(file = as.data.frame(formattedFile[i]), figureValue = figureValueToPlot, savePath = "./OutputFiles/Figures", Y = cogVar, subset = subset)
  }
} # generate figures for each subset and for each cognitive variable

