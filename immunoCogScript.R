# 0. Preamble ----
###
# This file contains is the analysis script for the immunology, cognition and brain mediation analysis of the UKBB.
# The following script was computed on a system running Ubuntu 22.04.1 LTS with an x86 platform.
# Prepared by Daniel Mendelson with help from Katie M. Lavigne and Joshua Unrau.
# Work performed under the supervision of Dr. Martin Lepage.
# 14 December, 2022
###

## Load packages ----

library(DescTools)
library(DiagrammeR)
library(DiagrammeRsvg)
library(dplyr)
library(effsize)
library(factoextra)
library(ggplot2)
library(ggseg)
library(ggthemes)
library(glue)
library(lubridate) 
library(psych)
library(rcompanion)
library(readr)
library(RMediation)
library(rsvg)
library(stats)
library(stringr)
library(tidyr)

# 1. Data preperation -----------
## 1.A Functions ----
getDate <- function(){
  require(stringr)
  return(str_c(format(Sys.time(), "%m"), "_", format(Sys.time(), "%d"), "_", format(Sys.time(), "%Y")))
}

filter_rows <- function(df, ...) {
  require(dplyr)
  vars = as.list(substitute(list(...)))[-1L]
  for(arg in vars) {
    df_new <- df %>% filter(!!arg)
    rows_filtered <- nrow(df) - nrow(df_new)
    cat(sprintf('Filtered out %s rows using: %s\n', rows_filtered, deparse(arg)))
    df = df_new
  }
  return(df_new)
} # Filters rows according to criteria and prints number of cases filtered per criteria. https://stackoverflow.com/questions/41021976/print-number-of-rows-filtered-out-by-dplyrs-filter-function
### Check data types and run non-parametric analyses ----
numNotFactor <- function(x){
  if(is.numeric(x) == TRUE && is.factor(x) == FALSE && is.logical(x) == FALSE){
    return(x)
  }
} # This function is made to be used by mapply(). 'x' - a column from a df.

normalCheck <- function(x, colName, colNum){
  require(psych)
  plotList <- c()
  
  degFree <- length(x)
  skew <- skew(x)
  se_skew <- skew/length(x)
  t_skew <- skew/se_skew
  p_skew <- pt(q = abs(t_skew), df = degFree, lower.tail = F)
  skewReport <- paste("skew=", round(skew, 2), "(t=", format(round(t_skew, 2), scientific = T), ", p=", format(round(p_skew, 2), scientific = T), ")")
  
  kurt <- kurtosi(x)
  se_kurt <- sqrt(24/degFree)
  t_kurt <- kurt/se_kurt
  p_kurt <- pt(q = abs(t_kurt), df = degFree, lower.tail = F)
  kurtReport <- paste("kurt=", round(kurt, 2), "(t=", format(round(t_kurt, 2), scientific = T), ", p=", format(round(p_kurt, 2), scientific = T), ")")
  
  layout(matrix(c(1,2), ncol = 2, nrow = 1))
  hist <- hist(x, xlab = "",main = paste("Hist - ", colNum, colName), sub = paste("colNum: ", colNum, skewReport, "\n", "colNum: ", colNum, kurtReport))
  qqplot <- qqnorm(x, main = paste("QQplot - ", colNum, colName))
  qqline(x)
  iPlots <- list(hist, qqplot)
  plotList <- append(plotList, iPlots)
  return(plotList)
} # This function is made to be used by mapply(). 'x' - a column from a df; 'colName' - a string specifying the name of this column; 'colNum' - a number specifying the column number of this variable in df

runWilcox <- function(varsOfInterest, df, by, outputDir, outputName){
  require(rcompanion)
  require(stringr)
  require(glue)
  
  date <- getDate() # set todays date for easier output filenaming
  # wilcoxOutput <- data.frame("value type" = c("statistic","p", "n1", "n2", "r(effSize)"), row.names = 1) # initialize dataframe for the output
  wilcoxOutput <- data.frame("value type" = c("statistic","p","n_subset", "n_exclude", "r(effSize)"), row.names = 1) # initialize dataframe for the output
  i <- integer(0)
  for(i in varsOfInterest){ # iterate through the variables of interest. Run wilcoxon-paired rank test. Put results into a dataframe
    
    n1 <- length(df[[i]][df$subset == T]) - sum(is.na(df[[i]][df$subset == T]))
    n2 <- length(df[[i]][df$subset == F]) - sum(is.na(df[[i]][df$subset == F]))
    
    x <- as.numeric(unlist(df[i])) # put variable {i} for completers into a temporary vector, make numeric.
    
    if(n1 <= 1 | n2 <= 1){
      print(paste("Not enough observations to compare ", i, ". This variable was skipped."))
    } else {
      wilcoxRaw <- wilcox.test(x ~ unlist(df[by]), alternative = c("two.sided"), paired = F) # run wilcoxon paired-rank
      wilcoxStat <- format(wilcoxRaw$statistic, scientific = T, digits = 3) # Extract W statistic of test. N.b., the 'W' statistic reported here is equivalent to the U statistic: https://stats.stackexchange.com/questions/79843/is-the-w-statistic-output-by-wilcox-test-in-r-the-same-as-the-u-statistic
      wilcoxP <- format(wilcoxRaw$p.value, scientific = T, digits = 3) # Extract P statistic of test
      r <- wilcoxonR(x = x, g = unlist(df[by]))
      wilcoxInterest <- c(wilcoxStat, wilcoxP, n1, n2, r) # combine W statistic, p, n1 and n2 into a vector
      wilcoxOutput[[i]] <- with(wilcoxOutput, wilcoxInterest) # Add the results for variable {i} to the output dataframe
      # print(wilcoxInterest)
    }
  }
  # return(wilcoxonOutput)
  write.csv(wilcoxOutput, file = glue("{outputDir}/{outputName}_{getDate()}.csv")) #exports wilcoxOutput to a CSV file
  print(paste("The Wilcoxon comparison has been saved as:", glue("{outputName}_{getDate()}.csv")))
} # "varsOfInterest" - list of variable names to compare between groups; "df" data frame with all vars of interest and grouping variable; "by" - grouping variable; "outputDir" - path to output directory; "outputName" - desired name of output file
### PCA analyses (no longer used) -----
PCACompute <- function(df, vars){
  # PCA (see 536 lecture notes week 10)
  rejectVars <- c()
  for(i in vars){
    if(is.numeric(as.matrix(df[[i]])) == FALSE){ # checks that variables are numeric
      rejectVars <- append(rejectVars, i)
    } else {
      mean <- mean(as.matrix(df[i]), na.rm = T)
      sd <- sd(as.matrix(df[i]), na.rm = T)
      
      if(round(mean,4) != 0 || round(sd,4) != 1){ # check if standardized
        df[i] <- scale(df[i])
        # print(c("Not standard.",colnames(df[i]), typeof(as.matrix(df[i]))))
        # print(c(mean(as.matrix(df[i]), na.rm = T), sd(as.matrix(df[i]), na.rm = T)))
      } else {
        # print(c("Yes.",colnames(df[i]), typeof(as.matrix(df[i]))))   
        # print(c(mean(as.matrix(df[i]), na.rm = T), sd(as.matrix(df[i]),na.rm = T)))
      }
    }
  }
  
  if(length(rejectVars) > 0){
    cat("The following variables are not numeric and cannot be included in PCA. \n\t")
    for(i in colnames(df[rejectVars])){
      cat(i, "\t")
    }
    vars <- vars[vars %in% rejectVars == FALSE] # removes the reject variables from the 'vars' list
  }
  
  out <- prcomp(as.matrix(df[,vars]), center = TRUE, scale = TRUE) # prints output
  return(out)
  
} # 'df' - dataframe; 'vars' - name of vars to conduct PCA on. N.B. these variables must be standardized before being entered here.

PCAPlot <- function(PCAOut, name){
  require(factoextra)
  
  screePlot <- fviz_eig(PCAOut,
                        choice = c("eigenvalue"),
                        geom = c("line"),
                        linecolor = "black",
                        ncp = 10,
                        addlabels = TRUE,
                        hjust = 0,
                        main = paste("Scree plot - ", name),
                        xlab = "Principal component",
                        ylab = NULL) +
    geom_hline(yintercept = 1, linetype = 2)
  
  return(screePlot)
  
} # 'PCAOut' - an object of class PCA (i.e., from the PCA compute function); 'name' - name of this PCA analysis

PCAApply <- function(df, PCAOut, componentsToRetain, prefix){
  
  varCols <- which(colnames(df) %in% names(PCAOut$rotation[,1]))
  df_z <- df %>%
    as_tibble() %>%
    mutate_at(varCols, ~(scale(.) %>% as.vector))
  
  for(i in 1:componentsToRetain){
    rowName <- paste(prefix, "PC", i, sep="_")
    df_z <- df_z %>%
      mutate("{{rowName}}" := as.matrix(df_z[varCols]) %*% PCAOut$rotation[,i])
  }
  
  return(df_z)
  
} # Gives each observation a score for each PCA component retained. 'df' - dataframe; 'PCAOut' - output from the PCACompute function; 'componentsToRetain' - an integer value specifying the number of principal components to retain (value must be between 1 and total number of components)); 'prefix' - a character string to be used as a prefix for the new column names.

PCASave <- function(PCAOut, fileName){
  fileName <- paste(fileName, "_", getDate(), ".csv", sep="")
  capture.output(PCAOut$rotation, file = fileName)
}

## 1.B Import data ----
### Remove EIDs who withdrew consent ----
df_all <- read_csv("./Data/Daniel_2022-05-02.csv") # import df
eidToRemove <- read_csv("./Data/eidToRemove-w45551_20220222.csv", col_names = "eid") # file with EID of participants who withdrew consent from UKBB
df_all <- filter_rows(df_all, !(eid %in% eidToRemove)) # remove rows with eid in this file
write.csv(df_all, glue("./Data/Daniel_{getDate()}_ProperEID.csv"), row.names=FALSE)
rm(df_all, eidToRemove) # remove dataframes not in use

### Import data for analysis ----
dataDate <- getDate()
# dataDate <- "2022-05-02"
df_all <- read_csv(glue("./Data/Daniel_{dataDate}_ProperEID.csv"), guess_max = 10000) # import df
df_all <- filter_rows(df_all, !(is.na(df_all$date_assess2_t2))) # remove rows without timepoint 2 assessment date
# df_small <- df_all[1:10000,] # for testing purposes
df <- df_all # specifies what primary df to use below.
rm(df_all) # remove dataframes not in use
# sort(colnames(df))
# View(df)

## 1.C Standardize categorical variables -----
df <- df %>%
  transform(sleep_duration0_t0 = as.character(sleep_duration0_t0)) %>%
  transform(sleep_duration2_t2 = as.character(sleep_duration2_t2))

df <- df %>%
  mutate(diet_cookedVeg0_t0 = case_when(  # for field: 1289 ('diet_cookedVeg'), -10 = Less than one; -1 = Do not know; -3 = Prefer not to answer;  
    diet_cookedVeg0_t0 == "Less than one" ~ "0.5",       
    diet_cookedVeg0_t0 == "Do not know" ~ "NA",       
    diet_cookedVeg0_t0 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_cookedVeg0_t0
  )) %>%
  mutate(diet_cookedVeg2_t2 = case_when(  # for field: 1289 ('diet_cookedVeg'), -10 = Less than one; -1 = Do not know; -3 = Prefer not to answer;  
    diet_cookedVeg2_t2 == "Less than one" ~ "0.5",       
    diet_cookedVeg2_t2 == "Do not know" ~ "NA",       
    diet_cookedVeg2_t2 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_cookedVeg2_t2
  )) %>%
  mutate(diet_fruit0_t0 = case_when( # for field: 1309 ('diet_fruit'), -10 = Less than one; -1 = Do not know; -3 = Prefer not to answer;  
    diet_fruit0_t0 == "Less than one" ~ "0.5",       
    diet_fruit0_t0 == "Do not know" ~ "NA",       
    diet_fruit0_t0 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_fruit0_t0
  )) %>%
  mutate(diet_fruit2_t2 = case_when( # for field: 1309 ('diet_fruit'), -10 = Less than one; -1 = Do not know; -3 = Prefer not to answer;  
    diet_fruit2_t2 == "Less than one" ~ "0.5",       
    diet_fruit2_t2 == "Do not know" ~ "NA",       
    diet_fruit2_t2 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_fruit2_t2
  )) %>%
  mutate(diet_rawVeg_t0 = case_when( # for field: 1299 ('diet_rawVeg'), -10 = Less than one; -1 = Do not know; -3 = Prefer not to answer;
    diet_rawVeg_t0 == "Less than one" ~ "0.5",       
    diet_rawVeg_t0 == "Do not know" ~ "NA",       
    diet_rawVeg_t0 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_rawVeg_t0
  )) %>%
  mutate(diet_processedMeat0_t0 = case_when( # field 1349 ('diet_processedMeat'), -1 = Do not know; -3 = Prefer not to answer;
    diet_processedMeat0_t0 == "Do not know" ~ "NA",       
    diet_processedMeat0_t0 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_processedMeat0_t0
  )) %>%
  mutate(diet_processedMeat2_t2 = case_when( # field 1349 ('diet_processedMeat'), -1 = Do not know; -3 = Prefer not to answer;
    diet_processedMeat2_t2 == "Do not know" ~ "NA",       
    diet_processedMeat2_t2 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_processedMeat2_t2
  )) %>%
  mutate(diet_water_t0 = case_when( # field 1528 ('diet_water'), -10 = Less than one; -1 = Do not know; -3 = Prefer not to answer;
    diet_water_t0 == "Less than one" ~ "0.5",       
    diet_water_t0 == "Do not know" ~ "NA",       
    diet_water_t0 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_water_t0
  )) %>%
  mutate(diet_alc_freq0_t0 = case_when( # field 1558 ('diet_alc_freq'), -3 = Prefer not to answer;
    diet_alc_freq0_t0 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_alc_freq0_t0
  )) %>%
  mutate(diet_alc_freq2_t2 = case_when( # field 1558 ('diet_alc_freq'), -3 = Prefer not to answer;
    diet_alc_freq2_t2 == "Prefer not to answer" ~ "NA",
    TRUE ~ diet_alc_freq2_t2
  )) %>%
  mutate(menopause0_t0 = case_when( # field 2724 ('menopause'), -3 = Prefer not to answer;
    menopause0_t0 == "Prefer not to answer" ~ "NA",
    TRUE ~ menopause0_t0
  )) %>%
  mutate(menopause2_t2 = case_when( # field 2724 ('menopause'), -3 = Prefer not to answer;
    menopause2_t2 == "Prefer not to answer" ~ "NA",
    TRUE ~ menopause2_t2
  )) %>%
  mutate(sleep_duration0_t0 = case_when( # field 1160 ('sleep_duration'), -1 = Do not know; -3 = Prefer not to answer;
    sleep_duration0_t0 == "Prefer not to answer" ~ "NA",
    sleep_duration0_t0 == "Do not know" ~ "NA",
    TRUE ~ sleep_duration0_t0
  )) %>%
  mutate(sleep_duration2_t2 = case_when( # field 1160 ('sleep_duration'), -1 = Do not know; -3 = Prefer not to answer;
    sleep_duration2_t2 == "Prefer not to answer" ~ "NA",
    sleep_duration2_t2 == "Do not know" ~ "NA",
    TRUE ~ sleep_duration2_t2
  )) %>%
  mutate(smoke_currently0_t0 = case_when( # field 2724 ('smoke_currently'), -3 = Prefer not to answer;
    smoke_currently0_t0 == "Prefer not to answer" ~ "NA",
    smoke_currently0_t0 == "No" ~ smoke_currently0_t0,
    TRUE ~ "Yes"
  )) %>%
  mutate(smoke_currently2_t2 = case_when( # field 2724 ('smoke_currently'), -3 = Prefer not to answer;
    smoke_currently2_t2 == "Prefer not to answer" ~ "NA",
    smoke_currently2_t2 == "No" ~ smoke_currently2_t2,
    TRUE ~ "Yes"
  )) %>%
  mutate(smoke_packYears_t0 = case_when( # field 2724 ('smoke_currently'), -3 = Prefer not to answer;which(colnames(df) == "dx_AnyOfInterest")
    smoke_ever_t0 == "No" ~ 0,
    TRUE ~ smoke_packYears_t0
  )) %>%
  mutate(cog_TMT_numericDuration_t2 =  case_when( # fields: 6348 ('cog_TMT_numericDuration'), 6350 ('cog_TMT_alphanumDuration'), 0 = trail not completed
    cog_TMT_numericDuration_t2 == "trail not completed" ~ "NA",
    TRUE ~ cog_TMT_numericDuration_t2
  )) %>%
  mutate(cog_TMT_alphanumDuration_t2 =  case_when(
    cog_TMT_alphanumDuration_t2 == "trail not completed" ~ "NA",
    TRUE ~ cog_TMT_alphanumDuration_t2
  )) %>%
  mutate(cog_numMem_maxDigitRemem_t2 = case_when( # field 4282 ('cog_numMem_maxDigitRemem_t2'), -1 = abandoned
    cog_numMem_maxDigitRemem_t2 == "abandoned" ~ "NA",
    TRUE ~ cog_numMem_maxDigitRemem_t2
  )) %>%
  mutate(demo_ethnicity_t0 = case_when( # field 4282 ('cog_numMem_maxDigitRemem_t2'), -1 = abandoned
    demo_ethnicity_t0 == "Prefer not to answer" ~ "NA",
    demo_ethnicity_t0 == "Do not know" ~ "NA",
    demo_ethnicity_t0 == "British" ~ "White",
    demo_ethnicity_t0 == "Irish" ~ "White",
    demo_ethnicity_t0 == "Any other white background" ~ "White",
    TRUE ~ "Non White"
  ))

## 1.D Compute additional columns ----
### CRP log transform
df <- df %>%
  mutate(crp_log = log(crp_aliquot_t0))

### Diagnosis - any of interest
dxVars <- starts_with("dx_", vars = colnames(df))

df <- df %>%
  mutate("dx_AnyOfInterest" = factor(case_when(
    if_any(any_of(dxVars), function(x) (x == "TRUE")) ~ 1,
    if_all(all_of(dxVars), function(x) (x == "FALSE")) ~ 0)))
dxVars <- c(dxVars, which(colnames(df) == "dx_AnyOfInterest"))

dxVarsNoCases <- c()
for(i in dxVars){
  if(n_distinct(df[i]) == 1){
    print(paste("No participant has dx ", colnames(df[i]), "(column ", i, " ); levels: ", levels(as.factor(df[i]))))
    dxVarsNoCases <- append(dxVarsNoCases, i)
  }
}
if(n_distinct(df[i]) == 1){
  for(i in dxVarsNoCases){
    dxVars <- dxVars[-which(dxVars == i)]
  }
}
# table(df$anyDx)

### Medication - any of interest
medVars <- starts_with("med_", vars = colnames(df))
# colnames(df[medVars])
medVars0 <- ends_with("t0", vars = colnames(df[medVars]))
medVars0 <- which(colnames(df) %in% colnames(df[medVars[medVars0]]))

medVars2 <- ends_with("t2", vars = colnames(df[medVars]))
medVars2 <- which(colnames(df) %in% colnames(df[medVars[medVars2]]))

medVars0_excludingSSRI <- medVars0[-which(colnames(df[medVars0]) == "med_SSRI_t0")]
medVars2_excludingSSRI <- medVars2[-which(colnames(df[medVars2]) == "med_SSRI_t2")]
df <- df %>%
  mutate("med_AnyOfInterest0_excludingSSRI" = factor(case_when(
    if_any(any_of(medVars0_excludingSSRI), function(x) (x == "TRUE")) ~ 1,
    if_all(all_of(medVars0_excludingSSRI), function(x) (x == "FALSE")) ~ 0))) %>%
  mutate("med_AnyOfInterest2_excludingSSRI" = factor(case_when(
    if_any(any_of(medVars2_excludingSSRI), function(x) (x == "TRUE")) ~ 1,
    if_all(all_of(medVars2_excludingSSRI), function(x) (x == "FALSE")) ~ 0)))
medVars0 <- c(medVars0, which(colnames(df) == "med_AnyOfInterest0_excludingSSRI"))
medVars2 <- c(medVars2, which(colnames(df) == "med_AnyOfInterest2_excludingSSRI"))
df <- df %>%
  mutate("med_AnyOfInterest0" = factor(case_when(
    if_any(any_of(medVars0), function(x) (x == "TRUE")) ~ 1,
    TRUE ~ 0))) %>%
  mutate("med_AnyOfInterest2" = factor(case_when(
    if_any(any_of(medVars2), function(x) (x == "TRUE")) ~ 1,
    TRUE ~ 0)))

medVars0 <- c(medVars0, which(colnames(df) == "med_AnyOfInterest0"))
medVars2 <- c(medVars2, which(colnames(df) == "med_AnyOfInterest2"))
medVars0NoCases <- c()
for(i in medVars0){
  if(n_distinct(df[i]) == 1){
    print(paste("No participant takes med ", colnames(df[i]), "; levels: ", levels(as.factor(df[i]))))
    medVars0NoCases <- append(medVars0NoCases, i)
  }
}

medVars2 <- c(medVars2, which(colnames(df) == "med_AnyOfInterest2"))
medVars2NoCases <- c()
for(i in medVars2){
  if(n_distinct(df[i]) == FALSE){
    print(paste("No participant takes med ", colnames(df[i]), "; levels: ", levels(as.factor(df[i]))))
    medVars2NoCases <- append(medVars2NoCases, i)
  }
}

if(is.null(medVars0NoCases) == FALSE){
  for(i in medVars0NoCases){
    medVars0 <- medVars0[-which(medVars0 == i)]
  }
}
if(is.null(medVars2NoCases) == FALSE){
  for(i in medVars2NoCases){
    medVars2 <- medVars2[-which(medVars2 == i)]
  }
}
# table(df$anyMed)

### Waist:hip ratio
df <- df %>%
  mutate(weight_waistToHip_t0 = weight_waistCirc0_t0/weight_hipCirc0_t0) %>%
  mutate(weight_waistToHip_t2 = weight_waistCirc2_t2/weight_hipCirc2_t2)

### Num days between assessment 0 and 2
df <- df %>%
  mutate(demo_daysBtwAssess = as.numeric(df$date_assess2_t2 - df$date_assess0_t0, units = "days"))

# saveRDS(df, file = glue("./Data/Processed/df_preLobes_{getDate()}.rds"))
# rm(df)

### Lobar brain variables
# dataDate <- "04_15_2023"
# df <- readRDS(glue("./Data/Processed/df_preLobes_{dataDate}.rds"))

vars_BrainRegion <- read_csv("  ./Variables/brainVarsbyLobe.csv")
vars_BrainRegion$varNameDf <- glue("{vars_BrainRegion$VariableName}_t2")

cat("\n Create lobar brain variable columns.")

for(metric_counter in unique(vars_BrainRegion$metric)){
  # cat("\n\t Metric: ", metric_counter)
  if(metric_counter %in% c("area", "mThick", "vol")){
    for(lobe in unique(vars_BrainRegion$AssociatedLobe)){
      # cat("\n\t\t Lobe: ", lobe)
      
        vars_lobe_metric_L <- vars_BrainRegion %>%
          filter(AssociatedLobe == lobe) %>%
          filter(metric == metric_counter) %>%
          filter(hemisphere == "L")
        vars_lobe_metric_R <- vars_BrainRegion %>%
          filter(AssociatedLobe == lobe) %>%
          filter(metric == metric_counter) %>%
          filter(hemisphere == "R")
        
        if(metric_counter == "area"){
          metric_lobe_L <- rowSums(df[vars_lobe_metric_L$varNameDf], na.rm = F)
          metric_lobe_R <- rowSums(df[vars_lobe_metric_R$varNameDf], na.rm = F)
          metric_lobe_WB <- rowSums(data.frame(metric_lobe_L, metric_lobe_R), na.rm = F)
        } else if (metric_counter == "mThick"){
          metric_lobe_L <- rowMeans(df[vars_lobe_metric_L$varNameDf], na.rm = F)
          metric_lobe_R <- rowMeans(df[vars_lobe_metric_R$varNameDf], na.rm = F)
          metric_lobe_WB <- rowMeans(data.frame(metric_lobe_L, metric_lobe_R), na.rm = F)
        } else if (metric_counter == "vol"){
          metric_lobe_L <- rowSums(df[vars_lobe_metric_L$varNameDf], na.rm = F)
          metric_lobe_R <- rowSums(df[vars_lobe_metric_R$varNameDf], na.rm = F)
          metric_lobe_WB <- rowSums(data.frame(metric_lobe_L, metric_lobe_R), na.rm = F)
        }
        
        tempMatrix <- as.matrix(cbind(metric_lobe_L,metric_lobe_R, metric_lobe_WB))
        colnames(tempMatrix) <- c(glue("{metric_counter}_{lobe}_L_t2"),glue("{metric_counter}_{lobe}_R_t2"),glue("{metric_counter}_{lobe}_WB_t2"))
        df <- cbind(df, as.data.frame(tempMatrix)) # creates columns in df for later computation of lobar brain metrics
        rm(vars_lobe_metric_L,vars_lobe_metric_R, metric_lobe_L, metric_lobe_R, metric_lobe_WB, tempMatrix)
      
      } 
    
    } else {
      cat("\n\t\t Error. The metric `", metric_counter, "` defined in the brain variables by brain lobe dataframe is not one of \'area\', \'vol\' or \'mthick\'. \n Skipping computation of variables of this metric.")
  }
}

### Hour of day of assessments
df <- df %>%
  mutate(crp_hourCollected = hour(df$crp_timeCollected_t0)) %>%
  mutate(cog_hourCompleted = hour(df$cog_timeCompleted_t2)) %>%
  mutate(brain_hourCompleted = hour(df$brain_timeCompleted_t2)) # Takes variables with form: 'YYYY-MM-DDTHH:MM:SS' and return the hour
df <- df %>% mutate(timeDif_brainHourMinusCogHourCompleted = df[124] - df[123]) # number of hours difference between cognitive assessment and MRI
hourColNums <- c(contains("_hourCompleted", vars = colnames(df)),contains("_hourCollected", vars = colnames(df)))

## 1.E Make columns appropriate data types ----
df <- readr::type_convert(df) # automatically detect
factorVars <- which(colnames(df) %in% c("demo_sex_t0", "demo_ethnicity_t0", "cog_prospMem_result_t2", "hand_t0","smoke_currently0_t0", "smoke_currently2_t2", "menopause0_t0", "menopause2_t2", "anyMedOfInterest0", "anyMedOfInterest2", "anyDxOfInterest", colnames(df[dxVars]), colnames(df[medVars0]), colnames(df[medVars2]))) # Factors
# sort(colnames(df[factorVars]))
df[,factorVars] <- lapply(df[,factorVars], factor)

df <- df %>%
  mutate(crp_aliquot_fctr = cut(crp_aliquot_t0, # CRP cutoffs (see Pearson et al., 2003)
                                breaks = c(0,1,3,10),
                                labels = c("Low", "Medium", "High"),
                                include.lowest = T,
                                right = F), .after = crp_aliquot_t0
  ) %>%
  mutate(cog_prospMem_result_t2 = case_when( # field 20018 ('cog_prospMem_result_t2'), 0 = instruction not recalled, either skipped or incorrect; 1 = correct recall on first attempt; 2 = correct recall on second attempt (as in Cullen et al., 2017)
    cog_prospMem_result_t2 == "correct recall on first attempt" ~ "1",
    cog_prospMem_result_t2 == "NA" ~ "NA",
    TRUE ~ "0"
  )) %>% mutate(cog_prospMem_result_t2  = as.numeric(cog_prospMem_result_t2))

df <- df %>%
  mutate(smoke_currently0_t0 = factor(smoke_currently0_t0, order = F, levels = c(
    "No",
    "Yes"))) %>%
  mutate(smoke_currently2_t2 = factor(smoke_currently2_t2, order = F, levels = c(
    "No",
    "Yes"))) %>%
  mutate(diet_processedMeat0_t0 = factor(diet_processedMeat0_t0, order = T, levels = c(
    "Never",
    "Less than once a week",
    "Once a week",
    "2-4 times a week",
    "5-6 times a week",
    "Once or more daily"))) %>%
  mutate(diet_processedMeat2_t2 = factor(diet_processedMeat2_t2, order = T, levels = c(
    "Never",
    "Less than once a week",
    "Once a week",
    "2-4 times a week",
    "5-6 times a week",
    "Once or more daily"))) %>%
  mutate(diet_alc_freq0_t0 = factor(diet_alc_freq0_t0, order = T, levels = c(
    "Never",
    "Special occasions only",
    "One to three times a month",
    "Once or twice a week",
    "Three or four times a week",
    "Daily or almost daily"))) %>%
  mutate(diet_alc_freq2_t2 = factor(diet_alc_freq2_t2, order = T, levels = c(
    "Never",
    "Special occasions only",
    "One to three times a month",
    "Once or twice a week",
    "Three or four times a week",
    "Daily or almost daily"))) %>%
  mutate(exercise_IPAQActivityGroup_t0 = factor(exercise_IPAQActivityGroup_t0, order = T, levels = c(
    "low",
    "moderate",
    "high"))) %>%
  mutate(ses_houseIncome0_t0 = factor(ses_houseIncome0_t0, order = T, levels = c(
    "Less than 18,000",
    "18,000 to 30,999",
    "31,000 to 51,999",
    "52,000 to 100,000",
    "Greater than 100,000"))) %>%
  mutate(ses_houseIncome2_t2 = factor(ses_houseIncome2_t2, order = T, levels = c(
    "Less than 18,000",
    "18,000 to 30,999",
    "31,000 to 51,999",
    "52,000 to 100,000",
    "Greater than 100,000"))) %>%
  mutate(med_SSRI_t0 = factor(case_when(
    med_SSRI_t0 == "FALSE" ~ 0,
    med_SSRI_t0 == "TRUE" ~ 1))) %>%
  mutate(med_SSRI_t2 = factor(case_when(
    med_SSRI_t2 == "FALSE" ~ 0,
    med_SSRI_t2 == "TRUE" ~ 1))) %>%
  mutate(med_Statin_t0 = factor(case_when(
    med_Statin_t0 == "FALSE" ~ 0,
    med_Statin_t0 == "TRUE" ~ 1))) %>%
  mutate(med_Statin_t2 = factor(case_when(
    med_Statin_t2 == "FALSE" ~ 0,
    med_Statin_t2 == "TRUE" ~ 1))) %>%
  mutate(med_Statin_t02 = factor(case_when(
    med_Statin_t0 == 1 ~ 1,
    med_Statin_t2 == 1 ~ 1, 
    TRUE ~ 0))) %>% # if using statin at either time point, will code as 1, else 0  
  mutate(med_Antihypertensive_t0 = factor(case_when(
    med_Antihypertensive_t0 == "FALSE" ~ 0,
    med_Antihypertensive_t0 == "TRUE" ~ 1))) %>%
  mutate(med_Antihypertensive_t2 = factor(case_when(
    med_Antihypertensive_t2 == "FALSE" ~ 0,
    med_Antihypertensive_t2 == "TRUE" ~ 1))) %>%
  mutate(med_Antihypertensive_t02 = factor(case_when(
    med_Antihypertensive_t0 == 1 ~ 1,
    med_Antihypertensive_t2 == 1 ~ 1, 
    TRUE ~ 0))) %>% # if using Antihypertensive at either time point, will code as 1, else 0
  mutate(med_AnyOfInterest0 = factor(med_AnyOfInterest0)) %>%
  mutate(med_AnyOfInterest0_excludingSSRI = factor(med_AnyOfInterest0_excludingSSRI)) %>%
  mutate(med_AnyOfInterest2 = factor(med_AnyOfInterest2)) %>%
  mutate(med_AnyOfInterest2_excludingSSRI = factor(med_AnyOfInterest2_excludingSSRI))

df$crp_aliquotdelay <- as.duration(interval(df$crp_timeCollected_t0, df$crp_assaydate_t0)) %/% as.duration(months(1))

df <- readr::type_convert(df) # automatically detect

## 1.F Take mean of covariates at instances 0 and 2 ----
df <- df %>%
  mutate(age_mean02 = (demo_age_assess0_t0 + demo_age_assess2_t2)/2) %>%
  mutate(weight_waistToHip_mean02 = (weight_waistToHip_t0 + weight_waistToHip_t2)/2) %>%
  mutate(sleep_duration_mean02 = (sleep_duration0_t0 + sleep_duration2_t2)/2)

## 1.G Save cleaned df -----
df_interim <- df
saveRDS(df_interim, file = glue("./Data/Processed/dfFormatted_full_{getDate()}.rds"))
rm(df_interim)

dataDate <- getDate()
# dataDate <- "12_14_2022" # For testing purposes. Format: MM_DD_YYYY
df <- readRDS(glue("./Data/Processed/dfFormatted_full_{dataDate}.rds"))

## 1.H Remove missing cases ----
### Define variables of interest ----
demoVars <- starts_with("demo_", vars = colnames(df))
sesVars <- starts_with("ses_", vars = colnames(df))
crpAliqVars <- which(colnames(df) %in% c("crp_aliquot_t0", "crp_log", "crp_aliquot_fctr"))
crpToZ <- which(colnames(df) %in% c("crp_aliquot_t0", "crp_log"))

brainPrefix <- c("area_","gmVol_", "GWcont_","mIntensity","mThick_", "vol_", "weigted_mFA_", "weigted_mICVF_", "weigted_mISOVF_", "weigted_mL1_", "weigted_mL2_", "weigted_mL3_", "weigted_mMD_", "weigted_mMO_", "weigted_mOD_")
brainMorphVars <- c()
for(i in brainPrefix){
  iVars <- starts_with(i, vars = colnames(df))
  # cat(paste(i, ":", length(iVars), "\n"), sep = "")
  brainMorphVars <- c(brainMorphVars, iVars)
}
# length(brainMorphVars)
# sort(colnames(df[brainMorphVars]))
# colnames(df[demoVars])

dietVars <- starts_with("diet_", vars = colnames(df))
# colnames(df[dietVars])
dietVarsToCor <- which(colnames(df) %in% c("diet_cookedVeg0_t0", "diet_cookedVeg2_t2", "diet_rawVeg_t0", "diet_fruit0_t0", "diet_fruit2_t2", "diet_processedMeat0_t0", "diet_processedMeat2_t2", "diet_water_t0"))

smokingVars <- starts_with("smoke_", vars = colnames(df))
weightVars <- c(starts_with("weight_", vars = colnames(df)), which(colnames(df) == "waistToHip"))
# colnames(df[weightVars])

exerciseVars <- starts_with("exercise_", vars = colnames(df))

completionTimeVars <- c(contains("timeCompleted", vars = colnames(df)), contains("timeCollected", vars = colnames(df)))
# colnames(df[completionTimeVars])
# cogPCAVars <- which(colnames(df) %in% c("cog_numMem_maxDigitRemem_t2", "cog_TMT_numericDuration_t2", "cog_TMT_numericErrors_t2", "cog_TMT_alphanumDuration_t2", "cog_TMT_alphanumErrors_t2", "cog_matrix_cor_t2", "cog_reactiontime_mean_t2", "cog_tower_cor_t2", "cog_digsub_cor_t2"))

cogVars <- starts_with("cog", vars = colnames(df))
cogVars_t2 <- cogVars[ends_with("_t2", vars = colnames(df[cogVars]))]
cogVarsToZ <- cogVars_t2[which(!colnames(df[c(cogVars_t2)]) %in% c("cog_prospMem_result_t2", "cog_timeCompleted_t2", "cog_hourCompleted"))]
# colnames(df[cogVarsToZ])

covarsToZ <- which(colnames(df) %in% c("demo_age_assess0_t0", "demo_age_assess2_t2", "age_mean02", "demo_daysBtwAssess", "ses_townsend_t0", "weight_waistToHip_t0", "weight_waistToHip_t2", "weight_waistToHip_mean02", "sleep_duration0_t0", "sleep_duration2_t2", "sleep_duration_mean02"))
varsToZ <- c(crpToZ, cogVarsToZ, covarsToZ)

df_varsToZ <- df[varsToZ] %>%
  mutate(across(.fns = scale, .names = "{colnames(df[varsToZ])}_z"))
# str(df_varsToZ)
vars_Z <- which(! colnames(df_varsToZ) %in% colnames(df[varsToZ])) # matrix of z-scores made in above line to add to main df\
df <- cbind(df, df_varsToZ[vars_Z])

df$cog_prospMem_result_t2 <- as.factor(df$cog_prospMem_result_t2)

finalCovars <- which(colnames(df) %in% c("demo_sex_t0", "demo_ethnicity_t0", "demo_age_assess0_t0_z", "demo_daysBtwAssess_z", "ses_townsend_t0_z", "weight_waistToHip_mean02_z", "sleep_duration_mean02_z", "exercise_IPAQActivityGroup_t0", "med_Antihypertensive_t02", "med_Statin_t02"))
# sort(colnames(df))
# colnames(df[sesVars])
# View(df[1:100, c(1,cogPCAVars)])

### Final list of covariates ----
essentialVars <- c(finalCovars, completionTimeVars, crpAliqVars, brainMorphVars) # list of variable names that must be complete in order for case to be retained
# essentialVars <- readRDS("revisedEssentialVarsAnalysis.rds")

### Remove cases missing data in essential variables ---- 
df <- df %>%
  mutate(missingEssentialVar = case_when(
    if_any(.cols = essentialVars, function(x)(is.na(x))) ~ TRUE,
    TRUE ~ FALSE)) %>%
  mutate(missingEssentialVar = factor(missingEssentialVar))

saveRDS(df, file = glue("./Data/Processed/df_PreSubset_{getDate()}.rds"))
# sort(colnames(df))
# table(df$missingEssentialVar)
# table(df$missingEssentialVarsanalysis)

cat("df_removed: ")
df_removed <- filter_rows(df, missingEssentialVar == TRUE)
# df_removed <- df %>%
#   filter(if_any(any_of(essentialVars), function(x)(is.na(x)))) # make new df of cases with missing values

df_retained <- filter_rows(df, missingEssentialVar == FALSE) # remove rows of those with missing values
# df_retained <- df %>%
#   filter(!(eid %in% df_removed$eid)) # remove rows of those with missing values

### Detail number of missing values per variable ----
numMissing <- matrix(ncol=2, nrow = length(essentialVars), dimnames = list(colnames(df_removed[essentialVars]), c("NACount", "NotMissingOnThisVar"))) # initialize df to hold number of missing cases per variable 

for(i in 1:length(essentialVars)){
  numMissing[i,] <- c(sum(is.na(df_removed[essentialVars[i]])), sum(!is.na(df_removed[essentialVars[i]])))
}

numMissing_df <- as.data.frame(numMissing) %>%  arrange(desc(NACount))

# kable(numMissing_df, caption = paste("Number of missing cases per variable. Number of rows with at least one missing value: ", nrow(df_removed), " out of ", sum(df$missingEssentialVar == FALSE)+nrow(df_removed), " total cases (", round(nrow(df_removed)/(sum(df$missingEssentialVar == FALSE)+nrow(df_removed))*100, 2), "%).", sep = ""))
write.csv(numMissing_df, file = glue("./outputs/descriptive/sourceOfNAValues_{getDate()}.csv"))

### Remove cases with CRP > 10 ----
# hist(df$crp_aliquot_t0)
# hist(df$crp_aliquot_t0[df$crp_aliquot_t0 > 10])
# max(df_retained$crp_aliquot_t0)
df_retained <- filter_rows(df_retained, crp_aliquot_t0 <= 10) # remove CRP > 10

### Save and describe retained cases ----
saveRDS(df_retained, file = glue("./Data/Processed/dfRetained_{getDate()}.rds"))
write.csv(df_retained, file = glue("./Data/Processed/dfRetained_{getDate()}.csv"))
write.csv(psych::describe(df), file = glue("./outputs/descriptive/describedRetainedCases_{getDate()}.csv"))

for(i in c(medVars0NoCases, medVars2NoCases, dxVarsNoCases)){
  if(is.null(i) == FALSE){
    if(i %in% varsToSummarise){
      colNum <- which(colnames(df) == colnames(df[i]))
      varsToSummarise <- varsToSummarise[-which(varsToSummarise == colNum)]
    }
  }
}

varsToSummarise <- c(1:ncol(df_retained))
numVars<- lapply(df_retained, numNotFactor)
df_numVars <- as.data.frame(do.call(cbind, numVars))
numVarCols <- (which(colnames(df_retained) %in% colnames(df_numVars)))
numVars <- c((which(colnames(df_retained) %in% colnames(df_numVars))), which(colnames(df_retained) == "date_assess0_t0"), which(colnames(df_retained) == "date_assess2_t2"), which(colnames(df_retained) == "cog_timeCompleted_t2"), which(colnames(df_retained) == "brain_timeCompleted_t2"), which(colnames(df_retained) == "crp_timeCollected_t0"), which(colnames(df_retained) == "crp_assaydate_t0"))

# write.csv(describeBy(df[numVars], group = df$crp_aliquot_fctr)[1], file = paste("dfRetained_Summary_byCRP_numVars_", names(describeBy(df, group = df$crp_aliquot_fctr)[1]), "_", getDate(), ".csv", sep = "")) # save counts of categorical variables
# write.csv(summary(df[-numVars]), file = paste("dfRetained_Summary_fctrVars", getDate(), ".csv", sep = "")) # save counts of categorical variables
# 
# write.csv(describeBy(df[numVars], group = df$crp_aliquot_fctr)[1], file = paste("dfRetained_Summary_byCRP_numVars_", names(describeBy(df, group = df$crp_aliquot_fctr)[1]), "_", getDate(), ".csv", sep = "")) # save counts of categorical variables
# write.csv(describeBy(df[numVars], group = df$crp_aliquot_fctr)[2], file = paste("dfRetained_Summary_byCRP_numVars_", names(describeBy(df, group = df$crp_aliquot_fctr)[2]), "_", getDate(), ".csv", sep = "")) # save counts of categorical variables
# write.csv(describeBy(df[numVars], group = df$crp_aliquot_fctr)[3], file = paste("dfRetained_Summary_byCRP_numVars_", names(describeBy(df, group = df$crp_aliquot_fctr)[3]), "_", getDate(), ".csv", sep = "")) # save counts of categorical variables

## 1.I Clean environment to avoid crashes -----
itemsToKeep <- c("df", "df_retained", "getDate", "filter_rows", "numNotFactor", "runWilcox", "medVars0NoCases", "medVars2NoCases", "dxVarsNoCases", "normalCheck")
itemstoRemove <- ls()[c(ls() %in% itemsToKeep) == FALSE]
rm(list = itemstoRemove)

## 1.J Create subsets ----
dataDate <- getDate()
# dataDate <- "12_14_2022"

df <- readRDS(glue("./Data/Processed/df_PreSubset_{dataDate}.rds"))
df_retained <- readRDS(glue("./Data/Processed/dfRetained_{dataDate}.rds"))

subsetData_outputPath <- "./Data/Processed/Subsets"
subsetNames <- c("noMedNoDx", "All", "oldest","onlyDx","noDx","onlyMed","noMed","onlySSRI","youngest","noDxNoSSRI","noSSRI")

df <- df %>%
  mutate(subset_all = factor(case_when(
    missingEssentialVar == TRUE ~ FALSE,
    TRUE ~ TRUE))) %>% 
  mutate(subset_oldest = factor(case_when(
    missingEssentialVar == TRUE ~ FALSE,
    demo_age_assess0_t0 >= quantile(df_retained$demo_age_assess0_t0, c(.66))  ~ TRUE,
    TRUE ~ FALSE))) %>% 
  mutate(subset_youngest = factor(case_when(
     missingEssentialVar == TRUE ~ FALSE,
     demo_age_assess0_t0 <= quantile(df_retained$demo_age_assess0_t0, c(.33))  ~ TRUE,
     TRUE ~ FALSE)))` %>% 
  mutate(subset_noMed = factor(case_when(
     missingEssentialVar == TRUE ~ FALSE,
     med_AnyOfInterest0 == "1" ~ FALSE,
     med_AnyOfInterest2 == "1" ~ FALSE,
     TRUE ~ TRUE))) %>%
  mutate(subset_onlyMed = factor(case_when(
     missingEssentialVar == TRUE ~ FALSE,
     med_AnyOfInterest0 == "1" ~ TRUE,
     med_AnyOfInterest2 == "1" ~ TRUE,
     TRUE ~ FALSE))) %>%
  mutate(subset_noSSRI = factor(case_when(
     missingEssentialVar == TRUE ~ FALSE,
     med_SSRI_t0 == "1" ~ FALSE, 
     med_SSRI_t2 == "1" ~ FALSE,
     TRUE ~ TRUE))) %>% 
  mutate(subset_onlySSRI = factor(case_when(
     missingEssentialVar == TRUE ~ FALSE,
     med_SSRI_t0 == "1" ~ TRUE,
     med_SSRI_t2 == "1" ~ TRUE,
     TRUE ~ FALSE))) %>% 
  mutate(subset_noDx = factor(case_when(
     missingEssentialVar == TRUE ~ FALSE,
     dx_AnyOfInterest == "0" ~ TRUE,
     TRUE ~ FALSE))) %>% 
  mutate(subset_onlyDx = factor(case_when(
     missingEssentialVar == TRUE ~ FALSE,
     dx_AnyOfInterest == "1" ~ TRUE,
     TRUE ~ FALSE))) %>% 
  mutate(subset_noMedNoDx = factor(case_when(
     missingEssentialVar == TRUE ~ FALSE,
     dx_AnyOfInterest == "1" ~ FALSE, 
     med_AnyOfInterest0 == "1" ~ FALSE,
     med_AnyOfInterest2 == "1" ~ FALSE,
     TRUE ~ TRUE))) %>%
  mutate(subset_noDxNoSSRI = factor(case_when(
     missingEssentialVar == TRUE ~ FALSE,
     dx_AnyOfInterest == "1" ~ FALSE,
     med_SSRI_t0 == "1" ~ FALSE,
     med_SSRI_t2 == "1" ~ FALSE,
     TRUE ~ TRUE))) # create one variable per subset to identify case's belonging to each subset

# table(df$subset_all)
# table(df$subset_oldest)
# table(df$subset_youngest)
# table(df$subset_noMed)
# table(df$subset_onlyMed)
# table(df$subset_noSSRI)
# table(df$subset_onlySSRI)
# table(df$subset_noDx)
# table(df$subset_onlyDx)
# table(df$subset_noMedNoDx)
# table(df$subset_noDxNoSSRI)

write.csv(df, file = glue("./Data/Processed/df_IndicatedSubsets_{getDate()}.csv")) # exports cleaned dataframe with all variables ready for analysis
saveRDS(df, file = glue("./Data/Processed/df_IndicatedSubsets_{getDate()}.rds"))

# Save
# subsetNames <- c("noMedNoDx")

for(subset in subsetNames){
  cat("Saving subset: ", subset, " \n\t")
  subsetColName <- glue("subset_{subset}")
  subsetCol <- ends_with(subsetColName, vars = colnames(df))
  df_subset <- filter_rows(df, df[subsetCol] == "TRUE") # filter by column identifying subset belonging
  write.csv(df_subset, file = glue("{subsetData_outputPath}/{subset}_{getDate()}.csv")) # exports cleaned dataframe with all variables ready for analysis
  saveRDS(df_subset, file = glue("{subsetData_outputPath}/{subset}_{getDate()}.rds"))
}

## 1.K Summarise variables ----
dataDate <- getDate()
# dataDate <- "12_14_2022"

df <- readRDS(glue("./Data/Processed/dfFormatted_full_{dataDate}.rds"))

dataPath <- "./Data/Processed/Subsets"
# subsetNames <- "noDxNoSSRI"
dataDate <- getDate()
# dataDate <- "12_14_2022" # MM_DD_YYYY

genDescriptionOutputPath <- "./outputs/descriptive"
nonParaDescriptionOutputPath <- "./outputs/descriptive"
WilcoxOutputPath <- "./outputs/descriptive/wilcox"

nonParaVars_names <- c("ses_townsend_t0", "weight_BMI0_t0", "weight_BMI2_t2", "sleep_duration0_t0", "sleep_duration2_t2", "diet_alcohol_yesterdayIntake0_t0", "diet_alcohol_yesterdayIntake2_t2", "crp_aliquot_t0")

for(subset in subsetNames){
  
  fileName <- glue("{dataPath}/{subset}_{dataDate}.rds")
  subset_df <- readRDS(fileName)
  print(subset)
  
  write.csv(psych::describe(subset_df), file = glue("{genDescriptionOutputPath}/described_{subset}_{getDate()}.csv")) # exports description file
  
  df <- df %>%  mutate(subset = case_when(
    df$eid %in% subset_df$eid ~ T, # Find cases in this subset
    TRUE ~ F   # Find complementary set of cases to the subset
  )) # Create new column in 'df' indicating if the case is in this subset
  
  if(all.equal(df$eid[df$subset == T], subset_df$eid) != TRUE){
    cat("Error. Unable to find subset cases in the complete unsubsetted data. Skipping wilcox test.")
  } else {
    runWilcox(varsOfInterest = nonParaVars_names, df = df, by = "subset", outputDir = WilcoxOutputPath, outputName = glue("{subset}_Wilcox"))# wilcoxOutput to a CSV file
  }
  
  nonParaVars_summary <- TOne(as.data.frame(df[nonParaVars_names]), grp = df$subset, colnames = c("Excluded", subset), add.length = T, total = T,
                              FUN = function(x) gettextf("%s (%s)",
                                                         Format(median(x, na.rm = T), digits = 2, big.mark = ","),
                                                         Format(IQR(x, na.rm = T), digits = 2, big.mark = ",")),
                              TEST = NA)
  write.csv(nonParaVars_summary, file = glue("{nonParaDescriptionOutputPath}/{subset}_nonParaSummary_{getDate()}.csv"))
  print(paste("The description of non-normal variables of interest has been saved as:", glue("{subset}_nonParaSummary_{getDate()}.csv")))
}
rm(subset_df)
## 1.L Compare removed cases to retained ----
varsToCompare <- c("demo_age_assess0_t0","demo_age_assess2_t2","demo_sex_t0","demo_ethnicity_t0","ses_townsend_t0","weight_waistToHip_t0","weight_BMI0_t0","weight_waistToHip_t2","weight_BMI2_t2","sleep_duration0_t0","sleep_duration2_t2","diet_alcohol_yesterdayIntake0_t0","diet_alcohol_yesterdayIntake2_t2","crp_aliquot_t0","crp_log","dx_AnyOfInterest","dx_Dementia","dx_SSD","dx_MoodDisorder","med_AnyOfInterest0","med_Antihypertensive_t0","med_SSRI_t0","med_Statin_t0","med_AnyOfInterest2","med_Antihypertensive_t2","med_SSRI_t2","med_Statin_t2","exercise_IPAQActivityGroup_t0","smoke_currently0_t0","smoke_currently2_t2","menopause0_t0","menopause2_t2","hand_t0","demo_birthYear_t0","age_mean02","sleep_duration_mean02","diet_alc_freq0_t0","diet_alc_freq2_t2","crp_hourCollected","crp_timeFasting_t0","crp_aliquotdelay","med_AnyOfInterest0_excludingSSRI","med_AnyOfInterest2_excludingSSRI","med_Antihypertensive_t02", "med_Statin_t02", "weight_waistToHip_mean02","diet_processedMeat0_t0","diet_processedMeat2_t2","diet_cookedVeg0_t0","diet_cookedVeg2_t2","diet_rawVeg_t0","diet_fruit0_t0","diet_fruit2_t2","diet_water_t0")

varsToCompare_ColNum <- which(colnames(df) %in% varsToCompare)

for(i in c(medVars0NoCases, medVars2NoCases, dxVarsNoCases)){
  if(is.null(i) == FALSE){
    if(i %in% varsToCompare_ColNum){
      colNum <- which(colnames(df) == colnames(df[i]))
      varsToCompare_ColNum <- varsToCompare_ColNum[-which(varsToCompare_ColNum == colNum)]
    }
  }
}

#### Define non-factor numeric variables for normality check ----
numVarsToCompare <- lapply(df[varsToCompare_ColNum], numNotFactor)
df_numVarsToCompare <- as.data.frame(do.call(cbind, numVarsToCompare))
numVarColsToCompare <- (which(colnames(df) %in% colnames(df_numVarsToCompare)))

#### Define non-normal variables for non-parametric descriptive analyses ----
normalityPlots <- mapply(normalCheck, df[c(numVarColsToCompare)], colName = colnames(df[c(numVarColsToCompare)]), colNum = c(numVarColsToCompare))
cat("Scroll through the plots to determine normality of variables. \n List non-normal variables on line 863.")
# N.b. QQ plot should follow a straight line
nonParaToCompare <- c(3,10,11,14,15,18,19,20,21,22,25,47,48,49,50,713,719,720,757,763)

#### Make Table 1 for those excluded based on missing variables
summaryMissing <- list()

# subsetNames <- "noDxNoSSRI"
for(subset in subsetNames){
  groupingVar <- glue("subset_{subset}")
  groupingVarCol <- which(colnames(df) == groupingVar)
  comparisonFormula <- glue("x ~ {groupingVar}")
  comparisonDroppedCases <- TOne(df[varsToCompare], grp = df[[groupingVarCol]], add.length = T, total = T,
                                 FUN = function(x) gettextf("%s (%s)",
                                                            Format(mean(x, na.rm = T), digits = 2, big.mark = ","),
                                                            Format(sd(x, na.rm = T), digits = 2, big.mark = ",")),
                                 TEST = list(
                                   num  = list(fun = function(x, g){paste(
                                     "F(",
                                     summary(aov(x ~ g))[[1]][1, "Df"], ",",
                                     summary(aov(x ~ g))[[1]][2, "Df"], ") = ",
                                     formatC(summary(aov(x ~ g))[[1]][1, "F value"], digits = 5),
                                     ", p = ", formatC(summary(aov(x ~ g))[[1]][1, "Pr(>F)"], format = "f", digits = 3),
                                     ", g = ", formatC(effsize::cohen.d(x ~ g, data = df, pooled = T, hedges.correction = T)$estimate, format = "f", digits = 3), # calculate g. See Durlak 2009, J. Ped. Psyc.
                                     sep = "")},
                                     lbl = "ANOVA"),
                                   cat  = list(fun = function(x, g){paste(
                                     "χ²(", chisq.test(table(x, g))$parameter, ") = ",
                                     formatC(chisq.test(table(x, g))$statistic, digits = 5),
                                     ", p = ", formatC(chisq.test(table(x, g))$p.val, format = "f", digits = 3),
                                     ", v = ", formatC(sqrt(chisq.test(table(x, g))$statistic)/(length(x)*chisq.test(table(x, g))$parameter), format = "g", digits = 2),
                                     sep = "")},
                                     lbl = "Chi-Square test"),
                                   dich = list(fun = function(x, g){paste(
                                     "χ²(", chisq.test(table(x, g))$parameter, ") =",
                                     formatC(chisq.test(table(x, g))$statistic, digits = 2),
                                     ", p = ",  formatC(chisq.test(table(x, g))$p.val, format = "f", digits = 3),
                                     ", v = ", formatC(sqrt(chisq.test(table(x, g))$statistic)/(length(x)*chisq.test(table(x, g))$parameter), format = "g", digits = 2),
                                     sep = "")},
                                     lbl = "Chi-Square test")),
                                 fmt = list(abs  = Fmt("abs"),
                                            num  = Fmt("num"),
                                            per  = Fmt("per"),
                                            pval = as.fmt(fmt = "p*")))
  write.csv(comparisonDroppedCases, file = glue("./outputs/descriptive/TOne_{subset}_{getDate()}.csv"))
}

## 1.M Preliminary correlation analyses ----
### CRP and blood draw characteristics  ----
cat("Hour collected: \n\t median:", median(hour(df$crp_timeCollected_t0), na.rm = T), "\n \t IQR: \t", IQR(hour(df$crp_timeCollected_t0), na.rm = T))
cat("Correlation between CRP and hour collected: ")
cor.test(rank(df$crp_aliquot_t0), rank(hour(df$crp_timeCollected_t0)), na.action = na.omit, conf.level = .95)

cat("Fasting time: \n\t median:", median(df$crp_timeFasting_t0, na.rm = T),"\n \t IQR: \t", IQR(df$crp_timeFasting_t0))
cat("Correlation between CRP and hours fasted before blood draw: ")
cor.test(rank(df$crp_aliquot_t0), rank(df$crp_timeFasting_t0), na.action = na.omit)

cat("CRP aliquots were performed between", paste(min(df$crp_assaydate_t0))," and ", paste(max(df$crp_assaydate_t0)))
df$crp_assaydate_t0 <- as.POSIXct(df$crp_assaydate_t0)

cat("the mean time between blood sample and CRP measurement was (in years): \n\t", "mean:", mean(df$crp_aliquotdelay, na.rm = T), "SD: ", sd(df$crp_aliquotdelay, na.rm = T), "\n\t Median: ", median(df$crp_aliquotdelay, na.rm = T), "IQR: ", IQR(df$crp_aliquotdelay))
cor.test(rank(df$crp_aliquotdelay), rank(df$crp_aliquot_t0))

# 2. Statistical Analyses -----
## 2.A Correlation ----
### 2.A.i. Functions ----
createLMExpression <- function(varNames, X, Y, M, C){
  colNumsPredictors <- c()
  for(i in c(X, M, C)){
    if(i != "ignore"){
      colNum <- which(varNames == i)
      # print(i)
      colNumsPredictors <- append(colNumsPredictors, colNum)
    }
  } # get column numbers for the relevant predictor variables
  
  LMpredictors <- paste(sapply(varNames[colNumsPredictors], paste), collapse = " + ")
  model <- paste(Y, " ~ ", LMpredictors, sep = "")
  return(model)
} # Input: 'varNames' - list of variable names in the dataframe to be analysed. Is given by the function colnames(df); 'Y' - name of the dependent variable; 'X' - name of the independent variable; 'M' - name of mediators mediators; 'C' - list of covariate variable names;

runCor <- function(df, X, Y, M, covars, contrasts, YisFactor, name){
  expression <- createLMExpression(varNames = colnames(df), X = X, Y = Y, M = M, C = covars) # total model: predict Y from X and C  
  
  if(YisFactor){
    model <- glm(formula = eval(parse(text = expression)), data = df, contrasts = contrasts, family = "binomial", na.action = na.omit)
  } else {
    model <- lm(eval(parse(text = expression)), data = df, contrasts = contrasts, na.action = na.omit)
  }
  
  capture.output(
    cat(" -- Linear Model Analysis Output -- "),
    cat("\n ------------------ \n Model expression: \n\t", expression),
    print(summary(model)),
    file = glue("{name}_{getDate()}.txt"))
  
  return(model)
}

corTxtFormat <- function(filePath, files, X, listOfCovars, logReg){  
  require(stringr)
  
  df_main <- matrix(ncol = 13, dimnames = list(c(), c("subset", "Y", "X", "b", "SE", "t", "p_predictor", "F", "df1", "df2", "p_model", "R2", "R2adj")))
  
  df_covarForModel_tmp <- c()
  
  for(file in files){
    specificFilePath <- glue("{filePath}/{file}.txt")
    fileObj <- readLines(specificFilePath)
    
    date <- str_match(file, "t2_z_(.*2)")[2]
    subsetName <- gsub("_","",str_match(file, "cor_*(.*?)_cog"))[2]
    Y <- str_match(file, "cog_*(.*?)_t2_z")[,2]
    modelDetails <- c(subsetName, Y)
    
    cat("\n\t", glue("{Y}"))
    
    # Results of interest for model
    effectSizes <- grep("Multiple R-squared",fileObj, value = TRUE)
    effectSizes <- unlist(strsplit(effectSizes, "\\s+"))
    R2 <- substr(effectSizes[3], 1, nchar(effectSizes[3])-1)
    R2adj <- effectSizes[6]
    
    FTest <- grep("F-statistic",fileObj, value = TRUE)
    FTest <- unlist(strsplit(FTest, "\\s+"))
    FStat <- FTest[2]
    df1 <- FTest[4]
    df2 <- FTest[6]
    
    p_model <- FTest[which(FTest == grep("p-value",FTest, value = TRUE)):length(FTest)]
    # print(p_model)
    if(TRUE %in% grepl("<", p_model) && nchar(p_model[grep("<",p_model)]) == 1){
      p_model <- glue("{p_model[[2]]}{p_model[[3]]}")
    } else {
      p_model <- p_model[2]
    }
    
    modelResults <- c(FStat, df1, df2, p_model, R2, R2adj)
    
    for(predictor in c(X,listOfCovars)){
      result <- grep(predictor, fileObj, value = T)
      results_seperated <- unlist(strsplit(result, "\\s+"))
      # print(results_seperated)
      # results of interest for each predictor
      b <- results_seperated[2]
      SE <-results_seperated[3] 
      t <- results_seperated[4]
      if("<" %in% results_seperated == TRUE){
        p_predictor <- glue("{results_seperated[[5]]}{results_seperated[[6]]}")
      } else {
        p_predictor <- results_seperated[5]
      }
      
      df_main <- rbind(df_main, c(modelDetails, predictor, b, SE, t, p_predictor, modelResults))
      # cat(paste(name, depVarValuesNoSigChar, sep = "\n\t"))
    } 
  }
  
  df_main <- df_main[-1,]
  write.csv(df_main, file = glue("{subsetName}_correlSummary_{date}.csv"))
} # 'filePath' - path to folder housing files, 'files' - name of files (without extensions), 'x' - independent variable (usually 'CRP_log_z'), 'listOfCovars' - list of covariate variable names as listed in LM output (include level of contrast for categorical variables). If no covars, define as '""'. 'logReg': logical specifying if file is logistic regression. N.B. takes subset name and dependent variable from file name.

corCSVFormat <- function(modelName, corResultPath, dataFilePrefix, subsets, cogVars, dataDate, depVar, listOfCovars, outputPrefix, outputPath){  
  df_main <- matrix(ncol = 15, dimnames = list(c(), c("subset",  "X","Y", "b", "b_SE", "b_95CI_Lo", "b_95CI_hi", "b_t", "b_p", "R2", "R2adj", "F", "df1", "df2", "p")))
  df_covarForModel <- matrix(ncol = length(listOfCovars)+2)
  df_covarForModel_tmp <- c()
  
  for(subset in subsets){
    for(cogVar in cogVars){
      name <- glue("{subset}_{cogVar}")
      fileName <- glue("{dataFilePrefix}_{name}_{dataDate}.txt")
      # cat(fileName, "\n")
      if(fileName %in% list.files(corResultPath)){
        fileObj <- readLines(glue("{corResultPath}/{fileName}"))
        
        depVarValues <- grep(depVar, fileObj, value = TRUE)
        depVarValues <- unlist(strsplit(depVarValues, "\\s+"))

        if(TRUE %in% grepl("<", depVarValues)){
          pAndSymbol <- glue("{depVarValues[[5]]}{depVarValues[[6]]}")
          depVarValuesNoSigChar <- c(depVarValues[1:4], pAndSymbol)
        } else {
          depVarValuesNoSigChar <- depVarValues[1:5]
        }
        # cat(paste(name, depVarValuesNoSigChar, sep = "\n\t"))

        CI_lo <- as.numeric(depVarValuesNoSigChar[2]) - abs(qnorm(.025))*as.numeric(depVarValuesNoSigChar[3]) # b - 1.96*SE
        CI_hi <- as.numeric(depVarValuesNoSigChar[2]) + abs(qnorm(.025))*as.numeric(depVarValuesNoSigChar[3]) # b + 1.96*SE
        # cat("95%CI: ", CI_lo, CI_hi)
        
        effectSizes <- grep("Multiple R-squared",fileObj, value = TRUE)
        effectSizes <- unlist(strsplit(effectSizes, "\\s+"))
        R2 <- substr(effectSizes[3], 1, nchar(effectSizes[3])-1)
        R2adj <- effectSizes[6]
        
        FTest <- grep("F-statistic",fileObj, value = TRUE)
        FTest <- unlist(strsplit(FTest, "\\s+"))
        FStat <- FTest[2]
        df1 <- FTest[4]
        df2 <- FTest[6]
        p <- FTest[which(FTest == grep("p-value",FTest, value = TRUE)):length(FTest)]
        
        if(TRUE %in% grepl("<", p)){
          p <- glue("{FTest[[9]]}{FTest[[10]]}")
          # p <- paste(FTest[[9]], FTest[[10]], sep = "")
        } else {
          p <- FTest[which(FTest == grep("p-value",FTest, value = TRUE))+1]
        }
        
        df_mainRow <- c(subset, depVarValuesNoSigChar[1], cogVar, depVarValuesNoSigChar[2:3], CI_lo, CI_hi, depVarValuesNoSigChar[4:length(depVarValuesNoSigChar)], R2, R2adj, FStat, df1, df2, p)  # "Model name", "indep var", "dep var", "parameter estimate", "SE", "t", "p", "R2", "R2adj", "F", "df1", "df2", "p"
        df_main <- rbind(df_main, " " = df_mainRow)
        
      } else {
        cat("\n Warning. The file ", fileName, " is not in the provided path. This analysis is not added to the summary file.")
      }
    }
  }
  df_main <- df_main[-1,]
  write.csv(df_main, file = glue("{outputPath}/{outputPrefix}_{getDate()}.csv"))
} # 'modelName' - name of model for labeling output rows (usually 'cor' or 'med'), 'corResultPath' - path to folder housing correlation text files, 'dataFilePrefix' - prefix of correlation files, 'subsets' - list of subset names, 'cogVars' - list of cognitive variables correlation analyses performed on, 'depVar' - dependent variable (usually CRP_log_z), 'dataDate' - date of correlational analysis (as in file name), 'listOfCovars' - list of covariate variable names as listed in LM output (include level of contrast for categorical variables). If no covars, define as '""', 'outputPrefix' - string indicating the prefix for the output file. Note, independent variable is assumed to be the variable beginning with "cog_", 'outputPath' - string defining path to save summary file.
### 2.A.ii. Define variables of interest ----
# sort(colnames(df))
crp_fctrCol <- which(colnames(df) == "crp_aliquot_fctr")
crp_log <-  which(colnames(df) == "crp_log")
crp_log_z <- which(colnames(df) == "crp_log_z")
crp <- c(which(colnames(df) == "crp_aliquot_t0"), crp_fctrCol, crp_log, crp_log_z)

cogVars <- starts_with("cog", vars = colnames(df))
cogVars_t2_z <- cogVars[ends_with("_t2_z", vars = colnames(df[cogVars]))]
cogOutcomes_colNum <- c(cogVars_t2_z)
cogOutcomes_colNames <- colnames(df[cogOutcomes_colNum]) # list of all numeric cognitive variables to assess
proceduralCogVarsToRemoveFromAnal <- c("cog_fluidIntel_QsAttempted_t2_z", "cog_matrix_numViewed_t2_z", "cog_digsub_numAttempted_t2_z", "cog_prospMem_timeDelay_t2_z")
cogOutcomes_colNames <- cogOutcomes_colNames[!cogOutcomes_colNames %in% proceduralCogVarsToRemoveFromAnal]
prospMemCol <- which(colnames(df) == "cog_prospMem_result_t2")
cogOutcomes_fctr <- c(prospMemCol)
cogOutcomes_fctr_names <- colnames(df[cogOutcomes_fctr])

df$cog_hourCompleted_z <- scale(df$cog_hourCompleted)
covars_colNames <- c("demo_sex_t0", "demo_ethnicity_t0", "demo_age_assess0_t0_z", "demo_daysBtwAssess_z", "ses_townsend_t0_z", "weight_waistToHip_mean02_z", "sleep_duration_mean02_z", "smoke_currently0_t0", "smoke_currently2_t2","exercise_IPAQActivityGroup_t0", "med_Antihypertensive_t02", "med_Statin_t02")
covars_colNum <-  which(colnames(df) %in% covars_colNames)

covars_brain_colNames <- c("hand_t0", "brain_headScale_t2")
covars_brain_colNum <- which(colnames(df) %in% covars_brain_colNames)

essentialVarsAnalysis_colNum <- c(crp, cogOutcomes_colNum, cogOutcomes_fctr, brainMorphVars_Z_colNum, covars_colNum, covars_brain_colNum) # list of variables used in correlation/mediation analyses
essentialVarsAnalysis_colNames <- colnames(df[essentialVarsAnalysis_colNum]) 
# essentialVarsAnalysisNum <- c(crp, crp_log_z, brainMorphVars_specific_z, cogPCOutputs_z, revisedCovars, covars_brain)
# saveRDS(essentialVarsAnalysis, file = "revisedEssentialVarsAnalysis.rds")
eid <- which(colnames(df) == "eid")
missingEssVar <- which(colnames(df) == "missingEssentialVar")

df <- df[c(eid, essentialVarsAnalysis_colNum, missingEssVar)]
# colnames(df_small)

ordinalFactorContrasts <- list(exercise_IPAQActivityGroup_t0 = "contr.treatment")

### 2.A.iii. Perform analyses -----
#### Correlation ----
subsetDate <- getDate() # by default will take today's date
# subsetDate <- "11_05_2022" # specify date of prior analysis if desired. Form: MM_DD_YYYY
subsetDF_path <- "./Data/Processed/Subsets"
# subsetNames <- c("NoMedNoDx") # for testing

covars <- c(covars_colNames, covars_brain_colNames)
listOfCovarLevels <- c("demo_sex_t0Male","smoke_currently0_t0Yes","smoke_currently0_t0Yes","demo_ethnicity_t0White","exercise_IPAQActivityGroup_t0low","exercise_IPAQActivityGroup_t0moderate","ses_townsend_t0_z","demo_age_assess0_t0_z","demo_daysBtwAssess_z","weight_waistToHip_mean02_z","sleep_duration_mean02_z", "hand_t0LH","hand_t0Not","hand_t0RH", "brain_headScale_t2", "med_Antihypertensive_t021", "med_Statin_t021")

outputPath <- "./outputs/cor"

for(subset in subsetNames){
  cat("\n", subset, ": performing correlation analysis and saving text files... \n")
  
  X <- "crp_log_z"
  outputFileNames <- c()
  varsOfInterest <- c(X, covars, cogOutcomes_colNames,cogOutcomes_fctr_names) 
  # dfSubset <- read.csv(glue("{subsetDF_path}/{subset}_{subsetDate}.csv"))
  dfSubset <- readRDS(glue("{subsetDF_path}/{subset}_{subsetDate}.rds"))
  dfSubset <- dfSubset[varsOfInterest]
  analysisName <- glue("cor_{subset}")
  
  # for numeric variables
  for(i in cogOutcomes_colNames){
    analysisName <- glue("{analysisName}_{i}")
    individualAnalysis <- runCor(df = dfSubset, X = X, Y = i, covars = covars, M = "", YisFactor = F, name = analysisName, contrasts = ordinalFactorContrasts)
    fileName <- glue("{analysisName}_{getDate()}")
    (saveRDS(individualAnalysis, file = glue("./outputs/cor/model/{fileName}")) + cat("Correlation model for subset ", subset, "and cog var", i, "saved as ", glue("./outputs/cor/model/{fileName}")))
    capture.output(summary(individualAnalysis), file = glue("{outputPath}/{fileName}.txt"))
    outputFileNames <- c(outputFileNames, fileName)
    analysisName <- glue("cor_{subset}")
  }
  
  for(outputFile in outputFileNames){
    path <- getwd()
    corTxtFormat(filePath = path, files = outputFileNames, X = X, listOfCovars = listOfCovarLevels)
  } # summarise outputs
  
  for(i in cogOutcomes_fctr_names){
    analysisName <- glue("{analysisName}_{i}")
    individualAnalysis <- runCor(df = dfSubset, X = X, Y = i, covars = colnames(df[covars]), M = "", YisFactor = T, name = analysisName, contrasts = ordinalFactorContrasts)
    fileName <- glue("{analysisName}_{getDate()}")
    capture.output(summary(individualAnalysis), file = glue("{outputPath}/{fileName}.txt"))
    outputFileNames <- c(outputFileNames, fileName)
    analysisName <- glue("cor_c_{subset}")
  }
  
  outputFileNames <- c()
  
}

#### Format ----
list_covarAllModels <- list() # this will store all dataframes with covariate information
corResultPath <- "./outputs/cor"
dataDate <- getDate()
# dataDate <- "12_14_2022" # MM_DD_YYYY
outputPath <- "./outputs/cor/summary"
# subsetNames <- "noDxNoSSRI" # for testing

corCSVFormat(depVar = "crp_log_z", modelName = "cor", outputPrefix = "cor_summary", corResultPath = corResultPath, listOfCovars = covars, dataFilePrefix = "cor", subsets = subsetNames, dataDate = dataDate, cogVars = cogOutcomes_colNames, outputPath = outputPath)
#### Correct p-values ----
require(stats)

dataDate <- getDate()
# dataDate <- "12_14_2022" # MM_DD_YYYY

dataPath <- "./outputs/cor/summary" # path to raw mediation output summary files
fileName <- "cor_summary"
outputPath <- "./outputs/cor/summary"
outputName <- "cor_summary_pCor"

fileName <- glue("{dataPath}/{file}_dataDate.csv")
summaryDF <- read_csv(fileName, show_col_types = F)
outputDF <- matrix() # ncol = 16, dimnames = list(c(), c(colnames(summaryDF)))

for(subset in subsetNames){
  outputDF_subset <- filter_rows(summaryDF, subset == subset)
  p_cor <- p.adjust(outputDF_subset$b_p, method = "fdr")
  outputDF_subset <- outputDF_subset %>%   
    mutate("b_p_cor" = p_cor, .after = "b_p")
  outputDF <- rbind(outputDF, outputDF_subset) 
}

saveName <- glue("{outputName}_{getDate()}.csv")
(write.csv(outputDF, file = glue("{outputPath}/{saveName}")) + 
    cat("\n Correlation summary with FDR correction for each subset saved as: ", saveName))


### 2.A.iv. Correl assumption checks ----
# Correctly specified
# avPlots(PC1_allC, pt.wts = T) # if slope is 0 then the predictor is not predicting unique variance in the model. Slope of lines should be equivalent to the regression coefs in the model.
# avPlots(mainCor_C_cog2, pt.wts = T)
# avPlots(mainCor_C_cog3, pt.wts = T)
# cat("\n \t\t Remove variable that is not predicting unique variance. Do so one at a time--excluding the one with a slope closest to 0 if applicable--and do so only at the end!")
# !! only remove predictors at the end !!
###
# plot(df_complete$demo_ethnicity_t0, df_complete$cog_PC_1_z, xlab = "", las = 3, res = 30000)
# ggplot(df_complete, aes(reorder(demo_ethnicity_t0, -cog_PC_1_z, sum), cog_PC_1_z)) + 
#   geom_boxplot() + 
#   geom_jitter(position=position_jitter(.01)) +
#   theme(axis.text.x = element_text(angle = 90), axis.line.x.top = element_line(1)) +
#   geom_point(alpha = .15, size = 1, color = "black")
# 
# ggplot(data = df_complete, aes(x=demo_daysBtwAssess, y=cog_PC_1_z)) + 
#                 geom_point(alpha = .15, size = 1, color = "black") +
#                 theme(legend.position = "none")


## 2.B Mediation analyses ----
### 2.B.i. Functions ----
# where M: continuous mediator; X: independent variable; Y: continuous dependent variable; C: covariates
runRMediation <- function(df, X, Y, M, covars, contrasts, YisFactor, name, alpha, outputPath){
  require(RMediation)
  relevantVars <- c()
  for(i in c(X, Y, M, covars)){
    if(i != ""){
      relevantVars <- c(relevantVars, i)
    }
  }
  df <- df[complete.cases(df[c(relevantVars)]),]
  totModelExpression <- createLMExpression(varNames = colnames(df), X = X, Y = Y, M = "", C = covars) # total model: predict Y from X and C
  aPathExpression <- createLMExpression(varNames = colnames(df), X = X, Y = M, M = "", C = covars) # a path model: predict M from X and C
  bPathExpression <- createLMExpression(varNames = colnames(df), X = X, Y = Y, M = M, C = covars) # b path & direct path model: predict Y from X, M, and C
  
  # total effect
  if(YisFactor == TRUE){
    totModel <- glm(eval(parse(text = totModelExpression)), data = df, contrasts = contrasts, family = "binomial")
  } else{
    totModel <- lm(eval(parse(text = totModelExpression)), data = df, contrasts = contrasts)
  }
  # a path -- predict M from X
  aModel <- lm(eval(parse(text = aPathExpression)), data = df,  contrasts = contrasts)
  
  # b path -- predict Y from X, M and C
  if(YisFactor == TRUE){
    bModel <- glm(eval(parse(text = bPathExpression)), data = df, contrasts = contrasts, family = "binomial")
  } else{
    bModel <- lm(eval(parse(text = bPathExpression)), data = df,  contrasts = contrasts)
  }
  
  a <- coef(aModel)[which(names(coef(aModel)) == X)] # coef for effect of X on M
  a.se <- coef(summary(aModel))[which(names(coef(aModel)) == X),2] # S.E. of the regression coef of X on M. N.b. index is [row of X, column of standard errors of estimate]
  
  b <- coef(bModel)[which(names(coef(bModel)) == M)] # coef for effect of M on Y
  b.se <- coefficients(summary(bModel))[which(names(coef(bModel)) == M),2] # coef for effect of M on Y
  
  
  fileName <- glue("{name}_{getDate()}")
  
  med <- medci(mu.x = a, mu.y = b, se.x = a.se, se.y = b.se, rho = 0, alpha = alpha, plot = F, type = "asymp")
  capture.output(
    cat(" -- Linear Model Analysis Output -- "),
    cat("\n ------------------ \n 'total model' linear model expression: \n\t", totModelExpression),
    print(summary(totModel)),
    cat("\n ------------------ \n 'a path' linear model expression: \n\t", aPathExpression),
    print(summary(aModel)),
    cat("\n ------------------ \n 'b path' linear model expression: \n\t", bPathExpression),
    print(summary(bModel)),
    cat("\n ------------------ \n Mediated Effect \n\t "),
    cat("Estimate: \t", med$Estimate, "\n\t SE:\t\t",  med$SE, "\n\t", 100*(1-alpha), "%CILow:\t", med$`97.5% CI`[1], "\n\t", 100*(1-alpha),"%CIHi:\t", med$`97.5% CI`[2],"\n"),
    file = glue("{outputPath}/{fileName}.txt"))
  
  return(fileName)
} # input: 'df' - dataframe containing all relevant variables; 'X' - name of the predictor variable; 'Y' - name of the outcome variable; 'M' - name of the mediator variable; 'C' - list of names of covariate variables; 'contrasts' - a list of ordinal variables specifying their contrasts; 'YisFactor' - if outcome is a factor, logical. If yes, logistic regression will be performed.; 'name' - name of analysis to be used for output filename; 'alpha' - critical p-value. Used for calculating confidence interval; 'outputPath' - path to save temporary text file.

formatMedOutput <- function(outputDF, dataPath, fileName, analysisName, logReg){
  specificFilePath <- glue("{dataPath}/{fileName}.txt")
  fileObj <- readLines(specificFilePath)
  
  splits <- grep("~", fileObj)
  split1 <- splits[2]-2
  split2 <- splits[3]-2
  medOutput <- grep("Mediated Effect", fileObj)
  
  text_totEffect <- fileObj[2:split1]
  text_aPath <- fileObj[(split1+1):split2]
  text_bPath <- fileObj[(split2+1):(medOutput-1)]
  text_mediatedEffect <- fileObj[medOutput:length(fileObj)]
  models <- list(text_totEffect, text_aPath, text_bPath)
  
  variables <- c()
  values <- c()
  
  for(model in 1:(length(models))){ # iterate through all three models
    sectionText <- unlist(models[model])
    lmExpression <- grep("~", sectionText, value = TRUE) # from this line, extract analysis name
    lmExpression <- unlist(strsplit(lmExpression, "\\s+"))
    Y <- lmExpression[which(lmExpression == "~")-1] # Y is not necessarily the same Y of mediation analyses
    if(model != 3){
      X <- lmExpression[which(lmExpression == "~")+1]
    } else {
      X <- lmExpression[which(lmExpression == "~")+3]
    }
    
    X_values  <- grep(X, sectionText, value = TRUE)
    if(length(X_values) == 2){
      X_values <- X_values[2]
    }
    
    X_values <- unlist(strsplit(X_values, "\\s+"))
    
    if(TRUE %in% grepl("<", X_values) && nchar(X_values[grep("<",X_values)]) == 1){
      pAndSymbol <- glue("{X_values[[5]]}{X_values[[6]]}")
      X_values_NoSigChar <- c(X_values[2:4], pAndSymbol)
    } else {
      X_values_NoSigChar <- X_values[2:5]
    }
    
    CI_lo <- as.numeric(X_values_NoSigChar[1]) - qnorm(.975)*as.numeric(X_values_NoSigChar[2]) # calculate 95%CI of parameter estimate (X_values_NoSigChar[1]) using its standard error (X_values_NoSigChar[2])
    CI_hi <- as.numeric(X_values_NoSigChar[1]) + qnorm(.975)*as.numeric(X_values_NoSigChar[2])
  
    if(logReg == F || model == 2){ # the a path will always be linear regression as all brain morphology variables are continuous
      effectSizes <- grep("Multiple R-squared",sectionText, value = TRUE)
      effectSizes <- unlist(strsplit(effectSizes, "\\s+"))
      R2 <- substr(effectSizes[3], 1, nchar(effectSizes[3])-1)
      R2adj <- effectSizes[6]
      
      FTest <- grep("F-statistic",sectionText, value = TRUE)
      FTest <- unlist(strsplit(FTest, "\\s+"))
      FStat <- FTest[2]
      df1 <- FTest[4]
      df2 <- FTest[6]
      
      p <- FTest[which(FTest == grep("p-value", FTest, value = TRUE)):length(FTest)]
      if(TRUE %in% grepl("<", p) && nchar(p[grep("<",p)]) == 1){
        p <- glue("{FTest[[9]]}{FTest[[10]]}")
        # p <- paste(FTest[[9]], FTest[[10]], sep = "")
      } else {
        p <- FTest[which(FTest == grep("p-value", FTest, value = TRUE))+1]
      }
      
      variables <- append(variables, c(X, Y))
      
      if(model == 2){
        values <- append(values, c(X_values_NoSigChar[1:2], CI_lo, CI_hi,X_values_NoSigChar[3:4], R2, R2adj, FStat, df1, df2, p))
      } else {
        AIC <- "NA"
        values <- append(values, c(X_values_NoSigChar[1:2], CI_lo, CI_hi,X_values_NoSigChar[3:4], R2, R2adj, FStat, df1, df2, p, AIC))
      }
      
    } else{
      R2 <- "NA"
      R2adj <- "NA"
      FStat <- "NA"
      
      df <- grep("Null deviance:",sectionText, value = TRUE)
      df <- unlist(strsplit(df, "\\s+"))
      # print(df)
      df1  <- df[4] # null deviance df
      df2 <- df[6]
      
      p <- "NA"
      AIC <- grep("AIC:",sectionText, value = TRUE)
      AIC <- unlist(strsplit(AIC , "\\s+"))
      # print(AIC)
      AIC  <-  AIC[2]
      
      variables <- append(variables, c(X, Y))
      values <- append(values, c(X_values_NoSigChar[1:2], CI_lo, CI_hi,X_values_NoSigChar[3:4], R2, R2adj, FStat, df1, df2, p, AIC))
    }
  }
  
  X <- variables[3]
  Y <- variables[2]
  M <- variables[4]
  vars <- c(X, M, Y)
  name <- glue("{analysisName}_{M}_{Y}")
  
  # Extract indirect effect
  sectionText <- text_mediatedEffect
  
  medEff_est <- grep("Estimate",sectionText, value = TRUE)
  medEff_est <- unlist(strsplit(medEff_est, "\\s+"))
  medEff_est <-  medEff_est[3]
  
  medEff_SE <- grep("SE",sectionText, value = TRUE)
  medEff_SE <- unlist(strsplit(medEff_SE, "\\s+"))
  medEff_SE <-  medEff_SE[3]
  
  medEff_95Low <- grep("%CILow",sectionText, value = TRUE)
  medEff_95Low <- unlist(strsplit(medEff_95Low, "\\s+"))
  medEff_95Low <-  medEff_95Low[4]
  
  medEff_95Hi <- grep("%CIHi",sectionText, value = TRUE)
  medEff_95Hi <- unlist(strsplit(medEff_95Hi, "\\s+"))
  
  medEff_95Hi <- medEff_95Hi[4]
  
  indEff <- c(medEff_est, medEff_SE, medEff_95Low, medEff_95Hi)
  
  outputDF <- rbind(outputDF, c(name, vars, indEff, values)) # specify the indices of valueList to retain in main_df
  
  return(outputDF)
} # takes: outputDF: formatted matrix containing summary of each mediation analysis; filePath: location of text file containing mediation analysis output; fileName: name of text file containint mediation analysis; analysisName: name of analysis, e.g., 'med_noC'; logReg: logical indicates if the model is a logistic regression (i.e., discrete outcome variable)

getPMed <- function(dataPath, subsetNames, dataFileInfix, dataDate, keyCognition, keyMediators, pCor, sigOnly, alpha, round, outputPath) {
  outputDf_p_cor <- c()
  
  for (subset in subsetNames) {
    cat(subset, " get and correct p-values \n")
    file <- glue("{subset}_{dataFileInfix}_{dataDate}")
    
    fileName <- glue("{dataPath}/{file}.csv")
    outputDF <- read_csv(fileName, show_col_types = F)
    
    cat("\t File: ", file, "\n")
    
    # Sobel's test for computing mediation effect's p-value
    z <- outputDF$IndEff_est / outputDF$IndEff_SE 
    p <- pnorm(z) # compute p-value
    outputDF <- outputDF %>%
      mutate("IndEff_p" = p, .after = "IndEff_SE")
  
    if (pCor == TRUE) {
      fileNameAffix <- "pCor"
      
      for (cog in keyCognition) {
        outputDF_keyYM <- outputDF %>%
          filter(Y %in% cog) %>%
          filter(M %in% keyMediators)
        
        indEff_p_cor <- p.adjust(outputDF_keyYM$IndEff_p, method = "fdr")
        totEff_p_cor <- p.adjust(outputDF_keyYM$totEff_X_p, method = "fdr")
        a_p_cor <- p.adjust(outputDF_keyYM$aPath_X_p, method = "fdr")
        b_p_cor <- p.adjust(outputDF_keyYM$bPath_X_p, method = "fdr")
        
        outputDF_keyYM <- outputDF_keyYM %>%
          mutate("IndEff_p_cor" = indEff_p_cor, .after = "IndEff_p") %>% 
          mutate("totEff_X_p_cor" = totEff_p_cor, .after = "totEff_X_p") %>%
          mutate("aPath_X_p_cor" = a_p_cor, .after = "aPath_X_p") %>%
          mutate("bPath_X_p_cor" = b_p_cor, .after = "bPath_X_p")
        
        outputDF_keyYM <- outputDF_keyYM[order(outputDF_keyYM$IndEff_p_cor), ] # order by corrected p-value
        
        outputDF_keyYM_sig <- outputDF_keyYM %>% filter(IndEff_p_cor <= alpha)
        cat("\n\t Cog variable: ", cog, ". \n\t\t Significant mediators:")
        if (length(outputDF_keyYM_sig$M) == 0) {
          smallest_NonSigP <- round(min(outputDF_keyYM$IndEff_p_cor), 3)
          cat("\t None. All corrected p >= ", smallest_NonSigP)
        } else {
          cat("(num. = ", length(outputDF_keyYM_sig$M), ") ")
          for (sigMediator in outputDF_keyYM_sig$M) {
            cat("\n\t\t", sigMediator)
          }
          outputDF_keyYM_nonSig <- outputDF_keyYM %>% filter(IndEff_p_cor > alpha)
          smallest_NonSigP <- round(min(outputDF_keyYM_nonSig$IndEff_p_cor), 3)
          cat("\n Remaining indirect effects, corrected p >= ", smallest_NonSigP)
        }
        
        if (sigOnly == TRUE) {
          outputDF_keyYM <- outputDF_keyYM_sig
        }
        
        outputDf_p_cor <- rbind(outputDf_p_cor, outputDF_keyYM)
      }
      
    } else {
      fileNameAffix <- "p"
    }
    
    outputDf_p_cor <- as.data.frame(sapply(outputDf_p_cor, function(.)format(., scientific = FALSE)))
    
    outputDf_p_cor <- outputDf_p_cor %>% 
      mutate("subset" = subset, .after = "modelName") %>% 
      mutate("round" = round, .after = "subset")
    
    saveName <- glue("{subset}_{dataFileInfix}_{fileNameAffix}_round{round}_{getDate()}.csv")
    write.csv(outputDf_p_cor, file = glue("{outputPath}/{saveName}"))
    
    cat("Saved p-value file: ", glue("{outputPath}/{saveName}"), "\n---------------------------")
    
    outputDf_p_cor <- c()
  }
} # 'dataPath': path to mediation output summary files; 'subsetNames': list of subset names of interest; 'dataFileInfix': infix in summary mediation file name; 'dataDate': date of summary mediation files as shown in file name; 'pCor': logical specifying if p-values are to be corrected. If TRUE, then will correct considering all mediators listed in 'keyMediators'; 'sigOnly': logical specifying if output file should return only significant mediation models (uses alpha given by 'alpha'); 'round': number indicating the round of the mediation analysis; 'outputPath': string of path defining desired output directory.

### 2.B.ii. Define variables of interest -----
# df <- readRDS("Data/Processed/Subsets/All_12_14_2022.rds")

brainPrefix <- c("area_","mThick_", "vol_", "weigted_mFA_")
brainMorphVars <- c()
for(i in brainPrefix){
  iVars <- starts_with(i, vars = colnames(df))
  # cat(paste(i, ":", length(iVars), "\n"), sep = "")
  brainMorphVars <- c(brainMorphVars, iVars)
}
brainMorphVars_colNames <- colnames(df[brainMorphVars])
df_brainMorphVarsZ <- df[brainMorphVars] %>% 
  mutate(across(.fns = scale, .names = "{colnames(df[brainMorphVars])}_z")) # standardize brain morphology variables
df_brainMorphVarsZ <- df_brainMorphVarsZ[,which(! colnames(df_brainMorphVarsZ) %in% colnames(df[brainMorphVars]))] # matrix of z-scores made in above line to add to main df
brainMorphVars_Z_colNames <- colnames(df_brainMorphVarsZ)
df <- cbind(df, df_brainMorphVarsZ[brainMorphVars_Z_colNames])
rm(df_brainMorphVarsZ)
brainMorphVars_Z_colNum <- which(colnames(df) %in% brainMorphVars_Z_colNames)

generalMetrics <- c("vol_BrainSeg_WB_t2","vol_BrainSegNotVent_WB_t2","vol_BrainSegNotVentSurf_WB_t2","vol_SubCortGray_WB_t2","vol_TotalGray_WB_t2","vol_SupraTentorial_WB_t2","vol_SupraTentorialNotVent_WB_t2","vol_EstimatedTotalIntraCranial_WB_t2","vol_BrainStem_WB_t2","vol_CSF_WB_t2","vol_wmHyperintensities_WB_t2","vol_WMhypointensities_WB_t2","vol_Hippocampus_L_t2","vol_Hippocampus_R_t2","vol_CerebellumCortex_R_t2","vol_CerebellumCortex_L_t2","vol_CerebellumWhiteMatter_L_t2","vol_CerebellumWhiteMatter_R_t2")
aseg_brainVars <- c("vol_VentricleChoroid_WB_t2","vol_3rdVentricle_WB_t2","vol_4thVentricle_WB_t2","vol_5thVentricle_WB_t2","vol_nonWMhypointensities_WB_t2","vol_OpticChiasm_WB_t2","vol_CCPosterior_WB_t2","vol_CCMidPosterior_WB_t2","vol_CCCentral_WB_t2","vol_CCMidAnterior_WB_t2","vol_CCAnterior_WB_t2","area_Accumbensarea_L_t2","vol_Cortex_L_t2","vol_CerebralWhiteMatter_L_t2","vol_LateralVentricle_L_t2","vol_InfLatVent_L_t2","vol_thalamusProper_L_t2","vol_Caudate_L_t2","vol_Putamen_L_t2","vol_Pallidum_L_t2","vol_Amygdala_L_t2","vol_Accumbensarea_L_t2","vol_VentralDC_L_t2","vol_choroidplexus_L_t2","area_Accumbensarea_R_t2","vol_Cortex_R_t2","vol_CerebralWhiteMatter_R_t2","vol_LateralVentricle_R_t2","vol_InfLatVent_R_t2","vol_thalamusProper_R_t2","vol_Caudate_R_t2","vol_Putamen_R_t2","vol_Pallidum_R_t2","vol_Amygdala_R_t2","vol_Accumbensarea_R_t2","vol_VentralDC_R_t2","vol_vessel_L_t2","vol_vessel_R_t2","vol_choroidplexus_R_t2") # all Freesurfer ASEG variables in UKBB
dkt_brainVars_area <- c("area_Caudalanteriorcingulate_L_t2","area_Caudalmiddlefrontal_L_t2","area_Cuneus_L_t2","area_Entorhinal_L_t2","area_Fusiform_L_t2","area_Inferiorparietal_L_t2","area_Inferiortemporal_L_t2","area_Isthmuscingulate_L_t2","area_Lateraloccipital_L_t2","area_Lateralorbitofrontal_L_t2","area_Lingual_L_t2","area_Medialorbitofrontal_L_t2","area_Middletemporal_L_t2","area_Parahippocampal_L_t2","area_Paracentral_L_t2","area_Parsopercularis_L_t2","area_Parsorbitalis_L_t2","area_Parstriangularis_L_t2","area_Pericalcarine_L_t2","area_Postcentral_L_t2","area_Posteriorcingulate_L_t2","area_Precentral_L_t2","area_Precuneus_L_t2","area_Rostralanteriorcingulate_L_t2","area_Rostralmiddlefrontal_L_t2","area_Superiorfrontal_L_t2","area_Superiorparietal_L_t2","area_Superiortemporal_L_t2","area_Supramarginal_L_t2","area_Transversetemporal_L_t2","area_Insula_L_t2", "area_Caudalanteriorcingulate_R_t2","area_Caudalmiddlefrontal_R_t2","area_Cuneus_R_t2","area_Entorhinal_R_t2","area_Fusiform_R_t2","area_Inferiorparietal_R_t2","area_Inferiortemporal_R_t2","area_Isthmuscingulate_R_t2","area_Lateraloccipital_R_t2","area_Lateralorbitofrontal_R_t2","area_Lingual_R_t2","area_Medialorbitofrontal_R_t2","area_Middletemporal_R_t2","area_Parahippocampal_R_t2","area_Paracentral_R_t2","area_Parsopercularis_R_t2","area_Parsorbitalis_R_t2","area_Parstriangularis_R_t2","area_Pericalcarine_R_t2","area_Postcentral_R_t2","area_Posteriorcingulate_R_t2","area_Precentral_R_t2","area_Precuneus_R_t2","area_Rostralanteriorcingulate_R_t2","area_Rostralmiddlefrontal_R_t2","area_Superiorfrontal_R_t2","area_Superiorparietal_R_t2","area_Superiortemporal_R_t2","area_Supramarginal_R_t2","area_Transversetemporal_R_t2","area_Insula_R_t2")
dkt_brainVars_vol <- c("vol_Caudalanteriorcingulate_R_t2","vol_Caudalmiddlefrontal_R_t2","vol_Cuneus_R_t2","vol_Entorhinal_R_t2","vol_Fusiform_R_t2","vol_Inferiorparietal_R_t2","vol_Inferiortemporal_R_t2","vol_Isthmuscingulate_R_t2","vol_Lateraloccipital_R_t2","vol_Lateralorbitofrontal_R_t2","vol_Lingual_R_t2","vol_Medialorbitofrontal_R_t2","vol_Middletemporal_R_t2","vol_Parahippocampal_R_t2","vol_Paracentral_R_t2","vol_Parsopercularis_R_t2","vol_Parsorbitalis_R_t2","vol_Parstriangularis_R_t2","vol_Pericalcarine_R_t2","vol_Postcentral_R_t2","vol_Posteriorcingulate_R_t2","vol_Precentral_R_t2","vol_Precuneus_R_t2","vol_Rostralanteriorcingulate_R_t2","vol_Rostralmiddlefrontal_R_t2","vol_Superiorfrontal_R_t2","vol_Superiorparietal_R_t2","vol_Superiortemporal_R_t2","vol_Supramarginal_R_t2","vol_Transversetemporal_R_t2","vol_Insula_R_t2", "vol_Caudalanteriorcingulate_L_t2","vol_Caudalmiddlefrontal_L_t2","vol_Cuneus_L_t2","vol_Entorhinal_L_t2","vol_Fusiform_L_t2","vol_Inferiorparietal_L_t2","vol_Inferiortemporal_L_t2","vol_Isthmuscingulate_L_t2","vol_Lateraloccipital_L_t2","vol_Lateralorbitofrontal_L_t2","vol_Lingual_L_t2","vol_Medialorbitofrontal_L_t2","vol_Middletemporal_L_t2","vol_Parahippocampal_L_t2","vol_Paracentral_L_t2","vol_Parsopercularis_L_t2","vol_Parsorbitalis_L_t2","vol_Parstriangularis_L_t2","vol_Pericalcarine_L_t2","vol_Postcentral_L_t2","vol_Posteriorcingulate_L_t2","vol_Precentral_L_t2","vol_Precuneus_L_t2","vol_Rostralanteriorcingulate_L_t2","vol_Rostralmiddlefrontal_L_t2","vol_Superiorfrontal_L_t2","vol_Superiorparietal_L_t2","vol_Superiortemporal_L_t2","vol_Supramarginal_L_t2","vol_Transversetemporal_L_t2","vol_Insula_L_t2")
dkt_brainVars_mThick <- c("mThick_Caudalanteriorcingulate_L_t2","mThick_Caudalmiddlefrontal_L_t2","mThick_Cuneus_L_t2","mThick_Entorhinal_L_t2","mThick_Fusiform_L_t2","mThick_Inferiorparietal_L_t2","mThick_Inferiortemporal_L_t2","mThick_Isthmuscingulate_L_t2","mThick_Lateraloccipital_L_t2","mThick_Lateralorbitofrontal_L_t2","mThick_Lingual_L_t2","mThick_Medialorbitofrontal_L_t2","mThick_Middletemporal_L_t2","mThick_Parahippocampal_L_t2","mThick_Paracentral_L_t2","mThick_Parsopercularis_L_t2","mThick_Parsorbitalis_L_t2","mThick_Parstriangularis_L_t2","mThick_Pericalcarine_L_t2","mThick_Postcentral_L_t2","mThick_Posteriorcingulate_L_t2","mThick_Precentral_L_t2","mThick_Precuneus_L_t2","mThick_Rostralanteriorcingulate_L_t2","mThick_Rostralmiddlefrontal_L_t2","mThick_Superiorfrontal_L_t2","mThick_Superiorparietal_L_t2","mThick_Superiortemporal_L_t2","mThick_Supramarginal_L_t2","mThick_Transversetemporal_L_t2","mThick_Insula_L_t2","mThick_Caudalanteriorcingulate_R_t2","mThick_Caudalmiddlefrontal_R_t2","mThick_Cuneus_R_t2","mThick_Entorhinal_R_t2","mThick_Fusiform_R_t2","mThick_Inferiorparietal_R_t2","mThick_Inferiortemporal_R_t2","mThick_Isthmuscingulate_R_t2","mThick_Lateraloccipital_R_t2","mThick_Lateralorbitofrontal_R_t2","mThick_Lingual_R_t2","mThick_Medialorbitofrontal_R_t2","mThick_Middletemporal_R_t2","mThick_Parahippocampal_R_t2","mThick_Paracentral_R_t2","mThick_Parsopercularis_R_t2","mThick_Parsorbitalis_R_t2","mThick_Parstriangularis_R_t2","mThick_Pericalcarine_R_t2","mThick_Postcentral_R_t2","mThick_Posteriorcingulate_R_t2","mThick_Precentral_R_t2","mThick_Precuneus_R_t2","mThick_Rostralanteriorcingulate_R_t2","mThick_Rostralmiddlefrontal_R_t2","mThick_Superiorfrontal_R_t2","mThick_Superiorparietal_R_t2","mThick_Superiortemporal_R_t2","mThick_Supramarginal_R_t2","mThick_Transversetemporal_R_t2","mThick_Insula_R_t2") 

lobarBrainMetrics <- c("area_insula_WB_t2", "vol_insula_WB_t2", "mThick_insula_WB_t2")
for(metric in c("area", "mThick", "vol")){
  for(lobe in c("frontal", "parietal", "occipital", "temporal")){
    for(hemi in c("L", "R", "WB")){ 
      lobarBrainMetrics <- c(lobarBrainMetrics, glue("{metric}_{lobe}_{hemi}_t2"))
      # print(glue("{metric}_{lobe}_{hemi}"))
    }
  }
}

round1Mediators <- c(lobarBrainMetrics, generalMetrics)
round2Mediators <- c(aseg_brainVars, dkt_brainVars_area, dkt_brainVars_vol, dkt_brainVars_mThick)

mediators <- c(round1Mediators, round2Mediators)
# mediators <- c("area_Caudalanteriorcingulate_L_t2") # for testing purposes

### 2.B.iii. Perform analyses ----
cat("WARNING. You are about to compute mediation analyses. \n\t As one analysis is computed for each combination of subset, cognitive variable and brain variable, these analyses may take many hours to complete.")

#### Mediation ----
subsetDate <- getDate() # by default will take today's date
subsetDate <- "12_14_2022" # specify date of prior analysis if desired. Form: MM_DD_YYYY
subsetDF_path <- "./Data/Processed/Subsets"
# subsetNames <- c("NoMedNoDx") # for testing
tempOutputPath <- "./outputs/med"
outputPath <- "./outputs/med/subsets"

covarsMed <- c(covars_colNames, covars_brain_colNames) # Add  cognition only covars
# colnames(df[covarsMed])

# brainMorphVars_colNames <- mediators
# brainMorphVars_colNames <- c("brain_vol_brainSegNoVent_t2_z","brain_vol_hippocamp_L_t2_z","brain_vol_hippocamp_R_t2_z") # list of mediator variable names

# cogOutcomes_colNames <- c("cog_numMem_maxDigitRemem_t2_z", "cog_fluidIntel_score_t2_z")

individualAnalysis <- list()
listAnalyses_med_C <- list()
mainMed_C <- list()
ordinalFactorContrasts <- c()

for(subset in subsetNames){
  df <- readRDS(glue("{subsetDF_path}/{subset}_{subsetDate}.rds"))
  df_brainMorphVarsZ <- df[brainMorphVars_colNames] %>% 
    mutate(across(.fns = scale, .names = "{colnames(df[brainMorphVars])}_z")) # standardize brain morphology variables
  df_brainMorphVarsZ <- df_brainMorphVarsZ[,which(! colnames(df_brainMorphVarsZ) %in% colnames(df[brainMorphVars]))] # matrix of z-scores made in above line to add to main df
  brainMorphVars_Z_colNames <- colnames(df_brainMorphVarsZ)
  df <- cbind(df, df_brainMorphVarsZ[brainMorphVars_Z_colNames])
  mediators <- brainMorphVars_Z_colNames
  
  outputDF <- matrix(ncol = 46, dimnames = list(c(), c("modelName", "X", "M", "Y", "IndEff_est", "IndEff_SE", "IndEff_95CI_Lo","IndEff_95CI_Hi","totEff_b", "totEff_SE", "totEff_95CI_Lo","totEff_95CI_Hi", "totEff_X_t", "totEff_X_p", "totEff_R2", "totEff_R2adj", "totEff_F", "totEff_df1", "totEff_df2", "totEff_model_p", "totEff_AIC", "aPath_b", "aPath_SE", "aPath_95CI_Lo", "aPath_95CI_Hi", "aPath_X_t", "aPath_X_p", "aPath_model_R2", "aPath_model_R2adj", "aPath_model_F", "aPath_model_df1", "aPath_model_df2", "aPath_model_p", "bPath_b", "bPath_SE", "bPath_95CI_Lo", "bPath_95CI_Hi","bPath_X_t", "bPath_X_p", "bPath_model_R2", "bPath_model_R2adj", "bPath_model_F", "bPath_model_df1", "bPath_model_df2", "bPath_model_p", "bPath_AIC")))
  
  analysisName <- glue("{subset}_med")
   cat("\n Mediation analyses for ", subset)
  for(med in colnames(df[mediators])){
    cat("\t mediator: ", med, "\n")
    for(cog in cogOutcomes_colNames){
      analysisName <- glue("{analysisName}_{med}_{cog}")
      fileName <- runRMediation(df = df, X = "crp_log_z", 
                                Y = cog, M = med, 
                                covars = colnames(df[covarsMed]), 
                                contrasts = ordinalFactorContrasts, 
                                YisFactor = F, 
                                name = analysisName, alpha = .05,
                                outputPath = tempOutputPath)
      
      analysisName <- glue("{subset}_med")
      outputDF <- formatMedOutput(outputDF = outputDF, 
                                  fileName = fileName, 
                                  dataPath = tempOutputPath, 
                                  analysisName = analysisName,
                                  logReg = F) # reads text file and adds its data to a summary dataframe
      
      fileToRemove <- as.character(glue("{tempOutputPath}/{fileName}.txt")) #
      file.remove(fileToRemove) # remove raw text file
    }
    
    for(cog_fctr in cogOutcomes_fctr_names){
      analysisName <- glue("{analysisName}_{med}_{cog_fctr}")
      fileName <- runRMediation(df = df, X = "crp_log_z", 
                                Y = cog_fctr, M = med, 
                                covars = colnames(df[covarsMed]), 
                                contrasts = ordinalFactorContrasts, 
                                YisFactor = T, 
                                name = analysisName, alpha = .05,
                                outputPath = tempOutputPath)
      
      analysisName <- glue("{subset}_med_C")
      outputDF <- formatMedOutput(outputDF = outputDF, 
                                  fileName = fileName, 
                                  dataPath = tempOutputPath, 
                                  analysisName = analysisName,  
                                  logReg = T)
      
      fileToRemove <- as.character(glue("{tempOutputPath}/{fileName}.txt")) 
      file.remove(fileToRemove) # remove raw text file
    }
  }
  outputDF <- as.data.frame(outputDF)
  outputDF <- outputDF[-1,] # remove first row if as it is all NA values
  (write.csv(outputDF, file = glue("{outputPath}/{analysisName}_{getDate()}.csv")) + cat("\`", Subset, "\`summary file saved: ", glue("{analysisName}_{getDate()}.csv")))
}

#### Correct p-values -----
dataPath <- "./outputs/med/subsets" # Path for file containing analysis summary csv files. Update as necessary
infix <- "med_C"
dataDate <- getDate() # by default will take today's date
# dataDate <- "12_15_2022" # specify date of prior analysis if desired. Form: MM_DD_YYYY
outputPath <- "./outputs/med/pCor"

keyCognition <- c("cog_fluidIntel_score_t2_z", "cog_numMem_maxDigitRemem_t2_z") # name of cognition variables of interest for mediation analyses

# Correct p-values for general brain metrics
keyMediators_1 <- c(lobarBrainMetrics, generalMetrics)
keyMediators_1_z <- unlist(lapply(keyMediators_1, function(.)paste(., "_z", sep = "")))
getPMed(dataPath = dataPath, 
        subsetNames = subsetNames,
        dataFileInfix = infix,
        dataDate = dataDate,
        outputPath = glue("{outputPath}_round1"),
        keyMediators = keyMediators_1_z, 
        keyCognition = keyCognition, 
        pCor = T, sigOnly = F, alpha = .05, 
        round = 1) # for round 1, general brain metrics

# Correct p-values for specific brain metrics
keyMediators_2 <- c(aseg_brainVars, dkt_brainVars_area, dkt_brainVars_vol, dkt_brainVars_mThick) # list of mediators to correct p-values for round 2
keyMediators_2_z <- unlist(lapply(keyMediators_2, function(.)paste(., "_z", sep = "")))
getPMed(dataPath = dataPath, 
        subsetNames = subsetNames,
        dataFileInfix = infix,
        dataDate = dataDate,
        outputPath = glue("{outputPath}_round2"),
        keyMediators = keyMediators_2_z, 
        keyCognition = keyCognition, 
        pCor = T, sigOnly = F, alpha = .05, 
        round = 2) # for round 2, specific brain regions

# combine round 1 and 2 outputs across subsets 
roundsToCombine <- c(1,2) # rounds to combine
subsetsToCombine <- subsetNames

dataPath <- "./outputs/med" # Path for file containing analysis summary csv files. Update as necessary
infix <- "med_C_pCor_round"
dataDate <- getDate() # by default will take today's date
# dataDate <- "12_14_2022" # specify date of prior analysis if desired. Form: MM_DD_YYYY
outputPath <- "./outputs/med"
df_output <- c()

for(round in roundsToCombine){
  for(subset in subsetsToCombine){
    df_specific <- read.csv(file = glue("{dataPath}/pCor_round{round}/{subset}_{infix}{round}_{dataDate}.csv"))
    df_output <- unique(rbind(df_output, df_specific)) # remove duplicate rows, namely repeated column names
  }
}
(write.csv(x = df_output, file = glue("{outputPath}/MASTER_med_C_pCor_{dataDate}.csv")) + 
  cat("| Mediation files with corrected p-values for \n|\t rounds: ", roundsToCombine, "\n|\t subsets: ", subsetsToCombine, "\n| saved to ", glue("{outputPath}/MASTER_{infix}_{dataDate}.csv")))

# 4. Visualisation ----
## 4.A Figure 1: Mediation diagrams ----
### 4.A.i. Functions ----
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

round1_MedDiagrams <- function(dataPath, dataFileName, dataFileInfix, dataDate, subsets, rounds, rois, hemispheres, metrics, pCor, outputPath){
  require(glue)
  require(dplyr)
  require(stringr)
  require(readr)
  
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
  
  df <- read_csv(glue("{dataPath}/{dataFileName}{dataFileInfix}_{dataDate}.csv"), show_col_types = F)
  # cat("df read: ", glue("{dataPath}/{dataFileName}{dataFileInfix}_{dataDate}.csv")) # for debugging
  colNum_Y <- which(colnames(df) == "Y")
  colNum_M <- which(colnames(df) == "M")
  colNum_aPath_b <- which(colnames(df) == "aPath_b")
  colNum_bPath_b <- which(colnames(df) == "bPath_b")
  colNum_totEff_b <- which(colnames(df) == "totEff_b")
  colNum_indEff_b <- which(colnames(df) == "IndEff_est")
  
  if(pCor == T){
    colNum_indEff_p <- which(colnames(df) == "IndEff_p_cor")
    colNum_totEff_p <- which(colnames(df) == "totEff_X_p_cor")
    colNum_aPath_p <- which(colnames(df) == "aPath_X_p_cor")
    colNum_bPath_p <- which(colnames(df) == "bPath_X_p_cor")
  } else if(pCor == F){
    colNum_indEff_p <- which(colnames(df) == "IndEff_p")
    colNum_aPath_p <- which(colnames(df) == "aPath_X_p")
    colNum_totEff_p <- which(colnames(df) == "totEff_X_p")
    colNum_bPath_p <- which(colnames(df) == "bPath_X_p")
  } else {
    return(car("Invalid value for 'pCor'. Possible values are: \`T\` (corrected p is desired) or \`F\` (non-corrected p is desired)"))
  }
  

  for(subset_i in subsets){
    cat("Making mediation figures for: ", subset_i, "\n")
    for(round_j in rounds){
      for(roi in rois){ # For each roi, hemisphere, metric combination ----
        for(hemi in hemispheres){
          for(metric in metrics){
            medName <- glue("{metric}_{roi}_{hemi}") # for output file name
            cat(medName, "\n")
            dfMedName <- paste(glue("{medName}_t2_z")) # for finding the unique mediator in df
            # print(dfMedName)
  
            if(dfMedName %in% df$M) { # Check that this combination exists. If it doesn't skip
              df_relevant <- df %>% filter(subset == subset_i, round == round_j, M == dfMedName)
              # View(df_relevant)
              # cat("Filtered. Length: ", length(df_relevant)) # for debugging

              for(row in 1:nrow(df_relevant)) { # iterate through rows containing this mediator
                ## Determine labels for mediation diagram. ----
                ### Extract values for diagram ----
                Y <- unlist(df_relevant[row, colNum_Y]) # Extract value in column Y
                
                # format of parameter estimates: 0.###
                # cat("p-values: \n\tindEff: ", unlist(df_relevant[row, colNum_indEff_p]), "\n\t totEff: ", unlist(df_relevant[row, colNum_totEff_p]), "\n\t aPath: ", unlist(df_relevant[row, colNum_aPath_p]), "\n\t bPath: ", unlist(df_relevant[row, colNum_bPath_p]))
                
                a_eff <- format(round(unlist(df_relevant[row, colNum_aPath_b]), 3), nsmall = 3) # 'aPath_b'
                a_p <- pToStar(df_relevant[row, colNum_aPath_p]) # 'aPath_model_p'
                
                b_eff <- format(round(unlist(df_relevant[row, colNum_bPath_b]), 3), nsmall = 3) # 'bPath_b'
                b_p <- pToStar(df_relevant[row, colNum_bPath_p]) # 'bPath_model_p'
                
                c_eff <- format(
                  round(
                    unlist(
                      df_relevant[row, colNum_totEff_b]
                    ), 3
                  ), nsmall = 3
                ) # 'totEff_b'
                c_p <- pToStar(df_relevant[row, colNum_totEff_p]) # 'totEff_p'
                
                c_prime_eff <- format(round(unlist(df_relevant[row, colNum_totEff_b]) - unlist(df_relevant[row, colNum_indEff_b]), nsmall = 3) # direct effect computed by taking difference between totEff and indEff.
                
                ab_indEff <- format(round(unlist(df_relevant[row, colNum_indEff_b]), 3), nsmall = 3) # colName: "IndEff_est_correct"
                ab_p <- pToStar(df_relevant[row, colNum_indEff_p]) # colname: "IndEff_p_cor" or "IndEff_p" (defined above) 
                
                #### Make labels. ----
                ##### Line labels -----
                a <- glue("a = {a_eff}{a_p}") # Text for A path in form: "{ind eff estimate (0.###)} {significance ('*', '**', '***' or '')}"
                b <- glue("b = {b_eff}{b_p}") # Text for B path in form: "{ind eff estimate (0.###)} {significance ('*', '**', '***' or '')}"
                c <- glue("c = {c_eff}{c_p} \n ab = {ab_indEff}{ab_p} \n c\' = {c_prime_eff}") # Text for C, indirect effect (ab), and c prime in form: "c = {ind eff estimate (0.###)} {significance ('*', '**', '***' or '')} \n ab = {ind eff estimate (0.###)} {significance ('*', '**', '***' or '')} \n c\' = {dir eff estimate (0.###)}"
                
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
                } else if(roi == "BrainSegNotVent"){
                  medLabel = glue("Total brain {metricLabel}")
                } else if(roi == "insula"){
                  medLabel = glue("Insula {metricLabel} ({hemi})")
                } else if(roi == "CerebellumCortex"){
                  medLabel = glue("Cerebellum Crtx. {metricLabel} ({hemi})")
                } else if(roi == "wmHyperintensities"){
                  medLabel = glue("WM Hyperintensities {metricLabel} ({hemi})")
                } else {
                  medLabel = medName
                }
                
                ### Cognitive-variable box label. ----
                if(Y == "cog_fluidIntel_score_t2_z"){
                  lab_y <- "Fluid Intelligence Score"
                  cog <- "FluidIntelScore"
                } else if(Y == "cog_numMem_maxDigitRemem_t2_z"){
                  lab_y <- "Num. Mem. (Max dig.)"
                  cog <- "NumMemMaxDig"
                } else {
                  lab_y <- Y
                  cog <- Y
                }
                
                ## Make and save diagram. ----
                outputName <- glue("{outputPath}/medDia_{subset_i}_{cog}_{medName}_{getDate()}.png") # name of diagram when saved
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
            cat("Warning. \`", medName, "\` is an invalid combination; it is not a mediator in at least one of the provided files. \n This combination will be skipped without making a mediation diagram. \n")
          }
        }
        }
      }
    }    
  }
} # Automatically make and save diagrams for specifies brain regions and metrics. 'dataPath' string of path to folder containing output of mediation analyses !!with corrected p-values and ideally the combined MASTER document; 'dataFileName' string of name of df for specific data subset containing mediation analyses with corrected p-value for indirect effect. Df must contain the column names: 'M', 'Y', 'IndEff_est', ('IndEff_p_cor' or 'IndEff_p'), 'aPath_b', 'aPath_X_p' 'bPath_b', 'bPath_X_p'; 'dataFileInfix': infix of data file name, include any seperators between the file name and the infix; 'dataDate': date of data file; 'rois' list of brain regions to make mediation models for; 'hemispheres' hemispheres of interest. Can be 'L', 'R' or 'WB'; 'metrics' list of metrics of interest. Can be 'area', 'mThick', and/or 'vol'; 'pCor' boolean specifying if corrected p-value of ind effect is to be used; 'outputPath' string specifying where to save output diagrams

### 4.A.ii. Make diagrams -----
dataPath <- "./outputs/med" # Path for file containing analysis summary csv files. Update as necessary
dataFileName <- "MASTER"
infix <- "_med_C_pCor" # incllude any seperators
dataDate <- getDate() # by default will take today's date
# dataDate <- "12_15_2022" # specify date of prior analysis if desired. Form: MM_DD_YYYY

subsetsToDiagram <- c("noDxNoSSRI") # list of subset names to make diagrams for
roundsToDiagram <- c(1) # list of rounds to diagram

mediatorsToDiagram <- c("wmHyperintensities") # "occipital"
hemispheresTodiagram <- c("WB")  # possibilities: "WB", "L", "R"
metricsToDiagram <- c("vol") # possibilities: "vol", "area", "mThick"

outputPath <- "./outputs/figures/Figure1_medDia"

round1_MedDiagrams(dataPath = dataPath,
                   dataFileName = dataFileName, 
                   dataFileInfix = infix,
                   dataDate = dataDate,
                   subsets = subsetsToDiagram,
                   rounds = roundsToDiagram,
                   rois = mediatorsToDiagram,
                   hemispheres = hemispheresTodiagram,
                   metrics = metricsToDiagram,
                   pCor = T, 
                   outputPath = outputPath)

### Manual diagrams
# medDigram_a <- med_diagram(mediation_data = data.frame(
#   main = "",  
#   lab_x   = "C-Reactive Protein",
#   lab_m   = "L Cerebellum cortex volume",
#   lab_y   = "Fluid Intelligence Score",
#   coef_xm = "-0.0XX **",
#   coef_my = "-0.0XX",
#   coef_xy = "c\` = -0.0XX *** \n ab = -0.0XX"), save = T, outputName = "test.png")

## 4.B Figure 2: ggSegmentation ----
### 4.B.i. Functions ----
ggSegFormatData <- function(dataPath, dataFileName, dataInfix, dataDate, subset_i, round_j){
  require(readr)
  require(dplyr)
  require(glue)

  df <- tibble(read_csv(glue("{dataPath}/{dataFileName}{dataInfix}_{dataDate}.csv"), show_col_types = FALSE)) # Import data
  df <- df %>% filter(subset == subset_i, round == round_j) # remove irrelevant rows
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
  cat("Done formatting")
  return(df_formatted)
} # Format mediation summary file to appropriate format for 'createBrainFigure' function. I.e., one row per ROI, metric column, proper names for regions. Input: dataframe with results for all cognitive outcomes of interest. Returns list of dataframes, one per cognitive outcome of interest. 'dataPath' string of path to folder containing output of mediation analyses; 'dataFileName' string of name of df for specific data subset containing mediation analyses with corrected p-value for indirect effect. 'dataFileInfix': infix of data file name, include any seperators between the file name and the infix; 'dataDate': date of data file; 'subset' - name of subset to show results from; 'round' - round of mediators to show on figure. N.B. the following function only accepts DKT variables, like those specified in round 2.

createBrainFigure <- function(file, figureValue, savePath, Y, subset){
  require(stringr)
  require(ggplot2)
  require(ggseg)
  require(ggthemes)
  require(tidyr)
  
  # Extract name of cognitive performance variable ----
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
    value <- glue("{Y}.IndEff_est")
    figureName <- glue("{subset}_{Yname_file}_indEff_est")
    title <- "Indirect effect"
    
    lowBound <- -.002 # change as necessary
    hiBound <- .002 # change as necessary
    gradientCodes <- c("#062635", "#0A4C6A", "#12719E",  "#1696D2", "#46ABDB",  "#73BFE2", "#A2D4EC", "#CFE8F3", "#A2D4EC", "#73BFE2","#46ABDB", "#1696D2", "#12719E", "#0A4C6A", "#062635") # Higher is better; full colour scale: c("#062635", "#0A4C6A", "#12719E", "#1696D2", "#46ABDB","#73BFE2","#A2D4EC","#CFE8F3"),
    # breakValues <- c(-.1, -.05, -.01, -.001, 0, .001, .01, .05, .1) # !!! PROPER BOUNDS TO BE DETERMINED
    # breakValueLabels <- c("-.1", ".05", "-.01", "-.001", "0", ".001", ".01", ".05", ".1")
    
    gradientSettings <- scale_fill_gradientn(
      name = valueName,
      colours = gradientCodes,
      limits = c(lowBound,hiBound),
      guide = "colourbar",
      # values = c(lowBound, breakValues, hiBound),
      # breaks = breakValues,
      # labels = breakValueLabels # Make scale for break values defined in if statements above
    )
    
  } else {
    
    return(cat("Error: ", figureValue, " is invalid. \n Options are 'p', 'pCor', or 'indEff'"))
    
  }
  
  cat("\n --------- \n Value:", value, "\n")
  
  regionColNum <- which(colnames(file) == glue("{Y}.Region"))
  
  # I. Extract DKT brain regions ---------
  ## I.a. makes a list of all ROIs to keep for DKT figures. ----
  DKTRegions <- brain_labels(dk) # name of all DKT labels
  DKT_matches <- c()
  for(roi in DKTRegions){
    new_matches <- file[[regionColNum]][grepl(roi, file[[regionColNum]], ignore.case = TRUE, fixed = FALSE) == TRUE]
    DKT_matches <- c(DKT_matches, unique(new_matches))
  } 
  DKT_matches <- unique(DKT_matches)
  # str(matches)
  
  # II. Extract rows with value of column 'Region' specified in matches. ----
  reallyShort_DK <- subset(file, file[[regionColNum]] %in% DKT_matches) # Extracts all columns.
  # str(reallyShort_DK)
  
  # III. Make and format df with only necessary variables. ----
  valueColNum <- which(colnames(reallyShort_DK) == glue("{Y}.{value}"))
  regionColNum <- which(colnames(reallyShort_DK) == glue("{Y}.Region"))
  metricColNum <- which(colnames(reallyShort_DK) == glue("{Y}.Metric"))
  
  # cat("- Debugging - \n \t Structure: \n ", str(file[valueColNum]), "\n \t Table: \n ", table(file[valueColNum]))
  
  reallyShort_DK <- reallyShort_DK[,c(regionColNum, metricColNum, valueColNum)] 
  colnames(reallyShort_DK) <- c("label", "measure", "dfValue")
  reallyShort_DK$label <- as.factor(reallyShort_DK$label)
  # View(reallyShort_DK)
  # str(reallyShort_DK)
  
  # IV. Make figures ----
  cortical_pos <- c("left lateral", "left medial", "right medial", "right lateral")
  ## IV.a. DKT figure ----
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
      plot.background = element_rect(fill = "white", linewidth = 0),
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
  
  (cat("Saving DKT figure: ", glue("{figureName}_DK_{getDate()}.jpg"), "\n") + ggsave(glue("{savePath}/{figureName}_DK_{getDate()}.jpg"), dpi=450))
} # Creates ggSegmentation for DSK atlas. Takes: 'file' a data file formatted appropriately (see ROIS.csv for reference); 'figureValue' string specifying the values the figure is to represent (either 'p', 'indEff', or 'pCor'); 'savePath' path where figure should be saved; 'Y' dependent variable name

#### Call functions -----

dataFilePath <- "outputs/med" # to MASTER file
round = 2 # round 2 has the specific brain regions
dataFileName <- "MASTER"
infix <- "_med_C_pCor"
dataDate <- getDate() # by default will take today's date
# dataDate <- "11_05_2022" # specify date of prior analysis if desired. Form: MM_DD_YYYY
subsetNames <- c("noDxNoSSRI")

figureValueToPlot <- "indEff" # either 'p', 'indEff', or 'pCor'

outputPath <- "./outputs/figures/Figure2_ggSeg"

for(subset in subsetNames){
  cat("\n---- ", subset, "-----------------------\n")
  # load original file for that subset
  formattedFile <- ggSegFormatData(dataPath = dataFilePath, dataFileName = dataFileName, dataInfix = infix, dataDate = dataDate, subset_i = subset, round_j = round)
  
  for(i in 1:length(formattedFile)){
    cogVar <- names(formattedFile[i])
    createBrainFigure(file = as.data.frame(formattedFile[i]), figureValue = figureValueToPlot, savePath = outputPath, Y = cogVar, subset = subset)
  }
} # generate figures for each subset and for each cognitive variable

# 5. Varia ----
## Print session information ----
sessionInfo()
capture.output(sessionInfo(), file = paste("rSessionInfo_DataPrep_", getDate(), ".csv", sep = ""))
      
