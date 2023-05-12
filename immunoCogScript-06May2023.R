###
# This file contains is the analysis script for the immunology, cognition and brain mediation analysis of the UKBB.
# Prepared by Daniel Mendelson with help from Katie M. Lavigne and Joshua Unrau.
# Work performed under the supervision of Dr. Martin Lepage.
# 06 May 2023
###

# Preamble -----
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

library(arsenal)
library(finalfit)
library(Hmisc)
library(coin)
library(esc)
library(glmnet)
library(car)
library(gridExtra)

# 0 - Functions ----
getDate <- function(){
  require(stringr)
  return(str_c(format(Sys.time(), "%m"), "_", format(Sys.time(), "%d"), "_", format(Sys.time(), "%Y")))
}

## Assumption checks ----
### Normality ----
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

## Compare non-normal vairables
runWilcox <- function(varsOfInterest, df, by, outputDir, fileName){
  require(rcompanion)
  require(stringr)
  require(glue)
  
  # wilcoxOutput <- data.frame("value type" = c("statistic","p", "n1", "n2", "r(effSize)"), row.names = 1) # initialize dataframe for the output
  wilcoxOutput <- data.frame("value type" = c("statistic","p", "r(effSize)"), row.names = 1) # initialize dataframe for the output
  i <- integer(0)
  for(i in varsOfInterest){ # iterate through the variables of interest. Run wilcoxon-paired rank test. Put results into a dataframe 
    # print(i)

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
    wilcoxInterest <- c(wilcoxStat, wilcoxP, r) # combine W statistic, p, n1 and n2 into a vector
    wilcoxOutput[[i]] <- with(wilcoxOutput, wilcoxInterest) # Add the results for variable {i} to the output dataframe
    # print(wilcoxInterest)
    }
  }  
  # return(wilcoxonOutput)
  write.csv(wilcoxOutput, file = glue("{outputDir}/{fileName}_{getDate()}.csv")) #exports wilcoxOutput to a CSV file
  print(paste("The Wilcoxon comparison has been saved as:", glue("{fileName}_{getDate()}.csv")))
} # "varsOfInterest" - list of variable names to compare between groups; "df" data frame with all vars of interest and grouping variable; "by" - grouping variable; "outputDir" - path to output directory; "fileName" - desired name of output file

## Compare kept cases to excluded cases 

compareExcluded <- function(df_subset, varsToCompare, grp, subset, simplified){
  require(DescTools)
  require(effsize)
  date <- getDate() # set todays date for easier output filenaming
  
  if(simplified == T){
    output <- TOne(df_subset[varsToCompare], grp = df_subset[[grp]], add.length = T, total = T, 
                   FUN = function(x) gettextf("%s (%s)",
                                              Format(mean(x, na.rm = T), digits = 2, big.mark = ","),
                                              Format(sd(x, na.rm = T), digits = 2, big.mark = ",")),
                   TEST = list(
                     num  = list(fun = function(x, g){paste(
                       "F(", 
                       summary(aov(x ~ g))[[1]][1, "Df"], ",",
                       summary(aov(x ~ g))[[1]][2, "Df"], ") = ", 
                       formatC(summary(aov(x ~ g))[[1]][1, "F value"], digits = 5), 
                       if(summary(aov(x ~ g))[[1]][1, "Pr(>F)"] < .001){
                         ", p <.001"
                       } else {
                         paste(", p = ", formatC(summary(aov(x ~ g))[[1]][1, "Pr(>F)"], format = "g", digits = 3), sep = "")
                       }, 
                       ", g = ", formatC(effsize::cohen.d(x ~ g, data = df, pooled = T, hedges.correction = T)$estimate, format = "g", digits = 3), # calculate g. See Durlak 2009, J. Ped. Psyc.
                       sep = "")},
                       lbl = "ANOVA"),
                     cat  = list(fun = function(x, g){paste(
                       "χ²(", chisq.test(table(x, g))$parameter, ") = ", 
                       formatC(chisq.test(table(x, g))$statistic, digits = 5), 
                       if(chisq.test(table(x, g))$p.val < .001){
                         ", p <.001"
                       } else {
                         paste(", p = ", formatC(chisq.test(table(x, g))$p.val, format = "g", digits = 3), sep = "")
                       },
                       ", v = ", formatC(sqrt(chisq.test(table(x, g))$statistic)/(length(x)*chisq.test(table(x, g))$parameter), format = "g", digits = 2), 
                       sep = "")},
                       lbl = "Chi-Square test"),
                     dich = list(fun = function(x, g){paste(
                       "χ²(", chisq.test(table(x, g))$parameter, ") =", 
                       formatC(chisq.test(table(x, g))$statistic, digits = 2), 
                       if(chisq.test(table(x, g))$p.val < .001){
                         ", p <.001"
                       } else {
                         paste(", p = ", formatC(chisq.test(table(x, g))$p.val, format = "g", digits = 3), sep = "")
                       },
                       ", v = ", formatC(sqrt(chisq.test(table(x, g))$statistic)/(length(x)*chisq.test(table(x, g))$parameter), format = "g", digits = 2), 
                       sep = "")},
                       lbl = "Chi-Square test")), 
                   fmt = list(abs  = Fmt("abs"), 
                              num  = Fmt("num"), 
                              per  = Fmt("per"),
                              pval = as.fmt(fmt = "p*")))
    fileName <- glue("CompareMissingCases_TOne_{subset}_simple_{date}.csv")
  } else {
    output <- TOne(df_subset[varsToCompare], grp = df_subset[[grp]], add.length = T, total = T, 
                   FUN = function(x) gettextf("%s (%s); %s (%s)",
                                              Format(mean(x, na.rm = T), digits = 2, big.mark = ","),
                                              Format(sd(x, na.rm = T), digits = 2, big.mark = ","),
                                              Format(median(x, na.rm = T), digits = 2, big.mark = ","),
                                              Format(IQR(x, na.rm = T), digits = 2, big.mark = ",")),
                   TEST = list(
                     num  = list(fun = function(x, g){paste(
                       "F(", 
                       summary(aov(x ~ g))[[1]][1, "Df"], ",",
                       summary(aov(x ~ g))[[1]][2, "Df"], ") = ", 
                       if(summary(aov(x ~ g))[[1]][1, "Pr(>F)"] < .001){
                         ", p <.001"
                       } else {
                         paste(", p = ", formatC(summary(aov(x ~ g))[[1]][1, "Pr(>F)"], format = "g", digits = 3), sep = "")
                       }, 
                       ", p = ", formatC(summary(aov(x ~ g))[[1]][1, "Pr(>F)"], format = "g", digits = 3), 
                       ", g = ", formatC(effsize::cohen.d(x ~ g, data = df, pooled = T, hedges.correction = T)$estimate, format = "g", digits = 3), # calculate g. See Durlak 2009, J. Ped. Psyc.
                       sep = "")},
                       lbl = "ANOVA"),
                     cat  = list(fun = function(x, g){paste(
                       "χ²(", chisq.test(table(x, g))$parameter, ") = ", 
                       formatC(chisq.test(table(x, g))$statistic, digits = 5), 
                       if(chisq.test(table(x, g))$p.val < .001){
                         ", p <.001"
                       } else {
                         paste(", p = ", formatC(chisq.test(table(x, g))$p.val, format = "g", digits = 3), sep = "")
                       },
                       ", v = ", formatC(sqrt(chisq.test(table(x, g))$statistic)/(length(x)*chisq.test(table(x, g))$parameter), format = "g", digits = 2), 
                       sep = "")},
                       lbl = "Chi-Square test"),
                     dich = list(fun = function(x, g){paste(
                       "χ²(", chisq.test(table(x, g))$parameter, ") =", 
                       formatC(chisq.test(table(x, g))$statistic, digits = 2), 
                       if(chisq.test(table(x, g))$p.val < .001){
                         ", p <.001"
                       } else {
                         paste(", p = ", formatC(chisq.test(table(x, g))$p.val, format = "g", digits = 3), sep = "")
                       },
                       ", v = ", formatC(sqrt(chisq.test(table(x, g))$statistic)/(length(x)*chisq.test(table(x, g))$parameter), format = "g", digits = 2), 
                       sep = "")},
                       lbl = "Chi-Square test")), 
                   fmt = list(abs  = Fmt("abs"), 
                              num  = Fmt("num"), 
                              per  = Fmt("per"),
                              pval = as.fmt(fmt = "p*")))
    fileName <- glue("CompareMissingCases_TOne_{subset}_extensive_{date}.csv")
  }
  write.csv(output, file = fileName)
  print(paste("The file ", fileName, " has been saved.", sep = ""))
  return(output)
} # df: dataframe with all observations, all variables of interest (listed in 'varsToCompare') and the grouping variable; varsToCompare: list of variables that should be compared between groups; grp: name of grouping variable; subset: name of subset for use in output file naming; simplified: logical indicating if the comprehensive table or a shortened version should be produced.

setwd("C:/Users/katie/OneDrive - McGill University/publishing/_rev_Mendelson_UKBB/res")

# 1 - Biobank data analysis----

## Data Preparation -----
df_all <- read_csv("../data/RawData/Daniel_2022-05-02.csv") # import df
eidToRemove <- read_csv("../data/RawData/w45551_2023-04-25.csv", col_names = "eid") # file with EID of participants who withdrew consent from UKBB

df_all <- df_all %>% filter(!(eid %in% eidToRemove)) # remove rows with eid in this file
write.csv(df_all, "../data/RawData/Daniel_2022-05-02_ProperEID.csv", row.names=FALSE)
rm(df_all, eidToRemove) # remove dataframes not in use

date <- str_c(format(Sys.time(), "%m"), "_", format(Sys.time(), "%d"), "_", format(Sys.time(), "%Y")) # set todays date for easier output filenaming
df_all <- read_csv("../data/RawData/Daniel_2022-05-02_ProperEID.csv", guess_max = 10000) # import df

cat(paste("There are ", sum(is.na(df_all$date_assess2_t2)), "cases that have not completed time point 2. These cases will be removed."))
df_all <- df_all %>%
  filter(!(is.na(df_all$date_assess2_t2))) # remove rows without timepoint 2 assessment date

# df_small <- df_all[1:10000,]
df <- df_all # specifies what df to use. In final analyses will want to use df_all but want to use smaller df while developing code.
# rm(df_small) # remove dataframes not in use
# sort(colnames(df))
# View(df)

### Standardize categorical vars ----
# standardize categorical variable responses

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

### Add columns ----

# Waist:hip ratio
df <- df %>% 
  mutate(weight_waistToHip_t0 = weight_waistCirc0_t0/weight_hipCirc0_t0) %>% 
  mutate(weight_waistToHip_t2 = weight_waistCirc2_t2/weight_hipCirc2_t2)

# CRP log transform
df <- df %>% 
  mutate(crp_log = log(crp_aliquot_t0))

# Diagnosis - any of interest
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

# Medication - any of interest
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
  if(n_distinct(df[i]) == 1){
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

# Num days between assessment 2 and assessment 0
df <- df %>% 
  mutate(demo_daysBtwAssess = as.numeric(df$date_assess2_t2 - df$date_assess0_t0, units = "days"))

# saveRDS(df, file = glue("./Data/Processed/df_preLobes_{getDate()}.rds"))
# rm(df)
df <- df %>% 
  mutate(crp_hourCollected = hour(df$crp_timeCollected_t0)) %>% 
  mutate(cog_hourCompleted = hour(df$cog_timeCompleted_t2)) %>% 
  mutate(brain_hourCompleted = hour(df$brain_timeCompleted_t2))

### Make columns appropriate data types ----

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

### Take mean of covariates measured at both time points ----
df <- df %>% 
  mutate(age_mean02 = (demo_age_assess0_t0 + demo_age_assess2_t2)/2) %>% 
  mutate(weight_waistToHip_mean02 = (weight_waistToHip_t0 + weight_waistToHip_t2)/2) %>%
  mutate(sleep_duration_mean02 = (sleep_duration0_t0 + sleep_duration2_t2)/2)

### Compute and add lobular brain volumes ----

# dataDate <- "04_15_2023"
# df <- readRDS(glue("./Data/Processed/df_preLobes_{dataDate}.rds"))
vars_BrainRegion <- read_csv("./Variables/brainVarsbyLobe.csv")
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

### Save cleaned DF ----
df_interim <- df
saveRDS(df_interim, file = glue("./Data/Processed/df_full_{getDate()}.rds"))
rm(df_interim)
# df <- df_interim

dataDate <- getDate()
# dataDate <- "12_14_2022" # For testing purposes. Format: MM_DD_YYYY
df <- readRDS(glue("./Data/Processed/df_full_{dataDate}.rds"))

## Exclusions ----
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


### Remove missing ----
essentialVars <- c(finalCovars, completionTimeVars, crpAliqVars, brainMorphVars) # list of variable names that must be complete in order for case to be retained
# essentialVars <- readRDS("revisedEssentialVarsAnalysis.rds")

df <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    if_any(.cols = essentialVars, function(x)(
      is.na(x))) ~ TRUE,
    TRUE ~ FALSE)))

saveRDS(df, file = glue("./Data/Processed/df_PreSubset_{getDate()}.rds"))
df_removed <- df %>%
  filter(if_any(any_of(essentialVars), function(x)(is.na(x)))) # make new df of cases with missing values

numMissing <- matrix(ncol=2, nrow = length(essentialVars), dimnames = list(colnames(df_removed[essentialVars]), c("NACount", "NotMissingOnThisVar")))

for(i in 1:length(essentialVars)){
  numMissing[i,] <- c(sum(is.na(df_removed[essentialVars[i]])), sum(!is.na(df_removed[essentialVars[i]])))
}

numMissing_df <- as.data.frame(numMissing) %>% arrange(desc(NACount))
write.csv(numMissing_df, file = paste("UKBB_sourceOfNAValues_", date, ".csv", sep = ""))

### Remove CRP > 10 ----
cat(sum(df$crp_aliquot_t0 > 10, na.rm = T), " cases to remove given CRP > 10.", sep = "")
df <- df %>%
  mutate(missingEssentialVar = if_else(df$crp_aliquot_t0 > 10, TRUE, FALSE, missing = TRUE)) # mark CRP > 10 to remove


## Assumption checks ----
### Normality ----
varsToCompare <- c("demo_sex_t0", "demo_birthYear_t0", "demo_ethnicity_t0", "demo_age_assess0_t0", "demo_age_assess2_t2", "age_mean02", "ses_townsend_t0", "crp_aliquot_t0", "crp_log", "med_AnyOfInterest0", "med_AnyOfInterest2", "med_AnyOfInterest0_excludingSSRI", "med_AnyOfInterest2_excludingSSRI", "med_SSRI_t0","med_SSRI_t2", "dx_AnyOfInterest", "dx_Dementia", "dx_SSD", "dx_MoodDisorder", "smoke_currently0_t0", "smoke_currently2_t2", "exercise_IPAQActivityGroup_t0", "weight_BMI0_t0","weight_BMI2_t2", "weight_waistToHip_t0","weight_waistToHip_t2","weight_waistToHip_mean02", "sleep_duration0_t0", "sleep_duration2_t2", "sleep_duration_mean02", "menopause0_t0", "menopause2_t2", "hand_t0", "crp_hourCollected", "crp_timeFasting_t0", "crp_aliquotdelay",  "diet_water_t0", "diet_alc_freq0_t0", "diet_alc_freq2_t2", "diet_alcohol_yesterdayIntake0_t0", "diet_alcohol_yesterdayIntake2_t2") # Add variables to compare between excluded and retained participants

varsToCompare_ColNum <- which(colnames(df) %in% varsToCompare)

for(i in c(medVars0NoCases, medVars2NoCases, dxVarsNoCases)){
  if(is.null(i) == FALSE){
    if(i %in% varsToCompare_ColNum){
      colNum <- which(colnames(df) == colnames(df[i]))
      varsToCompare_ColNum <- varsToCompare_ColNum[-which(varsToCompare_ColNum == colNum)]
    }
  }
}

# create list of non-factor numeric variables for normality check
numVarsToCompare <- lapply(df[varsToCompare_ColNum], numNotFactor)
df_numVarsToCompare <- as.data.frame(do.call(cbind, numVarsToCompare))
numVarColsToCompare <- (which(colnames(df) %in% colnames(df_numVarsToCompare)))

# For each variable in this list, determine if normal
normalityPlots <- mapply(normalCheck, df[c(numVarColsToCompare)], colName = colnames(df[c(numVarColsToCompare)]), colNum = c(numVarColsToCompare))
# N.b. QQ plot should follow a straight line
# nonParametricVars <- c(3,4,5,9,15,16, 17,19,36,38,41,43,49,56,57,58,66,67,69,101) # vector specifying which values of dependent to run non-parametric tests on 
nonParaToCompare <- c(3,10,11,14,15,18,19,20,21,22,25,47,48,49,50,713,719,720,762)

## Subset ----
# df <- read_csv("./data/forAnalysis/UKBB_dfForAnal_All_05_05_2022.csv", lazy = T) # import df if desired
subsetData_outputPath <- "./Data/Processed/Subsets"
subsetNames <- c("noMedNoDx", "all", "oldest","onlyDx","noDx","onlyMed","noMed","onlySSRI","youngest","noDxNoSSRI","noSSRI")

### Subset 1: no further exclusions (only missing data and CRP > 10) ----
comparisonExcluded_all <- compareExcluded(df, varsToCompare, grp = "missingEssentialVar", subset = "all", simplified = T) 
#runWilcox(df = df, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "all_wilcox")
df_all <- df %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_all_", getDate(), sep = "")
write.csv(df_all, file = glue("{subsetData_outputPath}/{fileName}.csv")) # exports cleaned dataframe with all variables ready for analysis
cat("Subset only with complete cases and CRP < 10 saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_all))

### Subset 2: oldest tertile ----
tertCutOff <- quantile(df$demo_age_assess0_t0, c(.66))
df_oldest <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    demo_age_assess0_t0 >= tertCutOff ~ FALSE,
    TRUE ~ TRUE)))

# table(df_oldest$missingEssentialVar)
# str(df_oldest$missingEssentialVar)
comparisonExcluded_oldest <- compareExcluded(df_oldest, varsToCompare, grp = "missingEssentialVar", subset = "oldestTert", simplified = T) 
#runWilcox(df = df_oldest, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "oldest_wilcox")

df_oldest <- df_oldest %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_oldest_", getDate(), sep = "")
write.csv(df_oldest, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset with the oldest tertile of cases at instrance 0 (", tertCutOff, ") saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_oldest))

### Subset 3: youngest tertile ----
tertCutOff <- quantile(df$demo_age_assess0_t0, c(.33))
df_youngest <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    demo_age_assess0_t0 <= tertCutOff ~ FALSE,
    TRUE ~ TRUE)))
comparisonExcluded_youngest <- compareExcluded(df_youngest, varsToCompare, grp = "missingEssentialVar", subset = "youngestTert", simplified = T) 
#runWilcox(df = df_youngest, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "youngest_wilcox")

df_youngest <- df_youngest %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_youngest_", getDate(), sep = "")
write.csv(df_youngest, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset with the youngest tertile of cases at instrance 0 (", tertCutOff, ") saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_youngest))

### Subset 4: Exclude Diagnoses -------
df_noDx <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    dx_AnyOfInterest == 0 ~ FALSE,
    TRUE ~ TRUE)))
comparisonExcluded_noDx <- compareExcluded(df_noDx, varsToCompare, grp = "missingEssentialVar", subset = "NoDx", simplified = T)
#runWilcox(df = df_noDx, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "noDx_wilcox")

df_noDx <- df_noDx %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_noDx_", getDate(), sep = "")
write.csv(df_noDx, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset excluding those with diagnosis of interest saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_noDx))

### Subset 5: Only Diagnoses ------------
df_onlyDx <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    dx_AnyOfInterest == 1 ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_onlyDx <- compareExcluded(df_onlyDx, varsToCompare, grp = "missingEssentialVar", subset = "OnlyDx", simplified = T)
#runWilcox(df = df_onlyDx, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "onlyDx_wilcox")

df_onlyDx <- df_onlyDx %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_onlyDx_", getDate(), sep = "")
write.csv(df_onlyDx, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset only those with a diagnosis of interest saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_onlyDx))

### Subset 6: Exclude medications -----------
df_noMed <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    med_AnyOfInterest0 == 0 & med_AnyOfInterest2 == 0 ~ FALSE,
    TRUE ~ TRUE)))
comparisonExcluded_noMed <- compareExcluded(df_noMed, varsToCompare, grp = "missingEssentialVar", subset = "NoMed", simplified = T)
#runWilcox(df = df_noMed, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "noMed_wilcox")

df_noMed <- df_noMed %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_noMed_", getDate(), sep = "")
write.csv(df_noMed, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset excluding those using medications of interest saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_noMed))

### Subset 7: Only medications ----
df_onlyMed <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    med_AnyOfInterest0 == 1 | med_AnyOfInterest2 == 1 ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_onlyMed <- compareExcluded(df_onlyMed, varsToCompare, grp = "missingEssentialVar", subset = "OnlyMed", simplified = T)
#runWilcox(df = df_onlyMed, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "onlyMed_wilcox")

df_onlyMed <- df_onlyMed %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_onlyMed_", getDate(), sep = "")
write.csv(df_onlyMed, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset only with those using medications of interest saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_onlyMed))

### Subset 8: Exclude SSRI -----------
df_noSSRI <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    med_SSRI_t0 == 0 & med_SSRI_t2 == 0 ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_noSSRI <- compareExcluded(df_noSSRI, varsToCompare, grp = "missingEssentialVar", subset = "NoSSRI", simplified = T)
#runWilcox(df = df_noSSRI, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "noSSRI_wilcox")

df_noSSRI <- df_noSSRI %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_noSSRI_", getDate(), sep = "")
write.csv(df_noSSRI, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset excluding those using SSRIs saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_noSSRI))

### Subset 9: Only SSRI ----
df_onlySSRI <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    med_SSRI_t0 == 1 | med_SSRI_t2 == 1 ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_onlySSRI <- compareExcluded(df_onlySSRI, varsToCompare, grp = "missingEssentialVar", subset = "onlySSRI", simplified = T)
#runWilcox(df = df_onlySSRI, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "onlySSRI_wilcox")

df_noSSRI <- df_noSSRI %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_onlySSRI_", getDate(), sep = "")
write.csv(df_noSSRI, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset onlly with those using SSRIs saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_onlySSRI))

### Subset 10: exclude Med and Dx --------
df_noMedNoDx <- df %>%
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    dx_AnyOfInterest == 0 & med_AnyOfInterest0 == 0 & med_AnyOfInterest2 == 0 ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_noMedNoDx <- compareExcluded(df_noMedNoDx, varsToCompare, grp = "missingEssentialVar", subset = "NoMedNoDx", simplified = T)
#runWilcox(df = df_noMedNoDx, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "noMedNoDx_wilcox")

df_noMedNoDx <- df_noMedNoDx %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_noMedNoDx_", getDate(), sep = "")
write.csv(df_noMedNoDx, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset excluding those with diagnoses or using any medication of interest saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_noMedNoDx))

### Subset 11: exclude Dx and SSRI --------
df_noDxNoSSRI <- df %>%
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    dx_AnyOfInterest == 0 & med_SSRI_t0 == 0 & med_SSRI_t2 == 0 ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_noDxNoSSRI <- compareExcluded(df_noDxNoSSRI, varsToCompare, grp = "missingEssentialVar", subset = "NoDxNoSSRI", simplified = T)
#runWilcox(df = df_noDxNoSSRI, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "noDxNoSSRI_wilcox")

fileName <- paste("df_noDxNoSSRI_", getDate(), sep = "")
df_noDxNoSSRI <- df_noDxNoSSRI %>% filter(missingEssentialVar == FALSE)
write.csv(df_noDxNoSSRI, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset excluding those with diagnoses or using SSRIs saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_noDxNoSSRI))

# Describe Subsets ----
path <- subsetData_outputPath
#path <- "./Data/Subsets"
subsetNames <- c("noMedNoDx", "all", "oldest","onlyDx","noDx","onlyMed","noMed","onlySSRI","youngest","noDxNoSSRI","noSSRI")
dataDate <- getDate()
#dataDate <- ""

for(subset in subsetNames){
  fileName <- glue("{path}/df_{subset}_{dataDate}.csv")
  subset_df <- read_csv(fileName, lazy = T, show_col_types = F)
  print(subset)
  write.csv(psych::describe(subset_df), file = glue("./output/descriptive/{subset}_summary_{getDate()}.csv"))
  rm(subset_df)
}

### Summarise
write.csv(df_retained, file = paste("./data/processed/dfRetained_", getDate(), ".csv", sep = ""))
df <- df_retained
varsToSummarise <- c(1:ncol(df))

numVars<- lapply(df, numNotFactor)
df_numVars <- as.data.frame(do.call(cbind, numVars))
numVarCols <- (which(colnames(df) %in% colnames(df_numVars)))
numVars <- c((which(colnames(df) %in% colnames(df_numVars))), which(colnames(df) == "date_assess0_t0"), which(colnames(df) == "date_assess2_t2"), which(colnames(df) == "cog_timeCompleted_t2"), which(colnames(df) == "brain_timeCompleted_t2"), which(colnames(df) == "crp_timeCollected_t0"), which(colnames(df) == "crp_assaydate_t0"))

for(i in c(medVars0NoCases, medVars2NoCases, dxVarsNoCases)){
  if(is.null(i) == FALSE){
    if(i %in% varsToSummarise){
      colNum <- which(colnames(df) == colnames(df[i]))
      varsToSummarise <- varsToSummarise[-which(varsToSummarise == colNum)]
    }
  }
}

write.csv(summary(df[-numVars]), file = paste("./outputs/descriptive/df_all_summary_fctrVars", getDate(), ".csv", sep = "")) # save counts of categorical variables

# Varia correlations
## Correl of CRP vars 
print("Hour collected")
median(hour(df$crp_timeCollected_t0), na.rm = T)
iqr(hour(df$crp_timeCollected_t0), na.rm = T)
cat("Correlation between CRP and hour collected: ")
cor.test(rank(df$crp_aliquot_t0), rank(hour(df$crp_timeCollected_t0)), na.action = na.omit, conf.level = .95)

print("Fasting time")
median(df$crp_timeFasting_t0, na.rm = T)
iqr(df$crp_timeFasting_t0)
cat("Correlation between CRP and hours fasted before blood draw: ")
cor.test(rank(df$crp_aliquot_t0), rank(df$crp_timeFasting_t0), na.action = na.omit)

cat("CRP aliquots were performed between", paste(min(df$crp_assaydate_t0))," and ", paste(max(df$crp_assaydate_t0)))
df$crp_assaydate_t0 <- as.POSIXct(df$crp_assaydate_t0)

cat("the mean time between blood sample and CRP measurement was (in years): \n\t", "mean:", mean(df$crp_aliquotdelay, na.rm = T), "SD: ", sd(df$crp_aliquotdelay, na.rm = T), "\n\t Median: ", median(df$crp_aliquotdelay, na.rm = T), "IQR: ", iqr(df$crp_aliquotdelay))

cor.test(rank(df$crp_aliquotdelay), rank(df$crp_aliquot_t0))

# 2 - Analyses ----------
## Analysis Functions -----
### Linear model expression constructor ----
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

### Correlation Analysis functions ----
runCor <- function(df, X, Y, M, covars, contrasts, YisFactor, name){
  expression <- createLMExpression(varNames = colnames(df), X = X, Y = Y, M = M, C = covars) # total model: predict Y from X and C  
  
  if(YisFactor){
    model <- glm(formula = eval(parse(text = expression)), data = df, contrasts = contrasts, family = "binomial", na.action = na.omit)
  } else {
    model <- lm(eval(parse(text = expression)), data = df, contrasts = contrasts, na.action = na.omit)
  }
  
  fileName <- glue("{name}_{getDate()}.txt")

  capture.output(
    cat(" -- Linear Model Analysis Output -- "),
    print(c("\n ------------------ \n Model expression: \n\t", expression)),
    cat(summary(model)),
    file = fileName)
  
  return(model)
}

corTxtFormat <- function(filePath, files, X, listOfCovars, logReg, outputPath){  
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
  write.csv(df_main, file = glue("{outputPath}/{subsetName}_corSummary_{date}.csv"))
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

### Mediation Analysis functions ----
runRMediation <- function(df, X, Y, M, covars, contrasts, YisFactor, name, alpha){
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
  
  a <- coef(aModel)[which(names(coef(aModel)) == X)] # want coef for effect of X on M
  a.se <- coef(summary(aModel))[which(names(coef(aModel)) == X),2] # want S.E. of the regression coef of X on M. N.b. index is [row of X, column of standard errors of estimate]
  b <- coef(bModel)[which(names(coef(bModel)) == M)] # want coef for effect of M on Y
  b.se <- coefficients(summary(bModel))[which(names(coef(bModel)) == M),2] # want coef for effect of M on Y
  
  fileName <- paste(name, "_", getDate(), ".txt", sep = "")
  
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
    file = fileName)
  
  return(fileName)
} # input: 'df' - dataframe containing all relevant variables; 'X' - name of the predictor variable; 'Y' - name of the outcome variable; 'M' - name of the mediator variable; 'C' - list of names of covariate variables; 'contrasts' - a list of ordinal variables specifying their contrasts; 'YisFactor' - if outcome is a factor, logical. If yes, logistic regression will be performed.; 'name' - name of analysis to be used for output filename; 'alpha' - critical p-value. Used for calculating confidence interval.

### Format RMediation function Outputs ----
formatMedOutput <- function(outputDF, filePath, fileName, analysisName, logReg){
  specificFilePath <- glue("{filePath}/{fileName}")
  fileObj <- readLines(specificFilePath)
  
  splits <- grep("~", fileObj)
  split1 <- splits[2]-2
  split2 <- splits[3]-2
  medOutput <- grep("Mediated Effect", fileObj)
  
  text_directEffect <- fileObj[2:split1]
  text_aPath <- fileObj[(split1+1):split2]
  text_bPath <- fileObj[(split2+1):(medOutput-1)]
  text_mediatedEffect <- fileObj[medOutput:length(fileObj)]
  models <- list(text_directEffect, text_aPath, text_bPath)
  
  variables <- c()
  values <- c()
  
  for(i in 1:(length(models))){ # iterate through all three models
    # Find name of mediator, indep var
    sectionText <- unlist(models[i])
    lmExpression <- grep("~", sectionText, value = TRUE) # from this line, extract analysis name
    lmExpression <- unlist(strsplit(lmExpression, "\\s+"))
    indepVar <- lmExpression[which(lmExpression == "~")-1] # n.b. indep var of this model may not be indep var of mediation analysis
    depVar <- lmExpression[which(lmExpression == "~")+1]
    
    depVarValues <- grep(depVar, sectionText, value = TRUE)
    if(length(depVarValues) == 2){
      depVarValues <- depVarValues[2]
    }
    
    depVarValues <- unlist(strsplit(depVarValues, "\\s+"))
    
    if(TRUE %in% grepl("<", depVarValues) && nchar(depVarValues[grep("<",depVarValues)]) == 1){
      pAndSymbol <- glue("{depVarValues[[5]]}{depVarValues[[6]]}")
      depVarValuesNoSigChar <- c(depVarValues[2:4], pAndSymbol)
    } else {
      depVarValuesNoSigChar <- depVarValues[2:5]
    }
    if(logReg == F || i == 2){ # the b path will always be linear regression as all brain morphology variables are continuous
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
      
      if(i == 2){
        variables <- append(variables, c(indepVar, depVar))
        values <- append(values, c(depVarValuesNoSigChar, R2, R2adj, FStat, df1, df2, p))
      } else {
        AIC <- "NA"
        variables <- append(variables, c(indepVar, depVar))
        values <- append(values, c(depVarValuesNoSigChar, R2, R2adj, FStat, df1, df2, p, AIC))
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
      
      variables <- append(variables, c(indepVar, depVar))
      values <- append(values, c(depVarValuesNoSigChar, R2, R2adj, FStat, df1, df2, p, AIC))
    }
  }
  
  indepVar <- variables[2]
  medVar <- variables[3]
  depVar <- variables[1]
  vars <- c(indepVar, medVar, depVar)
  name <- glue("{analysisName}_{medVar}_{depVar}")
  
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

## Run analyses -----
### Import data -----
path <- subsetData_outputPath
#path <- ""
fileName <- "df_noMedNoDx"
dataDate <- getDate()
#dataDate <- ""

df <- read.csv("./{path}/{fileName}_{dataDate}.csv") # import df

### Vars of interest ----
# sort(colnames(df))
crp_log <-  which(colnames(df) == "crp_log")
crp_log_z <- which(colnames(df) == "crp_log_z")
crp <- c(which(colnames(df) == "crp_aliquot_t0"), crp_log, crp_log_z)

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

lobularBrainMetrics <- c("area_insula_WB_t2", "vol_insula_WB_t2", "mThick_insula_WB_t2")
for(i in c("area", "mThick", "vol")){
  for(j in c("frontal", "parietal", "occipital", "temporal")){
    for(k in c("L", "R", "WB")){ 
      lobularBrainMetrics <- c(lobularBrainMetrics, glue("{i}_{j}_{k}_t2"))
      # print(glue("{i}_{j}_{k}"))
    }
  }
}
generalMetrics <- c("vol_BrainSeg_WB_t2","vol_BrainSegNotVent_WB_t2","vol_BrainSegNotVentSurf_WB_t2","vol_SubCortGray_WB_t2","vol_TotalGray_WB_t2","vol_SupraTentorial_WB_t2","vol_SupraTentorialNotVent_WB_t2","vol_EstimatedTotalIntraCranial_WB_t2","vol_BrainStem_WB_t2","vol_CSF_WB_t2","vol_wmHyperintensities_WB_t2","vol_WMhypointensities_WB_t2","vol_Hippocampus_L_t2","vol_Hippocampus_R_t2","vol_CerebellumCortex_R_t2","vol_CerebellumCortex_L_t2","vol_CerebellumWhiteMatter_L_t2","vol_CerebellumWhiteMatter_R_t2")

aseg_brainVars <- c("vol_VentricleChoroid_WB_t2","vol_3rdVentricle_WB_t2","vol_4thVentricle_WB_t2","vol_5thVentricle_WB_t2","vol_nonWMhypointensities_WB_t2","vol_OpticChiasm_WB_t2","vol_CCPosterior_WB_t2","vol_CCMidPosterior_WB_t2","vol_CCCentral_WB_t2","vol_CCMidAnterior_WB_t2","vol_CCAnterior_WB_t2","area_Accumbensarea_L_t2","vol_Cortex_L_t2","vol_CerebralWhiteMatter_L_t2","vol_LateralVentricle_L_t2","vol_InfLatVent_L_t2","vol_thalamusProper_L_t2","vol_Caudate_L_t2","vol_Putamen_L_t2","vol_Pallidum_L_t2","vol_Amygdala_L_t2","vol_Accumbensarea_L_t2","vol_VentralDC_L_t2","vol_choroidplexus_L_t2","area_Accumbensarea_R_t2","vol_Cortex_R_t2","vol_CerebralWhiteMatter_R_t2","vol_LateralVentricle_R_t2","vol_InfLatVent_R_t2","vol_thalamusProper_R_t2","vol_Caudate_R_t2","vol_Putamen_R_t2","vol_Pallidum_R_t2","vol_Amygdala_R_t2","vol_Accumbensarea_R_t2","vol_VentralDC_R_t2","vol_vessel_L_t2","vol_vessel_R_t2","vol_choroidplexus_R_t2") # all Freesurfer ASEG variables in UKBB

dkt_brainVars_area <- c("area_Caudalanteriorcingulate_L_t2","area_Caudalmiddlefrontal_L_t2","area_Cuneus_L_t2","area_Entorhinal_L_t2","area_Fusiform_L_t2","area_Inferiorparietal_L_t2","area_Inferiortemporal_L_t2","area_Isthmuscingulate_L_t2","area_Lateraloccipital_L_t2","area_Lateralorbitofrontal_L_t2","area_Lingual_L_t2","area_Medialorbitofrontal_L_t2","area_Middletemporal_L_t2","area_Parahippocampal_L_t2","area_Paracentral_L_t2","area_Parsopercularis_L_t2","area_Parsorbitalis_L_t2","area_Parstriangularis_L_t2","area_Pericalcarine_L_t2","area_Postcentral_L_t2","area_Posteriorcingulate_L_t2","area_Precentral_L_t2","area_Precuneus_L_t2","area_Rostralanteriorcingulate_L_t2","area_Rostralmiddlefrontal_L_t2","area_Superiorfrontal_L_t2","area_Superiorparietal_L_t2","area_Superiortemporal_L_t2","area_Supramarginal_L_t2","area_Transversetemporal_L_t2","area_Insula_L_t2", "area_Caudalanteriorcingulate_R_t2","area_Caudalmiddlefrontal_R_t2","area_Cuneus_R_t2","area_Entorhinal_R_t2","area_Fusiform_R_t2","area_Inferiorparietal_R_t2","area_Inferiortemporal_R_t2","area_Isthmuscingulate_R_t2","area_Lateraloccipital_R_t2","area_Lateralorbitofrontal_R_t2","area_Lingual_R_t2","area_Medialorbitofrontal_R_t2","area_Middletemporal_R_t2","area_Parahippocampal_R_t2","area_Paracentral_R_t2","area_Parsopercularis_R_t2","area_Parsorbitalis_R_t2","area_Parstriangularis_R_t2","area_Pericalcarine_R_t2","area_Postcentral_R_t2","area_Posteriorcingulate_R_t2","area_Precentral_R_t2","area_Precuneus_R_t2","area_Rostralanteriorcingulate_R_t2","area_Rostralmiddlefrontal_R_t2","area_Superiorfrontal_R_t2","area_Superiorparietal_R_t2","area_Superiortemporal_R_t2","area_Supramarginal_R_t2","area_Transversetemporal_R_t2","area_Insula_R_t2")

dkt_brainVars_vol <- c("vol_Caudalanteriorcingulate_R_t2","vol_Caudalmiddlefrontal_R_t2","vol_Cuneus_R_t2","vol_Entorhinal_R_t2","vol_Fusiform_R_t2","vol_Inferiorparietal_R_t2","vol_Inferiortemporal_R_t2","vol_Isthmuscingulate_R_t2","vol_Lateraloccipital_R_t2","vol_Lateralorbitofrontal_R_t2","vol_Lingual_R_t2","vol_Medialorbitofrontal_R_t2","vol_Middletemporal_R_t2","vol_Parahippocampal_R_t2","vol_Paracentral_R_t2","vol_Parsopercularis_R_t2","vol_Parsorbitalis_R_t2","vol_Parstriangularis_R_t2","vol_Pericalcarine_R_t2","vol_Postcentral_R_t2","vol_Posteriorcingulate_R_t2","vol_Precentral_R_t2","vol_Precuneus_R_t2","vol_Rostralanteriorcingulate_R_t2","vol_Rostralmiddlefrontal_R_t2","vol_Superiorfrontal_R_t2","vol_Superiorparietal_R_t2","vol_Superiortemporal_R_t2","vol_Supramarginal_R_t2","vol_Transversetemporal_R_t2","vol_Insula_R_t2", "vol_Caudalanteriorcingulate_L_t2","vol_Caudalmiddlefrontal_L_t2","vol_Cuneus_L_t2","vol_Entorhinal_L_t2","vol_Fusiform_L_t2","vol_Inferiorparietal_L_t2","vol_Inferiortemporal_L_t2","vol_Isthmuscingulate_L_t2","vol_Lateraloccipital_L_t2","vol_Lateralorbitofrontal_L_t2","vol_Lingual_L_t2","vol_Medialorbitofrontal_L_t2","vol_Middletemporal_L_t2","vol_Parahippocampal_L_t2","vol_Paracentral_L_t2","vol_Parsopercularis_L_t2","vol_Parsorbitalis_L_t2","vol_Parstriangularis_L_t2","vol_Pericalcarine_L_t2","vol_Postcentral_L_t2","vol_Posteriorcingulate_L_t2","vol_Precentral_L_t2","vol_Precuneus_L_t2","vol_Rostralanteriorcingulate_L_t2","vol_Rostralmiddlefrontal_L_t2","vol_Superiorfrontal_L_t2","vol_Superiorparietal_L_t2","vol_Superiortemporal_L_t2","vol_Supramarginal_L_t2","vol_Transversetemporal_L_t2","vol_Insula_L_t2")

dkt_brainVars_mThick <- c("mThick_Caudalanteriorcingulate_L_t2","mThick_Caudalmiddlefrontal_L_t2","mThick_Cuneus_L_t2","mThick_Entorhinal_L_t2","mThick_Fusiform_L_t2","mThick_Inferiorparietal_L_t2","mThick_Inferiortemporal_L_t2","mThick_Isthmuscingulate_L_t2","mThick_Lateraloccipital_L_t2","mThick_Lateralorbitofrontal_L_t2","mThick_Lingual_L_t2","mThick_Medialorbitofrontal_L_t2","mThick_Middletemporal_L_t2","mThick_Parahippocampal_L_t2","mThick_Paracentral_L_t2","mThick_Parsopercularis_L_t2","mThick_Parsorbitalis_L_t2","mThick_Parstriangularis_L_t2","mThick_Pericalcarine_L_t2","mThick_Postcentral_L_t2","mThick_Posteriorcingulate_L_t2","mThick_Precentral_L_t2","mThick_Precuneus_L_t2","mThick_Rostralanteriorcingulate_L_t2","mThick_Rostralmiddlefrontal_L_t2","mThick_Superiorfrontal_L_t2","mThick_Superiorparietal_L_t2","mThick_Superiortemporal_L_t2","mThick_Supramarginal_L_t2","mThick_Transversetemporal_L_t2","mThick_Insula_L_t2","mThick_Caudalanteriorcingulate_R_t2","mThick_Caudalmiddlefrontal_R_t2","mThick_Cuneus_R_t2","mThick_Entorhinal_R_t2","mThick_Fusiform_R_t2","mThick_Inferiorparietal_R_t2","mThick_Inferiortemporal_R_t2","mThick_Isthmuscingulate_R_t2","mThick_Lateraloccipital_R_t2","mThick_Lateralorbitofrontal_R_t2","mThick_Lingual_R_t2","mThick_Medialorbitofrontal_R_t2","mThick_Middletemporal_R_t2","mThick_Parahippocampal_R_t2","mThick_Paracentral_R_t2","mThick_Parsopercularis_R_t2","mThick_Parsorbitalis_R_t2","mThick_Parstriangularis_R_t2","mThick_Pericalcarine_R_t2","mThick_Postcentral_R_t2","mThick_Posteriorcingulate_R_t2","mThick_Precentral_R_t2","mThick_Precuneus_R_t2","mThick_Rostralanteriorcingulate_R_t2","mThick_Rostralmiddlefrontal_R_t2","mThick_Superiorfrontal_R_t2","mThick_Superiorparietal_R_t2","mThick_Superiortemporal_R_t2","mThick_Supramarginal_R_t2","mThick_Transversetemporal_R_t2","mThick_Insula_R_t2") 

# all Freesurfer DKT volumes, area, mThickness variables in UKBB

round1Mediators <- c(lobularBrainMetrics, generalMetrics)
round2Mediators <- c(aseg_brainVars, dkt_brainVars_area, dkt_brainVars_vol, dkt_brainVars_mThick)

mediators <- c(round1Mediators, round2Mediators)

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
# save(essentialVarsAnalysis, file = "revisedEssentialVarsAnalysis.RData")
eid <- which(colnames(df) == "eid")
missingEssVar <- which(colnames(df) == "missingEssentialVar")

df <- df[c(eid, essentialVarsAnalysis_colNum, missingEssVar)]
# colnames(df_small)

ordinalFactorContrasts <- list(exercise_IPAQActivityGroup_t0 = "contr.treatment")

### Correlation analyses ----
individualAnalysis <- list()
output_path <- "./outputs/cor/model" 
subsetDF_path <- subsetData_outputPath
#subsetDF_path <- ""
#subsetNames <- c("noMedNoDx")

dataDate <- getDate()
#dataDate <- ""
X <- "crp_log_z"
finalCovars_colnames <- colnames(df[finalCovars]) 
varsOfInterest <- c(X, finalCovars_colnames, cogOutcomes_colNames,cogOutcomes_fctr_names) 
listOfCovars <- c("demo_sex_t0Male","smoke_currently0_t0Yes","smoke_currently0_t0Yes","demo_ethnicity_t0White","exercise_IPAQActivityGroup_t0low","exercise_IPAQActivityGroup_t0moderate","ses_townsend_t0_z","demo_age_assess0_t0_z","demo_daysBtwAssess_z","weight_waistToHip_mean02_z","sleep_duration_mean02_z", "hand_t0LH","hand_t0Not","hand_t0RH", "brain_headScale_t2", "med_Antihypertensive_t021", "med_Statin_t021")
outputFileNames <- c()

for(subset in subsetNames){
  cat("\n", subset, ": performing correlation analysis and saving text files... \n")

  fileName <- glue("df_{subset}_{dataDate}.csv")
  subsetDf <- read.csv(glue("{subsetDF_path}/{fileName}"))
  subsetDf <- subsetDf[varsOfInterest]
  analysisName <- glue("cor_{subset}")
  
  # for numeric variables
  for(i in cogOutcomes_colNames){
    analysisName <- glue("{analysisName}_{i}")
    fileName <- glue("{analysisName}_{getDate()}")
    individualAnalysis <- runCor(df = subsetDf, X = X, Y = i, covars = colnames(df[finalCovars_colnames]), M = "", YisFactor = F, name = analysisName, contrasts = ordinalFactorContrasts)
    capture.output(summary(individualAnalysis), file = glue("{output_path}/{fileName}.txt"))
    cat("Correlation model for subset ", subset, "and cog var", i, "saved as ", glue("./outputs/cor/model/{fileName}")
    outputFileNames <- c(outputFileNames, fileName)
    analysisName <- glue("cor_{subset}_")
  }

  corTxtFormat(filePath = output_path, files = outputFileNames, X = X, listOfCovars = listOfCovars, outputPath = "./outputs/cor/summary", logReg = FALSE) # summarise outputs
  
  # for factor variables
  for(i in cogOutcomes_fctr_names){
    analysisName <- glue("{analysisName}_{i}")
    fileName <- glue("{analysisName}_{getDate()}")
    individualAnalysis <- runCor(df = subsetDf, X = "crp_log_z", Y = i, covars = colnames(df[finalCovars_colnames]), M = "", YisFactor = T, name = analysisName, contrasts = ordinalFactorContrasts)
    capture.output(summary(individualAnalysis), file = glue("{output_path}/{fileName}.txt"))
    cat("Correlation model for subset ", subset, "and cog var", i, "saved as ", glue("./outputs/cor/model/{fileName}")
    outputFileNames <- c(outputFileNames, fileName)
    # print(summary(listAnalyses_cor_noC))
    analysisName <- glue("cor_{subset}_")
  }

  rm(subsetDf)
  outputFileNames <- c()
}

### Mediation analyses -----
covarsMed <- c(finalCovars_colnames, covars_brain_colNames) # Add  cognition only covars
# colnames(df[covarsMed ])

individualAnalysis <- list()
listAnalyses_med_C <- list()
mainMed_C <- list()

# mediators <- c("area_Caudalanteriorcingulate_L_t2") # for testing purposes

subsetPath <- subsetData_outputPath
#subsetNames <- c("noMedNoDx")# list of subset names to run through
#dataDate <- "05_11_2022" # date of subsets, to recreate file name


alldfs <- c("df_noMedNoDx")
for(j in alldfs){
  subsetFileName <- glue("{j}_{dataDate}.csv")
  subsetName <- unlist(strsplit(j, "_"))[2]
  subsetDf <- read.csv(glue("{subsetPath}/{subsetFileName}"))
  
  df_mediatorsZ <- subsetDf[mediators] %>% 
    mutate(across(.fns = scale, .names = "{colnames(subsetDf[mediators])}_z")) # standardize brain morphology variables
  df_mediatorsZ <- df_mediatorsZ[,which(! colnames(df_mediatorsZ) %in% colnames(subsetDf[mediators]))] # matrix of z-scores made in above line to add to main df
  mediators_Z_colNames <- colnames(df_mediatorsZ)
  subsetDf <- cbind(subsetDf, df_mediatorsZ[mediators_Z_colNames])
  
  outputDF <- matrix(ncol = 40, dimnames = list(c(), c("modelName", "X", "M", "Y", "IndEff_est", "IndEff_SE", "IndEff_95%CI-Lo","IndEff_95%CI-Hi","dirEff_b", "dirEff_SE", "dirEff_X_t", "dirEff_X_p", "dirEff_R^2", "dirEff_R^2adj", "dirEff_F", "dirEff_df1", "dirEff_df2", "dirEff_p", "dirEff_AIC", "aPath_b", "aPath_SE", "aPath_X_t", "aPath_X_p", "aPath_model_R^2", "aPath_model_R^2adj", "aPath_model_F", "aPath_model_df1", "aPath_model_df2", "aPath_model_p", "bPath_b", "bPath_SE", "bPath_X_t", "bPath_X_p", "bPath_model_R^2", "bPath_model_R^2adj", "bPath_model_F", "bPath_model_df1", "bPath_model_df2", "bPath_model_p", "bPath_AIC")))
  
  analysisName <- glue("{subsetName}_med_C")
  
  for(med in colnames(df[mediators_Z_colNames])){
    for(i in cogOutcomes_colNames){
      analysisName <- glue("UKBB_{analysisName}_{med}_{i}")
      fileName <- runRMediation(df = subsetDf, X = "crp_log_z", Y = i, M = med, covars = colnames(df[covarsMed]), contrasts = ordinalFactorContrasts, YisFactor = F, name = analysisName, alpha = .05)
      analysisName <- glue("{subsetName}_med_C")
      outputDF <- formatMedOutput(outputDF = outputDF, fileName = fileName, filePath = getwd(), analysisName = analysisName,  logReg = F)
      fileToRemove <- as.character(glue("./{fileName}")) 
      file.remove(fileToRemove) # remove file
    }
    for(i in cogOutcomes_fctr_names){
      analysisName <- glue("UKBB_{analysisName}_{med}_{i}")
      fileName <- runRMediation(df = subsetDf, X = "crp_log_z", Y = i, M = med, covars = colnames(df[covarsMed]), contrasts = ordinalFactorContrasts, YisFactor = T, name = analysisName, alpha = .05)
      analysisName <- glue("{subsetName}_med_C")
      outputDF <- formatMedOutput(outputDF = outputDF, fileName = fileName, filePath = getwd(), analysisName = analysisName,  logReg = T)
      fileToRemove <- as.character(glue("./{fileName}")) 
      file.remove(fileToRemove) # remove file
    }
  }
  outputDF <- as.data.frame(outputDF)
  outputDF <- outputDF[-1,] # remove first row if as it is all NA values
  write.csv(outputDF, file = glue("{analysisName}_results_{getDate()}.csv"))
}

# Format outputs -------

listOfCovars <- c("demo_sex_t0Male","smoke_currently0_t0Yes","smoke_currently0_t0Yes","demo_ethnicity_t0White","exercise_IPAQActivityGroup_t0low","exercise_IPAQActivityGroup_t0moderate","ses_townsend_t0_z","demo_age_assess0_t0_z","demo_daysBtwAssess_z","weight_waistToHip_mean02_z","sleep_duration_mean02_z", "hand_t0LH","hand_t0Not","hand_t0RH", "brain_headScale_t2", "med_Antihypertensive_t021", "med_Statin_t021")

list_covarAllModels <- list() # this will store all dataframes with covariate information

filePath <- "./outputs/cor"
files <- c("UKBB_Anal_cor_cog_digsub_cor_t2_z_05_05_2022","UKBB_Anal_cor_cog_digsub_numAttempted_t2_z_05_05_2022","UKBB_Anal_cor_cog_fluidIntel_QsAttempted_t2_z_05_05_2022","UKBB_Anal_cor_cog_fluidIntel_score_t2_z_05_05_2022","UKBB_Anal_cor_cog_matrix_cor_t2_z_05_05_2022","UKBB_Anal_cor_cog_matrix_numViewed_t2_z_05_05_2022","UKBB_Anal_cor_cog_numMem_maxDigitRemem_t2_z_05_05_2022","UKBB_Anal_cor_cog_prospMem_timeDelay_t2_z_05_05_2022","UKBB_Anal_cor_cog_reactiontime_mean_t2_z_05_05_2022","UKBB_Anal_cor_cog_TMT_alphanumDuration_t2_z_05_05_2022","UKBB_Anal_cor_cog_TMT_alphanumErrors_t2_z_05_05_2022","UKBB_Anal_cor_cog_TMT_numericDuration_t2_z_05_05_2022","UKBB_Anal_cor_cog_TMT_numericErrors_t2_z_05_05_2022","UKBB_Anal_cor_cog_tower_cor_t2_z_05_05_2022")
catFile <- c("UKBB_Anal_cor_cog_prospMem_result_t2_05_05_2022") # cor analysis output of categorical variables

corFormat(depVar = "crp_log_z", modelName = "cor_C", outputPrefix = "UKBB_cor_anal_crpCog_C", filePath = "/home/doodlefish/Documents/Research/LepageLab/immunologyAndSz/Analysis/immunoCognition-new/Outputs/Correl-WithC", listOfCovars = listOfCovars, files = files)

# !! below function is broken !! -----
corFormat <- function(modelName, filePath, files, depVar, listOfCovars, outputPrefix){  
  df_main <- matrix(ncol = 13, dimnames = list(c(), c("model Name", "indep var", "dep var", "parameter estimate", "SE", "t", "p", "R^2", "R^2adj", "F", "df1", "df2", "p")))
  df_covarForModel <- matrix(ncol = length(listOfCovars)+2)
  df_covarForModel_tmp <- c()
  
  for(file in files){
    specificFilePath <- glue("{filePath}/{file}.txt")
    fileObj <- readLines(specificFilePath)
    
    lmExpression <- grep("~", fileObj, value = TRUE) # from this line, extract analysis name
    lmExpression <- unlist(strsplit(lmExpression, "\\s+"))
    indepVar <- lmExpression[which(lmExpression == "~")-1] # n.b. indep var of this model may not be indep var of mediation analysis
    name <- glue("{modelName}_{indepVar}")
    cat(name)
    depVarValues <- grep(depVar, fileObj, value = TRUE)[2]
    depVarValues <- unlist(strsplit(depVarValues, "\\s+"))
    
    if(TRUE %in% grepl("<", depVarValues)){
      pAndSymbol <- glue("{depVarValues[[5]]}{depVarValues[[6]]}")
      depVarValuesNoSigChar <- c(depVarValues[1:4], pAndSymbol)
    } else {
      depVarValuesNoSigChar <- depVarValues[1:5]
    }
    # cat(paste(name, depVarValuesNoSigChar, sep = "\n\t"))
    
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
    
    # if(listOfCovars[1] != ""){
    #   for(covar in listOfCovars){
    #     
    #     covarValues <- grep(covar, fileObj, value = TRUE)
    #     if(is.null(covarValues) == F){
    #       if(length(covarValues) == 2){
    #         covarValues <- covarValues[2]
    #       }
    #       indexNumBegin <- length(unlist(strsplit(covar, "\\s+")))
    #       covarValues <- unlist(strsplit(covarValues, "\\s+"))
    #   
    #       indexOfFirstNum <- indexNumBegin + 1
    #       indexOfLastNum <- indexNumBegin + 4
    #       
    #       if(TRUE %in% grepl("<", covarValues)){
    #         indexPreP <- indexOfLastNum - 1
    #         covarValues <- c(covar, covarValues[indexOfFirstNum : indexPreP], glue("{covarValues[[indexOfLastNum]]}{covarValues[[indexOfLastNum+1]]}"))
    #       # p <- paste(FTest[[9]], FTest[[10]], sep = "")
    #       } else {
    #         covarValues <- c(covar, covarValues[indexOfFirstNum:indexOfLastNum])
    #       }
    #     } else{
    #       covarValues <- rep("Error, missing data for covariates.", 5)
    #     }
    #     df_covarForModel_tmp <- cbind(df_covarForModel, covarValues)
    #   }
    #   name <- glue("{name}_C")
    #   if(sum(is.na(df_covarForModel_tmp[,1])) == nrow(df_covarForModel_tmp)){
    #      df_covarForModel_tmp <-  df_covarForModel_tmp[,-1]
    #   }
    #   dimnames(df_covarForModel) <- list(c(), listOfCovars)
    #   df_covarForModel <- df_covarForModel_tmp
    #   df_covarForModel <- cbind("modelName" = name, "statistic" = c("parameter estimate", "SE", "t", "p"), df_covarForModel_tmp)
    #   if(nrow(df_covarForModel) == 1){
    #     df_covarForModel <- rbind(df_covarForModel[-1,], df_covarForModel_tmp)
    #   } else {
    #     df_covarForModel <- rbind(df_covarForModel, df_covarForModel_tmp)
    #   }
    #   df_covarForModel <- c()
    #   write.csv(df_covarForModel, file = glue("{outputPrefix}_CovarWeights_{getDate()}.csv"))
    # } else {
    #   name <- glue("{name}_noC")
    # }
    
    df_mainRow <- c(name, indepVar, depVarValuesNoSigChar, R2, R2adj, FStat, df1, df2, p)  # "Model name", "indep var", "dep var", "parameter estimate", "SE", "t", "p", "R^2", "R^2adj", "F", "df1", "df2", "p"
    df_main <- rbind(df_main, " " = df_mainRow)
  }
  
  df_main <- df_main[-1,]
  write.csv(df_main, file = glue("{outputPrefix}_CRPWeights_{getDate()}.csv"))
} # 'modelName' - name of model for labelling output rows (usually 'cor' or 'med'), 'filePath' - path to folder housing files, 'files' - name of files (without extensions), 'depVar' - dependent variable (usually CRP_log_z), 'listOfCovars' - list of covariate variable names as listed in LM output (include level of contrast for categorical variables). If no covars, define as '""', 'outputPrefix' - string indicating the prefix for the output file. Note, independent variable is assumed to be the variable beginning with "cog_"

# Correl analyses
# listOfFiles <- c("UKBB_Anal_cor_revC_cog_digsub_cor_t2_z_03_26_2022","UKBB_Anal_cor_revC_cog_fluidIntel_QsAttempted_t2_z_03_26_2022") # a list of all file names

# Correlational analyses -----
listOfCovars <- c("demo_sex_t0Male","smoke_currently0_t0Yes","smoke_currently0_t0Yes","demo_ethnicity_t0White","exercise_IPAQActivityGroup_t0low","exercise_IPAQActivityGroup_t0moderate","ses_townsend_t0_z","demo_age_assess0_t0_z","demo_daysBtwAssess_z","weight_waistToHip_mean02_z","sleep_duration_mean02_z")
dataDate <- getDate()
corOutput_path <- "./outputs/cor" 

fileNames <-list()
for(subset in subsetNames){
  temp <- glue("{subset}_{dataDate}.txt")
  fileNames <- c(fileNames,temp)
}

corFormat(depVar = "crp_log_z", modelName = "cor_", outputPrefix = "cor_", filePath = corOutput_path, listOfCovars = listOfCovars, files = fileNames)

# Mediation analyses
## !!!! NEED UPDATING !!!! -----
listOfMediators <- c("brain_vol_brainSegNoVent_t2_z","brain_vol_hippocamp_L_t2_z","brain_vol_hippocamp_R_t2_z") # list of mediator variable names
df_covarForModel <- matrix(nrow = 5)
df_covarForModel <- matrix(ncol = length(listOfCovars)+2)
list_covarAllModels <- list() # this will store all dataframes with covariate information
# initialize counters
file <- ""
covar <- ""

filePath <- "/home/doodlefish/Documents/Research/LepageLab/immunologyAndSz/Analysis/immunoCognition-new/Outputs/MedAnalyses/C/revC_04_01_2022_95CI"
files <- c("UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_digsub_cor_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_fluidIntel_QsAttempted_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_prospMem_timeDelay_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_fluidIntel_score_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_reactiontime_mean_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_matrix_cor_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_TMT_alphanumDuration_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_numMem_maxDigitRemem_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_TMT_alphanumErrors_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_PC_1_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_TMT_numericDuration_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_PC_2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_TMT_numericErrors_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_PC_3_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_tower_cor_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_digsub_cor_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_prospMem_timeDelay_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_fluidIntel_QsAttempted_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_reactiontime_mean_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_fluidIntel_score_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_TMT_alphanumDuration_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_matrix_cor_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_TMT_alphanumErrors_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_numMem_maxDigitRemem_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_TMT_numericDuration_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_PC_1_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_TMT_numericErrors_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_PC_2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_brainSegNoVent_t2_z_cog_tower_cor_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_PC_3_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_digsub_cor_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_fluidIntel_QsAttempted_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_prospMem_timeDelay_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_fluidIntel_score_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_reactiontime_mean_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_matrix_cor_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_TMT_alphanumDuration_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_numMem_maxDigitRemem_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_TMT_alphanumErrors_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_PC_1_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_TMT_numericDuration_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_PC_2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_TMT_numericErrors_t2_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_L_t2_z_cog_PC_3_z_04_01_2022","UKBB_Anal_med_revC_brain_vol_hippocamp_R_t2_z_cog_tower_cor_t2_z_04_01_2022")

listOfCovars <- ""

# filePath = "/home/doodlefish/Documents/Research/LepageLab/immunologyAndSz/Analysis/immunoCognition-new/Outputs/MedAnalyses/C"
# files = "UKBB_Anal_med_C_brain_vol_brainSegNoVent_t2_z_cog_PC_3_z_03_26_2022"
# covars <- c("brain_vol_brainSegNoVent_t2_z","demo_sex_t0Male","demo_ethnicity_t0East Asian","demo_ethnicity_t0Mixed","demo_ethnicity_t0Other ethnic group", "demo_ethnicity_t0South Asian", "demo_ethnicity_t0White","demo_age_assess0_t0","demo_daysBtwAssess_z", "hand_t0LH","hand_t0Not","hand_t0RH","brain_headScale_t2","cog_hourCompleted")

variables <- c()
values <- c()
valueList <- c()
df_main <- matrix(ncol = 38, dimnames = list(c(), c("modelName", "X", "M", "Y", "IndEff_est", "IndEff_SE", "IndEff_95%CI-Lo","IndEff_95%CI-Hi","DirEff_b", "DirEff_SE", "DirEff_X_t", "DirEff_X_p", "DirEff_model_R^2", "DirEff_model_R^2adj", "DirEff_model_F", "DirEff_model_df1", "DirEff_model_df2", "DirEff_model_p", "aPath_b", "aPath_SE", "aPath_X_t", "aPath_X_p", "aPath_model_R^2", "aPath_model_R^2adj", "aPath_model_F", "aPath_model_df1", "aPath_model_df2", "aPath_model_p", "bPath_b", "bPath_SE", "bPath_X_t", "bPath_X_p", "bPath_model_R^2", "bPath_model_R^2adj", "bPath_model_F", "bPath_model_df1", "bPath_model_df2", "bPath_model_p")))
name <- "med_revC"
bPathRows <- c()

for(file in files){
  specificFilePath <- glue("{filePath}/{file}.txt")
  fileObj <- readLines(specificFilePath)
  
  splits <- grep("~", fileObj)
  split1 <- splits[2]-2
  split2 <- splits[3]-2
  medOutput <- grep("Mediated Effect", fileObj)
  
  text_directEffect <- fileObj[2:split1]
  text_aPath <- fileObj[(split1+1):split2]
  text_bPath <- fileObj[(split2+1):(medOutput-1)]
  text_mediatedEffect <- fileObj[medOutput:length(fileObj)]
  models <- list(text_directEffect, text_aPath, text_bPath)
  
  for(i in 1:(length(models))){ # iterate through all three models
    # Find name of mediator, indep var
    sectionText <- unlist(models[i])
    lmExpression <- grep("~", sectionText, value = TRUE) # from this line, extract analysis name
    lmExpression <- unlist(strsplit(lmExpression, "\\s+"))
    indepVar <- lmExpression[which(lmExpression == "~")-1] # n.b. indep var of this model may not be indep var of mediation analysis
    depVar <- lmExpression[which(lmExpression == "~")+1]
    
    depVarValues <- grep(depVar, sectionText, value = TRUE)
    if(length(depVarValues) == 2){
      depVarValues <- depVarValues[2]
    }
    
    depVarValues <- unlist(strsplit(depVarValues, "\\s+"))
    
    if(TRUE %in% grepl("<", depVarValues) && nchar(depVarValues[grep("<",depVarValues)]) == 1){
      pAndSymbol <- glue("{depVarValues[[5]]}{depVarValues[[6]]}")
      depVarValuesNoSigChar <- c(depVarValues[2:4], pAndSymbol)
    } else {
      depVarValuesNoSigChar <- depVarValues[2:5]
    }
    
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
      p <- FTest[which(FTest == grep("p-value",FTest, value = TRUE))+1]
    }
    
    ########################################################33
    # bPathRows <- grep("b path", sectionText, value = TRUE)
    # if(is.null(bPathRows) == FALSE){ # extract med -> Y in b path model
    #   brainVarValues <- sectionText[starts_with("brain_", vars = sectionText)]
    #   for(i in length(brainVarValues)){
    #     splitRow <- strsplit(brainVarValues[i], "\\s+")
    #     for(j in splitRow){
    #       if(splitRow %in% listOfMediators == T){
    #         print("hi")
    #       }
    #     }
    #   }
    #   
    #   contains(brainVarValues, vars = listOfMediators)
    #   brainVarValues <- unlist(strsplit(brainVarValues, "\\s+"))
    #   if(TRUE %in% grepl("<", brainVarValues) && nchar(brainVarValues[grep("<", brainVarValues)]) == 1){
    #     pAndSymbol <- glue("{brainVarValues[[5]]}{brainVarValues[[6]]}")
    #     brainVarValuesNoSigChar <- c(brainVarValues[2:4], pAndSymbol)
    #   } else {
    #     brainVarValuesNoSigChar <- brainVarValues[2:5]
    #   }
    #   df_covarForModel <- cbind(df_covarForModel,  brainVarValuesNoSigChar)
    # }
    #####################################################################333
    
    # listOfCovars <- c(, )
    # if(listOfCovars[1] != ""){
    #   for(covar in listOfCovars){
    #   #   covarValues <- grep(covar, sectionText, value = TRUE)
    #   #   if(is.null(covarValues) == F){
    #   #     if(length(covarValues) == 2){
    #   #       covarValues <- covarValues[2]
    #   #     }
    #   #     indexNumBegin <- length(unlist(strsplit(covar, "\\s+")))
    #   #     covarValues <- unlist(strsplit(covarValues, "\\s+"))
    #   # 
    #   #     indexOfFirstNum <- indexNumBegin + 1
    #   #     indexOfLastNum <- indexNumBegin + 4
    #   #     
    #   #     if(TRUE %in% grepl("<", covarValues)){
    #   #       indexPreP <- indexOfLastNum - 1
    #   #       covarValues <- c(covar, covarValues[indexOfFirstNum : indexPreP], glue("{covarValues[[indexOfLastNum]]}{covarValues[[indexOfLastNum+1]]}"))
    #   #     # p <- paste(FTest[[9]], FTest[[10]], sep = "")
    #   #     } else {
    #   #       covarValues <- c(covar, covarValues[indexOfFirstNum:indexOfLastNum])
    #   #     }
    #   #   } else{
    #   #     covarValues <- rep("Error, missing data.", 5)
    #   #   }
    #   #   df_covarForModel <- cbind(df_covarForModel, covarValues)
    #   # }
    #   # name <- glue("{name}_C")
    #   # if(sum(is.na(df_covarForModel[,1])) == nrow(df_covarForModel)){
    #   #     df_covarForModel <-  df_covarForModel[,-1]
    #   # }
    #   # dimnames(df_covarForModel) <- list(c(), listOfCovars)
    #   # df_covarForModel_tmp <- df_covarForModel[-1,]
    #   # df_covarForModel_tmp <- cbind("modelName" = name, "statistic" = c("parameter estimate", "SE", "t", "p"), df_covarForModel_tmp)
    #   # if(nrow(df_covarForModel) == 1){
    #   #   df_covarForModel <- rbind(df_covarForModel[-1,], df_covarForModel_tmp)
    #   # } else {
    #   #   df_covarForModel <- rbind(df_covarForModel, df_covarForModel_tmp)
    #   # }
    #   # df_covarForModel <- c()
    #   # write.csv(df_covarForModel, file = glue("{outputPrefix}_CovarWeights_{getDate()}.csv"))
    #   }
    # } else {
    #   # name <- glue("{name}_noC")
    # }
    variables <- append(variables, c(indepVar, depVar))
    values <- append(values, c(depVarValuesNoSigChar, R2, R2adj, FStat, df1, df2, p))
  }
  
  indepVar <- variables[2]
  medVar <- variables[3]
  depVar <- variables[1]
  vars <- c(indepVar, medVar, depVar)
  name <- glue("{name}_{medVar}_{depVar}")
  
  # Extract indirect effect
  sectionText <- text_mediatedEffect
  
  medEff_est <- grep("Estimate",sectionText, value = TRUE)
  medEff_est <- unlist(strsplit(medEff_est, "\\s+"))
  medEff_est <-  medEff_est[3]
  
  medEff_SE <- grep("SE",sectionText, value = TRUE)
  medEff_SE <- unlist(strsplit(medEff_SE, "\\s+"))
  medEff_SE <-  medEff_SE[3]
  
  medEff_95Low <- grep("95%CILow",sectionText, value = TRUE)
  medEff_95Low <- unlist(strsplit(medEff_95Low, "\\s+"))
  medEff_95Low <-  medEff_95Low[3]
  
  medEff_95Hi<- grep("95%CIHi",sectionText, value = TRUE)
  medEff_95Hi <- unlist(strsplit(medEff_95Hi, "\\s+"))
  medEff_95Hi <-  medEff_95Hi[3]
  
  indEff <- c(medEff_est, medEff_SE, medEff_95Low, medEff_95Hi)
  
  df_main <- rbind(df_main, c(name, vars, indEff, values)) # specify the indices of valueList to retain in main_df
  
  variables <- c()
  values <- c()
  name <- "med_noC"
}
as.data.frame(df_main)
df_main <- df_main[-1,]
write.csv(df_main, file = glue("{outputPrefix}_CRPWeights_{getDate()}.csv"))

# 3 - Compute p-Values and FDR correction ----
path <- "./outputs/med/summary" # path to raw mediation output summary files
fileNames <- c("UKBB_med_C_results_CRPWeights_05_05_2022","UKBB_excludingDx_med_C_results_05_06_2022","UKBB_excludingSSRI_med_C_results_05_06_2022","UKBB_oldestTert_med_C_results_05_06_2022","UKBB_onlyDx_med_C_results_05_06_2022","UKBB_onlyMedUsers_med_C_results_05_06_2022","UKBB_youngestTert_med_C_results_05_06_2022", "UKBB_NoSSRINoDx_med_C_results_05_09_2022", "UKBB_NoMedNoDx_med_C_results_05_08_2022", "UKBB_noMed_med_C_results_05_08_2022") # name of summary mediation outcome files
fileNames <- c("UKBB_noDxNoSSRI_med_C_results_05_11_2022")

## Correlation analyses --------------
### goals: Determine significant CRP - cog var associations from mediation analysis output
### a) Make df with results from direct effect of CRP on cognition vars
### b) FDR correction on these associations

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

## Mediation analysis --------------
### Goals: 
### a) compute p-values for indirect effects
### b) FDR correction for indirect effects oof cognition variables that are retained

keyCognition <- c("cog_fluidIntel_score_t2_z", "cog_numMem_maxDigitRemem_t2_z") # name of cognition variables of interest for mediation analyses

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
    outputDF_keyCogOnly_TrendSig <- outputDF_keyCogOnly %>% filter(IndEff_p_cor < .05)
    
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
  write.csv(outputDf_p_cor, file = saveName)
  outputDf_p_cor <- c()
} # summarize significant mediators

# View(outputDf_p_cor)

for(file in fileNames){
  fileName <- glue("{path}/{file}.csv")
  outputDF <- read_csv(fileName, show_col_types = F)
  variableNameCols <- c("X", "Y", "M")
  colOfInt_names <- c("IndEff_est","IndEff_SE","IndEff_p","IndEff_p_cor","IndEff_95%CI-Lo","IndEff_95%CI-Hi")
  colOfInt_nums <- which(colnames(outputDF) %in% c(colOfInt_names,variableNameCols))
  trimmedDF_filtered <- trimmedDF %>% filter(IndEff_p_cor < .1)
  
  # trimmedDF <- trimmedDF[order(trimmedDF$dirEff_X_p),]
  # p_cor <- p.adjust(trimmedDF$dirEff_X_p, method = "fdr")
  # trimmedDF <- trimmedDF %>%   
  #   mutate("X_p_cor" = p_cor, .after = "dirEff_X_p")
  # 
  # saveName <- glue("{file}_cogSummary.csv")
  # write.csv(trimmedDF, file = saveName)
  # cat("\n Summary of associations between crp_log_z and cognition variables with corrected p-values saved as: ", saveName)
}

# 4 - Visualisation ----
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
                c_prime_eff <- format(round(unlist(df_relevant[row, colNum_totEff_b]) - unlist(df_relevant[row, colNum_indEff_b]), nsmall = 3)) # direct effect computed by taking difference between totEff and indEff.
                                      
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
infix <- "_med_C_pCor" # include any seperators
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
