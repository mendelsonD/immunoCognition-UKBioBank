#####
# This file contains is the analysis script for the immunology, cognition and 
# brain mediation analysis of the UKBB.
# Prepared by Daniel Mendelson with help from Katie M. Lavigne and Joshua Unrau.
# Work performed under the supervision of Dr. Martin Lepage.
# 19 May 2023
# RStudio v2023.03.1 Build 446
# R v4.3.0
#####

# Preamble -----
packages <- c("DescTools", "DiagrammeR", "DiagrammeRsvg", "dplyr", "effsize", 
              "factoextra", "ggplot2", "ggseg", "ggthemes", "glue", "lavaan",
              "lubridate", "psych", "rcompanion", "readr", "rsvg", "stats", 
              "stringr", "tidyr")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(packages)

setwd("C:/Users/katie/OneDrive - McGill University/publishing/_rev_Mendelson_UKBB/res")

# 0 - Functions
## Misc functions
getDate <- function(){
  require(stringr)
  return(str_c(format(Sys.time(), "%m"), "_", format(Sys.time(), "%d"), "_", format(Sys.time(), "%Y")))
}
## Assumption check functions ----
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

## Compare non-normal variables
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
  write.csv(wilcoxOutput, file = glue("{outputDir}/{fileName}_{dataDate}.csv")) #exports wilcoxOutput to a CSV file
  print(paste("The Wilcoxon comparison has been saved as:", glue("{fileName}_{dataDate}.csv")))
} # "varsOfInterest" - list of variable names to compare between groups; "df" data frame with all vars of interest and grouping variable; "by" - grouping variable; "outputDir" - path to output directory; "fileName" - desired name of output file

## Compare kept cases to excluded cases 
compareExcluded <- function(df_subset, varsToCompare, grp, subset, simplified){
  require(DescTools)
  require(effsize)
  date <- dataDate # set todays date for easier output filenaming
  
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

## Analysis Functions

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
  
  fileName <- glue("{name}_{dataDate}.txt")
  
  capture.output(
    cat(" -- Linear Model Analysis Output -- "),
    print(c("\n ------------------ \n Model expression: \n\t", expression)),
    print(summary(model)),
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
  write.csv(df_main, file = glue("{outputPath}/{outputPrefix}_{dataDate}.csv"))
} # 'modelName' - name of model for labeling output rows (usually 'cor' or 'med'), 'corResultPath' - path to folder housing correlation text files, 'dataFilePrefix' - prefix of correlation files, 'subsets' - list of subset names, 'cogVars' - list of cognitive variables correlation analyses performed on, 'depVar' - dependent variable (usually CRP_log_z), 'dataDate' - date of correlational analysis (as in file name), 'listOfCovars' - list of covariate variable names as listed in LM output (include level of contrast for categorical variables). If no covars, define as '""', 'outputPrefix' - string indicating the prefix for the output file. Note, independent variable is assumed to be the variable beginning with "cog_", 'outputPath' - string defining path to save summary file.

### Mediation functions ----
runMediation <- function(df, X, M, Y, covars){
  require(lavaan)
  cov = paste(covars, collapse=' + ')
  
  # 'M ~ a*X + cov
  # Y ~ b*M + c*X + cov
  # direct := c
  # indirect := a*b
  # total := c + (a*b)
  # perc := (indirect / total)*100'
  
  model <- sprintf('%s ~ a*%s + %s \n %s ~ b*%s + c*%s + %s \n direct := c \n indirect := a*b \n total := c + (a*b) \n perc := (indirect / total)*100', M, X, cov, Y, M, X, cov)
  
  fit <- sem(model, data=df, se='bootstrap', bootstrap=500)
  return(fit)
}

# 1 - UK BioBank data analysis ----
dataDate <- getDate()

## Data Preparation -----
df <- read_csv("../data/RawData/Daniel_2022-05-02.csv") # import df
eidToRemove <- read_csv("../data/RawData/w45551_2023-04-25.csv", col_names = "eid") # file with EID of participants who withdrew consent from UKBB

df <- df %>% filter(!(eid %in% eidToRemove)) # remove rows with eid in this file
write.csv(df, "../data/RawData/Daniel_{dataDate}_ProperEID.csv", row.names=FALSE)
rm(eidToRemove) # remove dataframe(s) not in use

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

# Num days between assessment 2 and assessment 0
df <- df %>% 
  mutate(demo_daysBtwAssess = as.numeric(df$date_assess2_t2 - df$date_assess0_t0, units = "days"))

df <- df %>% 
  mutate(crp_hourCollected = hour(df$crp_timeCollected_t0)) %>% 
  mutate(cog_hourCompleted = hour(df$cog_timeCompleted_t2)) %>% 
  mutate(brain_hourCompleted = hour(df$brain_timeCompleted_t2))

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

# Medication - any of interest
medVars <- starts_with("med_", vars = colnames(df))
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

### Make columns appropriate data types ----

df <- readr::type_convert(df) # automatically detect
factorVars <- which(colnames(df) %in% c("demo_sex_t0", "demo_ethnicity_t0", "cog_prospMem_result_t2", "hand_t0","smoke_currently0_t0", "smoke_currently2_t2", "menopause0_t0", "menopause2_t2", "anyMedOfInterest0", "anyMedOfInterest2", "anyDxOfInterest", colnames(df[dxVars]), colnames(df[medVars0]), colnames(df[medVars2]))) # Factors
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
vars_BrainRegion <- read_csv("./Variables/brainVarsbyLobe.csv")
vars_BrainRegion$varNameDf <- glue("{vars_BrainRegion$VariableName}_t2")
cat("\n Create lobar brain variable columns...")

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
saveRDS(df, file = glue("./Data/Processed/df_full_{dataDate}.rds"))

### Load cleaned df ----
df <- readRDS(glue("./Data/Processed/df_full_{dataDate}.rds"))
### Exclusions ----
demoVars <- starts_with("demo_", vars = colnames(df))
sesVars <- starts_with("ses_", vars = colnames(df))
crpAliqVars <- which(colnames(df) %in% c("crp_aliquot_t0", "crp_log", "crp_aliquot_fctr"))
crpToZ <- which(colnames(df) %in% c("crp_aliquot_t0", "crp_log"))

brainPrefix <- c("area_","gmVol_", "GWcont_","mIntensity","mThick_", "vol_", "weigted_mFA_", "weigted_mICVF_", "weigted_mISOVF_", "weigted_mL1_", "weigted_mL2_", "weigted_mL3_", "weigted_mMD_", "weigted_mMO_", "weigted_mOD_")
brainMorphVars <- c()
for(i in brainPrefix){
  iVars <- starts_with(i, vars = colnames(df))
  brainMorphVars <- c(brainMorphVars, iVars)
}

dietVars <- starts_with("diet_", vars = colnames(df))
dietVarsToCor <- which(colnames(df) %in% c("diet_cookedVeg0_t0", "diet_cookedVeg2_t2", "diet_rawVeg_t0", "diet_fruit0_t0", "diet_fruit2_t2", "diet_processedMeat0_t0", "diet_processedMeat2_t2", "diet_water_t0"))

smokingVars <- starts_with("smoke_", vars = colnames(df))
weightVars <- c(starts_with("weight_", vars = colnames(df)), which(colnames(df) == "waistToHip"))

exerciseVars <- starts_with("exercise_", vars = colnames(df))

completionTimeVars <- c(contains("timeCompleted", vars = colnames(df)), contains("timeCollected", vars = colnames(df)))

cogVars <- starts_with("cog", vars = colnames(df))
cogVars_t2 <- cogVars[ends_with("_t2", vars = colnames(df[cogVars]))]
cogVarsToZ <- cogVars_t2[which(!colnames(df[c(cogVars_t2)]) %in% c("cog_prospMem_result_t2", "cog_timeCompleted_t2", "cog_hourCompleted"))]

covarsToZ <- which(colnames(df) %in% c("demo_age_assess0_t0", "demo_age_assess2_t2", "age_mean02", "demo_daysBtwAssess", "ses_townsend_t0", "weight_waistToHip_t0", "weight_waistToHip_t2", "weight_waistToHip_mean02", "sleep_duration0_t0", "sleep_duration2_t2", "sleep_duration_mean02"))
varsToZ <- c(crpToZ, cogVarsToZ, covarsToZ)

df_varsToZ <- df[varsToZ] %>% 
  mutate(across(.fns = scale, .names = "{colnames(df[varsToZ])}_z"))
vars_Z <- which(! colnames(df_varsToZ) %in% colnames(df[varsToZ])) # matrix of z-scores made in above line to add to main df
df <- cbind(df, df_varsToZ[vars_Z])

df$cog_prospMem_result_t2 <- as.factor(df$cog_prospMem_result_t2)

finalCovars <- which(colnames(df) %in% c("demo_sex_t0", "demo_ethnicity_t0", "demo_age_assess0_t0_z", "demo_daysBtwAssess_z", "ses_townsend_t0_z", "weight_waistToHip_mean02_z", "sleep_duration_mean02_z", "exercise_IPAQActivityGroup_t0", "med_Antihypertensive_t02", "med_Statin_t02"))

### Remove missing ----
## Essential variables
essentialVars <- c(finalCovars, completionTimeVars, crpAliqVars, brainMorphVars) # list of variable names that must be complete in order for case to be retained
df <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    if_any(.cols = essentialVars, function(x)(
      is.na(x))) ~ TRUE,
    TRUE ~ FALSE))) # mark missing essential to remove

## CRP > 10
cat(sum(df$crp_aliquot_t0 > 10, na.rm = T), " cases to remove given CRP > 10.", sep = "")
df <- df %>%
  mutate(highCRP = if_else(df$crp_aliquot_t0 > 10, TRUE, FALSE, missing = FALSE)) # mark CRP > 10 to remove

saveRDS(df, file = glue("./Data/Processed/df_PreSubset_{dataDate}.rds"))

df_removed <- df %>%
  filter(if_any(any_of(essentialVars), function(x)(is.na(x)))) # make new df of cases with missing values

numMissing <- matrix(ncol=2, nrow = length(essentialVars), dimnames = list(colnames(df_removed[essentialVars]), c("NACount", "NotMissingOnThisVar")))

for(i in 1:length(essentialVars)){
  numMissing[i,] <- c(sum(is.na(df_removed[essentialVars[i]])), sum(!is.na(df_removed[essentialVars[i]])))
}

numMissing_df <- as.data.frame(numMissing) %>% arrange(desc(NACount))
write.csv(numMissing_df, file = paste("UKBB_sourceOfNAValues_", dataDate, ".csv", sep = ""))

# Appendix C. # missing info
missing_vars <- c("cog_TMT_alphanumDuration_t2", "cog_tower_cor_t2", "cog_TMT_numericDuration_t2", "cog_matrix_cor_t2", "cog_digsub_cor_t2", "cog_TMT_numericErrors_t2", "cog_TMT_alphanumErrors_t2", "cog_numMem_maxDigitRemem_t2", "cog_timeCompleted_t2", "exercise_IPAQActivityGroup_t0", "crp_aliquot_t0", "cog_reactiontime_mean_t2", "crp_aliquotdelay", "vol_EstimatedTotalIntraCranial_WB_t2", "crp_timeCollected_t0", "demo_ethnicity_t0", "ses_townsend_t0", "crp_timeFasting_t0")
missing_varsCols = c()
for(var in missing_vars){
  missing_varsCols <- c(missing_varsCols, which(colnames(df) == var))
}
na_count <-sapply(df[missing_varsCols], function(y) sum(length(which(is.na(y)))))
write.csv(na_count, file="appendix-C.csv")

### Load filtered df ----
df <- readRDS(glue("./Data/Processed/df_PreSubset_{dataDate}.rds"))
### Assumption checks ----
## Normality
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
subsetData_outputPath <- "./Data/Processed/Subsets"
subsetNames <- c("noMedNoDx", "all", "oldest","onlyDx","noDx","onlyMed","noMed","onlySSRI","youngest","noDxNoSSRI","noSSRI")

### Subset 1: no further exclusions (only missing data and CRP > 10) ----
comparisonExcluded_all <- compareExcluded(df, varsToCompare, grp = "missingEssentialVar", subset = "all", simplified = T) 
#runWilcox(df = df, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "all_wilcox")

## Filter data
df_all <- df %>% 
  filter(missingEssentialVar == FALSE) %>%
  filter(highCRP == FALSE)

fileName <- paste("df_all_", dataDate, sep = "")
write.csv(df_all, file = glue("{subsetData_outputPath}/{fileName}.csv")) # exports cleaned dataframe with all variables ready for analysis
cat("Subset only with complete cases and CRP < 10 saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_all))

### Subset 2: oldest tertile ----
tertCutOff <- quantile(df$demo_age_assess0_t0, c(.66))
df_oldest <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    demo_age_assess0_t0 <= tertCutOff ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_oldest <- compareExcluded(df_oldest, varsToCompare, grp = "missingEssentialVar", subset = "oldestTert", simplified = T) 
#runWilcox(df = df_oldest, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "oldest_wilcox")

df_oldest <- df_oldest %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_oldest_", dataDate, sep = "")
write.csv(df_oldest, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset with the oldest tertile of cases at instance 0 (", tertCutOff, ") saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_oldest))

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
fileName <- paste("df_youngest_", dataDate, sep = "")
write.csv(df_youngest, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset with the youngest tertile of cases at instance 0 (", tertCutOff, ") saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_youngest))

### Subset 4: Exclude Diagnoses -------
df_noDx <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    dx_AnyOfInterest == 0 ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_noDx <- compareExcluded(df_noDx, varsToCompare, grp = "missingEssentialVar", subset = "NoDx", simplified = T)
#runWilcox(df = df_noDx, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "noDx_wilcox")

df_noDx <- df_noDx %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_noDx_", dataDate, sep = "")
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
fileName <- paste("df_onlyDx_", dataDate, sep = "")
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
fileName <- paste("df_noMed_", dataDate, sep = "")
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
fileName <- paste("df_onlyMed_", dataDate, sep = "")
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
fileName <- paste("df_noSSRI_", dataDate, sep = "")
write.csv(df_noSSRI, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset excluding those using SSRIs saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_noSSRI))

### Subset 9: Only SSRI ----
df_onlySSRI <- df %>% 
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    med_SSRI_t0 == 1 & med_SSRI_t2 == 1 ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_onlySSRI <- compareExcluded(df_onlySSRI, varsToCompare, grp = "missingEssentialVar", subset = "onlySSRI", simplified = T)
#runWilcox(df = df_onlySSRI, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "onlySSRI_wilcox")

df_onlySSRI <- df_onlySSRI %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_onlySSRI_", dataDate, sep = "")
write.csv(df_noSSRI, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset only with those using SSRIs saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_onlySSRI))

### Subset 10: exclude Med and Dx --------
df_noMedNoDx <- df %>%
  mutate(missingEssentialVar = factor(case_when(
    missingEssentialVar == TRUE ~ TRUE,
    dx_AnyOfInterest == 0 & med_AnyOfInterest0 == 0 & med_AnyOfInterest2 == 0 ~ FALSE,
    TRUE ~ TRUE)))

comparisonExcluded_noMedNoDx <- compareExcluded(df_noMedNoDx, varsToCompare, grp = "missingEssentialVar", subset = "NoMedNoDx", simplified = T)
#runWilcox(df = df_noMedNoDx, varsOfInterest = , by = "missingEssentialVar", outputDir = "./outputs/descriptive/wilcox", fileName = "noMedNoDx_wilcox")

df_noMedNoDx <- df_noMedNoDx %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_noMedNoDx_", dataDate, sep = "")
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

df_noDxNoSSRI <- df_noDxNoSSRI %>% filter(missingEssentialVar == FALSE)
fileName <- paste("df_noDxNoSSRI_", dataDate, sep = "")
write.csv(df_noDxNoSSRI, file = glue("{subsetData_outputPath}/{fileName}.csv"))
cat("Subset excluding those with diagnoses or using SSRIs saved as: `", fileName, "`. \n\t complete cases: ", nrow(df_noDxNoSSRI))

### Describe Subsets ----
path <- subsetData_outputPath
subsetNames <- c("noDxNoSSRI", "noMedNoDx", "all", "oldest","onlyDx","noDx","onlyMed","noMed","onlySSRI","youngest","noSSRI")

for(subset in subsetNames){
  fileName <- glue("{path}/df_{subset}_{dataDate}.csv")
  subset_df <- read_csv(fileName, lazy = T, show_col_types = F)
  print(subset)
  write.csv(psych::describe(subset_df), file = glue("./outputs/descriptive/{subset}_summary_{dataDate}.csv"))
  rm(subset_df)
}

## Summarise
write.csv(df_all, file = paste("./data/processed/dfRetained_", dataDate, ".csv", sep = ""))
df <- df_all
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

write.csv(summary(df[-numVars]), file = paste("./outputs/descriptive/df_all_summary_fctrVars", dataDate, ".csv", sep = "")) # save counts of categorical variables


## Variable correlations----
### Correl of CRP vars 
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
## Run analyses -----
### Import data -----
path <- subsetData_outputPath
fileName <- "df_noDxNoSSRI"
df <- read.csv(glue("./{path}/{fileName}_{dataDate}.csv")) # import df
### Vars of interest ----
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

covars_colNames <- c("demo_sex_t0", "demo_ethnicity_t0", "demo_age_assess0_t0_z", "demo_daysBtwAssess_z", "ses_townsend_t0_z", "weight_waistToHip_mean02_z", "sleep_duration_mean02_z", "smoke_currently0_t0", "smoke_currently2_t2","exercise_IPAQActivityGroup_t0", "med_Antihypertensive_t02", "med_Statin_t02")
covars_colNum <-  which(colnames(df) %in% covars_colNames)

covars_brain_colNames <- c("hand_t0", "brain_headScale_t2")
covars_brain_colNum <- which(colnames(df) %in% covars_brain_colNames)

essentialVarsAnalysis_colNum <- c(crp, cogOutcomes_colNum, cogOutcomes_fctr, brainMorphVars_Z_colNum, covars_colNum, covars_brain_colNum) # list of variables used in correlation/mediation analyses
essentialVarsAnalysis_colNames <- colnames(df[essentialVarsAnalysis_colNum]) 
eid <- which(colnames(df) == "eid")
missingEssVar <- which(colnames(df) == "missingEssentialVar")

df <- df[c(eid, essentialVarsAnalysis_colNum, missingEssVar)]

ordinalFactorContrasts <- list(exercise_IPAQActivityGroup_t0 = "contr.treatment")

### Correlation analyses ----
individualAnalysis <- list()
output_path <- "./outputs/cor/model" 
subsetDF_path <- subsetData_outputPath

X <- "crp_log_z"
varsOfInterest <- c(X, covars_colNames, cogOutcomes_colNames,cogOutcomes_fctr_names) 
listOfCovars <- c("demo_sex_t0Male","smoke_currently0_t0Yes","smoke_currently2_t2Yes","demo_ethnicity_t0White","exercise_IPAQActivityGroup_t0low","exercise_IPAQActivityGroup_t0moderate","ses_townsend_t0_z","demo_age_assess0_t0_z","demo_daysBtwAssess_z","weight_waistToHip_mean02_z","sleep_duration_mean02_z", "hand_t0LH","hand_t0Not","hand_t0RH", "brain_headScale_t2", "med_Antihypertensive_t021", "med_Statin_t021")
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
    fileName <- glue("{analysisName}_{dataDate}")
    individualAnalysis <- runCor(df = subsetDf, X = X, Y = i, covars = covars_colNames, M = "", YisFactor = F, name = analysisName, contrasts = ordinalFactorContrasts)
    capture.output(summary(individualAnalysis), file = glue("{output_path}/{fileName}.txt"))
    #cat("Correlation model for subset ", subset, "and cog var", i, "saved as ", glue("./outputs/cor/model/{fileName}"))
    outputFileNames <- c(outputFileNames, fileName)
    analysisName <- glue("cor_{subset}_")
  }

  corTxtFormat(filePath = output_path, files = outputFileNames, X = X, listOfCovars = listOfCovars, outputPath = "./outputs/cor/summary", logReg = FALSE) # summarise outputs
  
  # for factor variables
  for(i in cogOutcomes_fctr_names){
    analysisName <- glue("{analysisName}_{i}")
    fileName <- glue("{analysisName}_{dataDate}")
    individualAnalysis <- runCor(df = subsetDf, X = "crp_log_z", Y = i, covars = covars_colNames, M = "", YisFactor = T, name = analysisName, contrasts = ordinalFactorContrasts)
    capture.output(summary(individualAnalysis), file = glue("{output_path}/{fileName}.txt"))
    cat("Correlation model for subset ", subset, "and cog var", i, "saved as ", glue("./outputs/cor/model/{fileName}"))
    outputFileNames <- c(outputFileNames, fileName)
    # print(summary(listAnalyses_cor_noC))
    analysisName <- glue("cor_{subset}_")
  }

  rm(subsetDf)
  outputFileNames <- c()
}

## Mediation ----
Ys <- c("cog_fluidIntel_score_t2_z", "cog_numMem_maxDigitRemem_t2_z")
X <- c("crp_log_z")
Ms <- mediators
covarsMed <- c(covars_colNames, covars_brain_colNames) # Add  cognition only covars
subsetNames <- c("noMedNoDx", "all", "oldest","onlyDx","noDx","onlyMed","noMed","onlySSRI","youngest","noSSRI")
#subsetNames <- c("noDxNoSSRI","noMedNoDx", "all", "oldest","onlyDx","noDx","onlyMed","noMed","onlySSRI","youngest","noSSRI")
#subsetNames <- c("noDxNoSSRI")

for(subset in subsetNames){
  df <- read.csv(glue("./{path}/df_{subset}_{dataDate}.csv")) # import df
  outdir <- glue('./outputs/med/subsets/{subset}')
  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  if(subset=="noMedNoDx" | subset=="noMed"){
    covarsMed <- covarsMed[covarsMed !="med_Antihypertensive_t02"]
    covarsMed <- covarsMed[covarsMed !="med_Statin_t02"]
  }
  for(Y in 1:length(Ys)){
    for(M in 1:length(Ms)){
      dfMed <- cbind(df[X], scale(df[Ms[M]]), df[Ys[Y]], df[covarsMed]) # Reduce df & standardize brain vars
      dfMed <- na.omit(dfMed) # Filter out missing data
      
      fit <- runMediation(dfMed, X, Ms[M], Ys[Y], covarsMed)
      #summary(fit, fit.measures=T, rsq=T, standardized = T, ci = T)
      # output fit indices
      medInfo <- data.frame("Subset" = subset, "X" = X, "M" = Ms[M], "Y" = Ys[Y])
      medFit <- as.data.frame(t(as.data.frame(fitmeasures(fit))))
      medFitOut <- cbind(medInfo, medFit)
      rownames(medFitOut) = NULL
      if (Y==1 && M==1) {
        headerTF=TRUE
      } else {
        headerTF=FALSE
      }
      write.table(medFitOut, glue('{outdir}/{subset}_X-{X}_M-{Ms[M]}_Y-{Ys[Y]}_med-fit-output.csv'),
                  append = TRUE,
                  sep = ",",
                  col.names = headerTF,
                  row.names = FALSE,
                  quote = FALSE)
      # output parameters
      medParam <- parameterestimates(fit, rsquare=TRUE)
      medParam <- medParam[(medParam$label!="" | medParam$op=="r2"),]
      rows <- nrow(medParam)
      medInfoParam <- data.frame("Subset" = subset, "X" = rep_len(X,rows), "M" = rep_len(Ms[M], rows), "Y" = rep_len(Ys[Y], rows))
      t <- cbind(medInfoParam,medParam)
      write.table(t, glue('{outdir}/{subset}_X-{X}_M-{Ms[M]}_Y-{Ys[Y]}_med-param-output.csv'),
                  append = TRUE,
                  sep = ",",
                  col.names = headerTF,
                  row.names = FALSE,
                  quote = FALSE)
    }
  }
  # Combine CSVs into master CSVs
  ## Fit indices
  file_names <- list.files(path = outdir,
                           pattern = "*med-fit-output.csv", full.names = TRUE, recursive = TRUE)
  fit_all <- do.call(rbind, lapply(file_names, read.csv, header = FALSE))
  colnames(fit_all) <- colnames(medFitOut)
  write.csv(fit_all,glue('{outdir}/{subset}_med-fit-all-output.csv'))
  
  ## Parameters
  file_names <- list.files(path = outdir,
                           pattern = "*med-param-output.csv", full.names = TRUE, recursive = TRUE)
  param_all <- do.call(rbind, lapply(file_names, read.csv, header = FALSE))
  colnames(param_all) <- colnames(t)
  
  ### Correct pvalues
  labels <- c("a", "b", "c", "direct", "indirect", "total")
  for(l in labels){
    tmp <- param_all %>% 
      filter(label == l)
    p_cor <- p.adjust(tmp$pvalue, method = "fdr")
    tmp <- tmp %>%   
      mutate("pvalue_cor" = p_cor, .after = "pvalue")
    if(l == "a"){
      param_cor <- tmp
    } else{
      param_cor <- rbind(param_cor, tmp)
    }
  }
  write.csv(param_cor,glue('{outdir}/{subset}_med-param-all-output-corrected.csv'))
}

### Combine all subsets
for(subset in subsetNames){
  ## Fit indices
  file_names <- list.files(path = './outputs/med/subsets/',
                           pattern = "*med-fit-output.csv", full.names = TRUE, recursive = TRUE)
  fit_all <- do.call(rbind, lapply(file_names, read.csv, header = FALSE))
  colnames(fit_all) <- colnames(medFitOut)
  write.csv(fit_all,glue('./outputs/med/med-fit-all-output.csv'))
  
  ## Parameters
  file_names <- list.files(path = './outputs/med/subsets/',
                           pattern = "*med-param-output.csv", full.names = TRUE, recursive = TRUE)
  param_all <- do.call(rbind, lapply(file_names, read.csv, header = FALSE))
  colnames(param_all) <- colnames(t)
  write.csv(param_all,glue('./outputs/med/med-param-all-output.csv'))
}

### Re-adjust p values by round
param_all <- read.csv(glue("./outputs/med/med-param-all-output.csv"))
param_all <- param_all %>%   
  mutate("round" = case_when(M %in% round1Mediators ~ "round1",
                             M %in% round2Mediators ~ "round2"),
         .after = "Subset")

param_all <- param_all %>% 
  group_by(Subset, round, label) %>%
  mutate("pvalue_cor" = case_when(label != "" ~ p.adjust(pvalue, method = "fdr")),
           .after = "pvalue")
write.csv(param_all,glue('./outputs/med/med-param-all-output-corrected.csv'))

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
round1_MedDiagrams <- function(df, dataDate, subsets, rounds, mediatorsToDiagram, cogsToDiagram, outputPath){
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
  
  for(subset_i in subsets){
    cat("Making mediation figures for: ", subset_i, "\n")
    for(round_j in rounds){
      for(meds in mediatorsToDiagram){ # For each mediator
        metric <- str_split(meds, "_")[[1]][1]
        roi <- str_split(meds, "_")[[1]][2]
        hemi <- str_split(meds, "_")[[1]][3]
        medName <- glue("{metric}_{roi}_{hemi}") # for output file name
        cat(medName, "\n")
        dfMedName <- paste(glue("{medName}_t2")) # for finding the unique mediator in df
        for(cogName in cogsToDiagram){
          if(dfMedName %in% df$M) { # Check that this combination exists. If it doesn't skip
            df_relevant <- df %>% filter(Subset == subset_i, round == round_j, M == dfMedName, Y == cogName)
            
            a_eff <- format(round(df_relevant[df_relevant$label == "a", ]$est, 3), nsmall = 3) # 'aPath_b'
            a_p <- pToStar(df_relevant[df_relevant$label == "a", ]$pvalue_cor) # 'aPath_model_p'
            
            b_eff <- format(round(unlist(df_relevant[df_relevant$label == "b", ]$est), 3), nsmall = 3) # 'bPath_b'
            b_p <- pToStar(df_relevant[df_relevant$label == "b", ]$pvalue_cor) # 'bPath_model_p'
            
            c_eff <- format(round(unlist(df_relevant[df_relevant$label == "total", ]$est), 3), nsmall = 3) # 'totEff_b'
            c_p <- pToStar(df_relevant[df_relevant$label == "total", ]$pvalue_cor) # 'totEff_p'
            
            c_prime_eff <- format(round(unlist(df_relevant[df_relevant$label == "direct", ]$est), 3), nsmall = 3) # direct effect computed by taking difference between totEff and indEff
            c_prime_p <- pToStar(df_relevant[df_relevant$label == "direct", ]$pvalue_cor) # 'direct effect p'
            
            ab_indEff <- format(round(unlist(df_relevant[df_relevant$label == "indirect", ]$est), 3), nsmall = 3) # colName: "IndEff_est_correct"
            ab_p <- pToStar(df_relevant[df_relevant$label == "indirect", ]$pvalue_cor) # colname: "IndEff_p_cor" or "IndEff_p" (defined above) 
            
            #### Make labels. ----
            ##### Line labels -----
            a <- glue("a = {a_eff}{a_p}") # Text for A path in form: "{ind eff estimate (0.###)} {significance ('*', '**', '***' or '')}"
            b <- glue("b = {b_eff}{b_p}") # Text for B path in form: "{ind eff estimate (0.###)} {significance ('*', '**', '***' or '')}"
            c <- glue("c = {c_eff}{c_p} \n ab = {ab_indEff}{ab_p} \n c\` = {c_prime_eff}{c_prime_p}") # Text for C, indirect effect (ab), and c prime in form: "c = {ind eff estimate (0.###)} {significance ('*', '**', '***' or '')} \n ab = {ind eff estimate (0.###)} {significance ('*', '**', '***' or '')} \n c\' = {dir eff estimate (0.###)}"
            
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
            if(cogName == "cog_fluidIntel_score_t2_z"){
              lab_y <- "Fluid Intelligence Score"
              cog <- "FluidIntelScore"
            } else if(cogName == "cog_numMem_maxDigitRemem_t2_z"){
              lab_y <- "Num. Mem. (Max dig.)"
              cog <- "NumMemMaxDig"
            } else {
              lab_y <- cogName
              cog <- cogName
            }
            
            ## Make and save diagram. ----
            outputName <- glue("{outputPath}/medDia_{subset_i}_{cog}_{medName}_{dataDate}.png") # name of diagram when saved
            
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
          } else {
            cat("Warning. \`", medName, "\` is an invalid combination; it is not a mediator in at least one of the provided files. \n This combination will be skipped without making a mediation diagram. \n")
          }
        }
      }
    }
  }
}# Automatically make and save diagrams for specifies brain regions and metrics. 'dataPath' string of path to folder containing output of mediation analyses !!with corrected p-values and ideally the combined MASTER document; 'dataFileName' string of name of df for specific data subset containing mediation analyses with corrected p-value for indirect effect. Df must contain the column names: 'M', 'Y', 'IndEff_est', ('IndEff_p_cor' or 'IndEff_p'), 'aPath_b', 'aPath_X_p' 'bPath_b', 'bPath_X_p'; 'dataFileInfix': infix of data file name, include any seperators between the file name and the infix; 'dataDate': date of data file; 'rois' list of brain regions to make mediation models for; 'hemispheres' hemispheres of interest. Can be 'L', 'R' or 'WB'; 'metrics' list of metrics of interest. Can be 'area', 'mThick', and/or 'vol'; 'pCor' boolean specifying if corrected p-value of ind effect is to be used; 'outputPath' string specifying where to save output diagrams

### 4.A.ii. Make diagrams -----

df <- read.csv(glue("./outputs/med/med-param-all-output-corrected.csv"))
outputPath <- "./outputs/figures/Figure1_medDia"

## Make and save diagram. ----
subsetsToDiagram <- c("noDxNoSSRI") # list of subset names to make diagrams for
roundsToDiagram <- c("round1") # list of rounds to diagram
mediatorsToDiagram <- round1Mediators
cogsToDiagram <- c("cog_fluidIntel_score_t2_z", "cog_numMem_maxDigitRemem_t2_z")

round1_MedDiagrams(df = df, 
                   dataDate = dataDate,
                   subsets = subsetsToDiagram,
                   rounds = roundsToDiagram,
                   mediatorsToDiagram = mediatorsToDiagram,
                   cogsToDiagram = cogsToDiagram, 
                   outputPath = outputPath)

## 4.B Figure 2: ggSegmentation ----
### 4.B.i. Functions ----
createBrainFigure <- function(df_subset, savePath, dataDate){
  require(stringr)
  require(ggplot2)
  require(ggseg)
  require(ggthemes)
  require(tidyr)
  
  Y <- df_subset$Y[1]
  subset <- df_subset$Subset[1]
  
  if(Y == "cog_fluidIntel_score_t2_z"){
    Yname_title <- "Fluid Intelligence - Score"
  } else if(Y == "cog_numMem_maxDigitRemem_t2_z") {
    Yname_title <- "Numeric Memory - Max. digits remembered"
  } else {
    return(cat("Error: ", Y, " is invalid. \n Options are 'cog_fluidIntel_score_t2_z', or 'cog_numMem_maxDigitRemem_t2_z'"))
  }
  
  figureName <- glue("{subset}_{Y}_indirect_est")
  
  valueName <- "Indirect effect parameter estimate"
  value <- glue("{Y}.IndEff_est")
  figureName <- glue("{subset}_{Yname_title}_indEff_est")
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
  
  cat("\n --------- \n Value:", value, "\n")
  
  regionColNum <- which(colnames(df_subset) == "Region")
  
  # I. Extract DKT brain regions ---------
  ## I.a. makes a list of all ROIs to keep for DKT figures. ----
  DKTRegions <- brain_labels(dk) # name of all DKT labels
  DKT_matches <- c()
  for(roi in DKTRegions){
    new_matches <- df_subset[[regionColNum]][grepl(roi, df_subset[[regionColNum]], ignore.case = TRUE, fixed = FALSE) == TRUE]
    DKT_matches <- c(DKT_matches, unique(new_matches))
  } 
  DKT_matches <- unique(DKT_matches)
  
  # II. Extract rows with value of column 'Region' specified in matches. ----
  reallyShort_DK <- subset(df_subset, df_subset[[regionColNum]] %in% DKT_matches) # Extracts all columns.
  
  # III. Make and format df with only necessary variables. ----
  valueColNum <- which(colnames(reallyShort_DK) == "est")
  metricColNum <- which(colnames(reallyShort_DK) == "Metric")
  
  reallyShort_DK <- reallyShort_DK[,c(regionColNum, metricColNum, valueColNum)] 
  colnames(reallyShort_DK) <- c("label", "measure", "dfValue")
  reallyShort_DK$label <- as.factor(reallyShort_DK$label)

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
  
  (cat("Saving DKT figure: ", glue("{figureName}_DK_{dataDate}.jpg"), "\n") + ggsave(glue("{savePath}/{figureName}_DK_{dataDate}.jpg"), dpi=450))
} # Creates ggSegmentation for DSK atlas. Takes: 'file' a data file formatted appropriately (see ROIS.csv for reference); 'figureValue' string specifying the values the figure is to represent (either 'p', 'indEff', or 'pCor'); 'savePath' path where figure should be saved; 'Y' dependent variable name

#### Call functions -----
df <- read.csv(glue("./outputs/med/med-param-all-output-corrected.csv"))
outputPath <- "./outputs/figures/Figure2_ggSeg"

subsetNames <- c("noDxNoSSRI")
rounds <- c("round2")
cogNames <- c("cog_fluidIntel_score_t2_z", "cog_numMem_maxDigitRemem_t2_z")

for(subset in subsetNames){
  cat("\n---- ", subset, "-----------------------\n")
  for(cogVar in cogNames){
    # Remove irrelevant rows
    df_subset <- df %>% filter(Subset == subset, round == rounds, Y == cogVar, label == "indirect")
    # Reformat mediator variables
    df_subset <- df_subset %>%
      tidyr::separate(col = 'M', sep = "_", into = c("Metric", "Region", "Hand", "t2")) %>% 
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
      dplyr::select(-c(t2, Hand))
    # Make figures
    createBrainFigure(df_subset, outputPath, dataDate)
  }
} # generate figures for each subset and for each cognitive variable