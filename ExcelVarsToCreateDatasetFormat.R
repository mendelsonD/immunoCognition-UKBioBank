# Author: Daniel Mendelson,
# Date: Feb 5, 2022

# Goal of this script: Convert variables of interest from an excel spread sheet to text formatted according to the requirements of create_dataset's config file.
# For details about the create_dataset script, see https://github.com/joshunrau/CognitiveSubtypes/tree/main/src/create_dataset .

# Input and output details:
  # Input: csv file with n rows (one per variable) and the columns:
    #   DataField: [int] datafield from the UKBB
    #   VarName: [str] desired name for this variable in final subset
    #   Coding: [str] 
      # If the string contains the substring ' = ', then there must be values and interpretations in the form:
      # ['value1 = Interpretation1; value2 = Interpretation2; ' ... 'valueN = InterpretationN']
      # N.B. each code value must be seperate by ' = ' from its meaning. Each value/meaning pair must be seperated by '; '
      # If the string does contain ' = ', then there is assumed to be no coding.
    #   InstanceNum: [int]
    #   NumArray: [int] max position of the array of interest (for variables not in array format, specify 0)
    #   Include: [boolena: 'T', 'True' or 'F', 'False'] specifies whether this variable should be included in subset.

  # output: text file with list for Josh's create_dataset's config.py file.
    #   'VarName': {
    #     'DataField': {DataField},
    #     'InstanceNum': {InstanceNum},
    #     'ArrayRange': range(0, {NumArray}),
    #     'Included': {Include:'True' or 'False'},
    #     'Coding': {
    #         "{codeValue1-string}": "{Interpretation1_String}",
    #         "{codeValueM-string}": "{InterpretationM_String}",
    #       } or ['None,']
    #   },

library(glue)
library(readxl)
library(stringr)
library(dplyr)


df <- read_xlsx("/home/doodlefish/Documents/Research/LepageLab/immunologyAndSz/Analysis/immunoCognition-new/data/Variables/VariableLists/UKBB_BrainVariables-05_02_2022.xlsx", sheet = 1)
# df <- read_xlsx("./UKBB_Variables-03_12_2022.xlsx", sheet = 1)

head(df)
dim(df)
getwd()
setwd("./Outputs/") # set working directory to desired output location
outputName <- str_c("Daniel-BrainVarList_Create-Dataset-", format(Sys.time(), "%m"), "_", format(Sys.time(), "%d"), "_", format(Sys.time(), "%Y"), ".txt")

df <- df %>%
  mutate_all(as.character)
# df <- df[-which(df$VarName == "eid"),]
output <- ""

for(i in 1:nrow(df)){
  VarName <- df[i, which(colnames(df) == "VarName")]
  DataField <- df[i, which(colnames(df) == "DataField")]
  Instance <- df[i, which(colnames(df) == "InstanceNum")]
  
  if(df[i, which(colnames(df) == "NumArray")] == 0){
    MaxArrayRange <- 1
  } else {
    MaxArrayRange <- df[i, which(colnames(df) == "NumArray")]
  }
  
  if(grepl(pattern = c("T|true|yes|y|1"), x = df[i, which(colnames(df) == "Include")], ignore.case = T)  == T){ 
    Include <- "True" 
  } else if(grepl(pattern = c("F|false|no|n|0"), x = df[i,which(colnames(df) == "Include")], ignore.case = T)  == T){
    Include <- "False" 
  } else { 
    cat(glue("Error. Variable, ", df[i,which(colnames(df) == "VarName")], " has an incorrect value for the 'Include' column. Accepted values are 'T' or 'F'. It is assumed to be False."))
    Include <- "False"
  }
  
  # coding
  if(grepl(" = ", df[i,which(colnames(df) == "Coding")], fixed = TRUE) == T){ # check if contains " = "
    Coding <- "{"
    individualValueMeaningPair <- strsplit(df[i,which(colnames(df) == "Coding")][[1]], "; ")[[1]] # separate each code/variable pair, split based on '; '
    for(j in 1:length(individualValueMeaningPair)){
      ValMean <- strsplit(individualValueMeaningPair[[j]][1], split = " = ")
      if(is.na(ValMean[[1]][2]) == F){ # If there is data specifying units, this argument will be ignored
        meaning <- ValMean[[1]][2] # extract meaning
        value <- ValMean[[1]][1] # extract code 
        Codej <- glue("\t\t\t \"{value}\": \"{meaning}\",")
        Coding <- paste(Coding, "\n\t", Codej)
      }
    }
    Coding <- paste(Coding, "\n\t\t\t}")
  } else {
    Coding <- "None,"
  }
  
  outputI <- glue("'{VarName}': {{ 
                    \t 'DataField': {DataField},
                    \t 'InstanceNum': {Instance},
                    \t 'ArrayRange': range(0, {MaxArrayRange}),
                    \t 'Included': {Include},
                    \t 'Coding': {Coding}
                  }},\n")
  
  output <- append(output, c(outputI, "\n")) # add current variable details to complete output
}

capture.output(writeLines(output), file = outputName) # save output to file
cat("The output file `", outputName, "` has been saved to `", getwd(), "`.") # print message specifying the new file's name and location.
