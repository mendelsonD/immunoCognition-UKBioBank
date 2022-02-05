# Daniel Mendelson, Feb 4, 2022

# Goal of this script: Convert variables of interest from an excel spread sheet to Josh's config file format.
# Input details: csv file with n rows (one per variable) and the columns:
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

# output details: text file with list for Josh's create_dataset's config.py file.
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

df <- read_xlsx("/home/doodlefish/Documents/Research/LepageLab/immunologyAndSz/Analysis/immunoCognition/UKBB_Variables-DanielVars.xlsx")
df <- df %>%
  mutate_all(as.character)
# View(df)
# nrow(df)
setwd("./") # set working directory to desired output location
outputName <- str_c("Daniel-VarList_Create-Dataset-", format(Sys.time(), "%m"), "_", format(Sys.time(), "%d"), "_", format(Sys.time(), "%Y"), ".txt")
output <- ""

for(i in 1:nrow(df)){
  VarName <- df[i, "VarName"]
  DataField <- df[i, "DataField"]
  Instance <- df[i, "InstanceNum"]
  
  if(df[i,"NumArray"] != "0"){
    MaxArrayRange <- df[i, "NumArray"]
  } else {
    MaxArrayRange <- 1
  }
  
  if(grepl(pattern = c("T|true|yes|y|1"), x = df[i,"Include"], ignore.case = T)  == T){ 
    Include <- "True" 
  } else if(grepl(pattern = c("F|false|no|n|0"), x = df[i,"Include"], ignore.case = T)  == T){
    Include <- "False" 
  } else { 
    cat(glue("Error. Variable, ", df[i,"VarName"], " has an incorrect value for the 'Include' column. Accepted values are 'T' or 'F'. It is assumed to be False."))
    Include <- "False"
  }
  
  # coding
  if(grepl(" = ", df[i,"Coding"], fixed = TRUE) == T){ # check if contains " = "
    Coding <- "{"
    individualValueMeaningPair <- strsplit(df[i,"Coding"][[1]], "; ")[[1]] # separate each code/variable pair, split based on '; '
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

# cat(output)

{ 
  capture.output(output, file = outputName, row.names = F) 
  cat("The output file ", outputName, "has been saved to: ", getwd(), ".") 
} # save output to file and print message specifying the new file's name and location.
