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
#       } or ['None']
#   },

library(glue)
library(readxl)

df <- read.csv("UKBB_Variables.csv")

for(i in 1:length(df)){
  VarName <- df[i, Varname]
  DataField <- df[i, Datafield]
  Instance <- df[i, InstanceNum]
  
  if(df[i,ArrayRange] != 0){
    MaxArrayRange <- df[i,ArrayRange]
  } else {
    MaxArrayRange <- 1
  }
  
  if(df[i,Include] == "T" | df[i,Include] == "True" | df[i,Include] == "Y" | df[i,Include] == "Yes"){ 
    Include <- "True" 
  } else if(df[i,Include] == "F" | df[i,Include] == "False" | df[i,Include] == "N" | df[i,Include] == "No"){
    Include <- "False" 
  } else { 
    cat("Error. Variable, ", df[i,VarName], "has an incorrect value for the 'Include' column. Accepted values are 'T' or 'F'. It is assumed to be False.") 
    Include <- "False"
  }
  
  # coding
  if(grepl(" = ", df[i,coding], fixed = TRUE) == T | grepl("= ", df[i,coding], fixed = TRUE) == T | grepl(" =", df[i,coding], fixed = TRUE) == T){ # check if contains " = "
    Coding <- "{"
    individualValueMeaningPair <- split() # seperate each code/variable pair, split based on '; '
    for(j in individualValueMeaningPair){
      value <- # extract code 
      meaning <- # extract meaning  
      Codej <- glue('\n\t "{value}": "{meaning}",')
      Coding <- append(Coding, Codej)
    }
    Coding <- append(Coding, "\n}")
  } else {
    Coding <- "None"
  }
  
  outputI <- glue("'{VarName}': { \n 
                    \t 'DataField' : {DataField}, \n
                    \t 'InstanceNum': 0, \n
                    \t 'ArrayRange': range(0, {MaxArrayRange}), \n
                    \t 'Included': True, \n
                    \t 'Coding': Coding \n
                  },")
  
  output <- append(output, outputI)
}

setwd() # set working directory to desired output location
write.table(output, "VarList_Create-Dataset.txt")