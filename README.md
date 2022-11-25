# immunoCognition-UKBioBank
Repository for analyses on UK Biobank data assessing C-Reactive Protein's relationship to cognitive performance. Also assesses metrics of structural brain imaging as potential mediators in that relationship.


Flow of scripts:
1 - DataPreperation: import data, prepare data for analyses, make subsets, perform descriptive statistics;
2 - Analyses: perform correlational and mediation analyses on data subsetted and prepared in 1;
3 - Format outputs: summarise output files from correlational analyses;
4 - Interpret results: perform FDR-correction;
5 - Visualise: Make figures.

Additional tools:
ExcelVarsToCreateDatasetFormat: takes list of UK Biobank variables and their details, prepares a text file for use as a library in the 'createDataset' program file to access UKBiobank data. The createDataset program was prepared by Mr. Joshua Unrau (see https://github.com/mendelsonD/CognitiveSubtypes/tree/main/src/create_dataset).
