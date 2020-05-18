#import libraries
library("BiocManager")
library("BiocVersion")
library("flowCore")
library("openCyto")
library("flowViz")
library("flowWorkspace")
library("ncdfFlow")
library(ggplot2)
library(ggcyto)
library(tibble)

#Source the analyse that creates plots, extracts RFP and GFP values and calculates RFP/GFP ratio for a specified dataset
source("functions/analyse_function.R")
#Source function that combines the tables for all mutants measured at the same time point
source("functions/makeTables_function.R")
#Source function that make a barplot of the mean RFP, GFP and ratio per mutants
source("functions/makeBarplot_function.R")

#Read in a list of fcs files fom the data diectory with the name structure 'YYYYMMDD_XXX_Zh.fcs'
file_list<-list.files("data", full.names = TRUE)

#Apply the 'analyse' function to all files in the list
lapply(file_list, analyse_function)

#For data with just one circumstance
makeTables_onetime_function()

#For data that has measurements on different timepoints (can be adjusted if the data has other different circumstances)
    #Make a list of all time points collected
    timepoints <- c("0h", "1h", "3h")

    #Apply the "makeTables" function to all time points
    lapply(timepoints,makeTables_function)

#For data with just one circumstance
makeBarplot_onetime_function()
    
#For data that has measurements on different time pints    
  #Apply the "makeHistograms" function to each time point
  lapply(timepoints, makeBarplot_function)



