library(tibble)
library(ggplot2)
#Create a function that will read in the combined table and create 
#3 barplots showing the mean RFP, GFP and ratio, respectively
makeBarplot_onetime_function <- function() {
        #read in file
        d<-read.csv(paste0("output_tables/.csv"), sep =" ")
        
        #create pdf file for RFP
        png(paste0("output_plots/bar_RFP.png"), width = 250, height = 250)
        
                #plot RFP, make bar colours red
                barplot(d$RFP, names.arg=rownames(d), col = "red", ylim = c(0,30000), ylab="RFP intensity", main = "RFP")
                dev.off()
        
        #create pdf file for GFP
        png(paste0("output_plots/bar_GFP.png"), width = 250, height = 250)
        
                #plot GFP, make bar colours green
                barplot(d$GFP, names.arg=rownames(d), col = "green", ylim = c(0,20000), ylab="GFP intensity", main = "GFP")
                dev.off()
        
        #create pdf file for Ratio
        png(paste0("output_plots/bar_Ratio.png"), width = 250, height = 250)
        
                #plot Ratios, make bar colours blue
                barplot(d$Ratio, names.arg=rownames(d), col = "blue", ylim = c(0,3), ylab="RFP/GFP", main = "Ratio")
                dev.off()
}

#Create a function that for a specified time point (inductiontime) will read in the combined table and create 
#3 barplots showing the mean RFP, GFP and ratio, respectively, for each mutant
makeBarplot_function <- function(inductiontime) {
  #read in file
  d<-read.csv(paste0("output_tables/",inductiontime,".csv"), sep =" ")
  
  #create pdf file
  png(paste0("output_plots/bar_RFP_",inductiontime,".png"), width = 250, height = 250)
  #plot RFP, make bar colours red
  barplot(d$RFP, names.arg=rownames(d), col = "red", ylim = c(0,30000), ylab="RFP intensity", main = "RFP")
  dev.off()
  
  #create pdf file for GFP
  png(paste0("output_plots/bar_GFP_",inductiontime,".png"), width = 250, height = 250)
  #plot GFP, make bar colours green
  barplot(d$GFP, names.arg=rownames(d), col = "green", ylim = c(0,20000), ylab="GFP intensity", main = "GFP")
  dev.off()
  #create pdf file for Ratio
  png(paste0("output_plots/bar_Ratio_",inductiontime,".png"), width = 250, height = 250)
  #plot Ratios, make bar colours blue
  barplot(d$Ratio, names.arg=rownames(d), col = "blue", ylim = c(0,3), ylab="RFP/GFP", main = "Ratio")
  dev.off()
}



