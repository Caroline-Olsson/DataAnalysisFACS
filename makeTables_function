library(tibble)
#Make a function that will combine all the output tables created by the analyse function 
makeTables_onetime_function <- function() {
        #read in a list of the files containing the tables for the specified timepoint
        tables<-list.files("output_tables/*.csv", full.names = TRUE)
        
        #create empty data frame
        df <-data.frame()
        for (i in 1:length(tables)) {
                #read in each table and transpose
                data <- t(read.csv(tables[i], header = FALSE))
                #add transposed table to the empty data frame
                df <- rbind(df, data.frame("Mutation" = data[1], "RFP" = data[2], "GFP" = data[3], "Ratio" = data[4]))
        }
        #change the first column to rownames
        table <- column_to_rownames(df, var = "Mutation")
        
        #Write table
        write.table(table, file=((paste0("output_tables/table.csv"))) )
}

#Make a function that will combine all the output tables created by the analyse function 
#for a specified time point (inductiontime)
makeTables_function <- function(inductiontime) {
  #read in a list of the files containing the tables for the specified timepoint
  tables<-list.files("output_tables", pattern = paste0("*",inductiontime,".csv"), full.names = TRUE)
  
  #create empty data frame
  df <-data.frame()
  for (i in 1:length(tables)) {
    #read in each table and transpose
    data <- t(read.csv(tables[i], header = FALSE))
    #add transposed table to the empty data frame
    df <- rbind(df, data.frame("Mutation" = data[1], "RFP" = data[2], "GFP" = data[3], "Ratio" = data[4]))
  }
  #change the first column to rownames
  table <- column_to_rownames(df, var = "Mutation")
  
  #Write table
  write.table(table, file=((paste0("output_tables/",inductiontime,".csv"))) )
}
