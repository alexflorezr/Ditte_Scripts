#====== Create input files using species database =====
rm(list=ls())
# Go to the directory with the databases
setwd("/Users/afr/Desktop/Temporal_test/")


#### TO DELETE BEFORE READ THE FUNCTION ########
ABC_table_species <-  "species_constant.txt"     ###
Database <- "DATABASE.txt"                   ###
s <- species_vector[1]
s <- "Saiga_tatarica"
################################################


ABC_infile_single <- function(ABC_table_species, Database){
  #### Only one population will be defined
  #### Define the variables for the function
  ABC_table <- read.delim(ABC_table_species ,header=T, sep="\t", stringsAsFactors=FALSE, na.strings="NA")
  species_vector <- as.vector(ABC_table$Species)
  for (s in seq_along(species_vector)){
    for (i in seq_along(colnames(ABC_table))){
      assign(colnames(ABC_table)[i], ABC_table[s,i])
    }
    # Create a conection (e.i. an empty file) to write the text
    fileConn <-file(paste(species_vector[s], "single", paste(Event_type,".par", sep=""), sep="-"), "w")
    # The input files contain 12 blocks, see http://web.stanford.edu/group/hadlylab/ssc/index.html
    # Before every block there is a description comment, which start with //
    ####################################################################################################
    ########                1st block, used to define the number of populations                 ########
    ####################################################################################################
    writeLines(paste("//", species_vector[s], "input parameters for BayesSSC"), fileConn)
    writeLines(paste(Npop, "population(s) with ancient data"), fileConn)
    ####################################################################################################
    ########                  2nd block, used to define the demes (populations) sizes           ########
    ####################################################################################################
    writeLines(paste("//", "Deme sizes (haploid number of genes)"), fileConn)
    writeLines(as.character(Spop), fileConn)
    ####################################################################################################
    ########                  3rd block, used to define the sample sizes                        ########
    ####################################################################################################
    writeLines(paste("//", "Sample sizes: # samples, age in generations, deme, stat group"), fileConn)
    # Read the database
    ABC_DB <- as.data.frame(read.delim(Database,header=T, sep="\t", stringsAsFactors=FALSE, na.strings="NA"))
    # Create a temporal table with the frequencies per record for each species
    temp_single <- as.data.frame(table(ABC_DB$Median_Age[which(ABC_DB$Species == species_vector[s] & nchar(ABC_DB$Sequence) > 2)]), stringsAsFactors=F)
    # Create an empty matrix to save the values for the sampling groups
    ABC_sampling <- matrix(data=as.numeric(), nrow=dim(temp_single)[1], ncol=4)
    # Send the sequence sampling to the first column
    ABC_sampling[,1] <- as.numeric(temp_single$Freq)
    # Send the sequence dating to the second column
    ABC_sampling[,2] <- round(as.numeric(temp_single$Var1)/Gen_time)
    # Only one population 0
    ABC_sampling[,3] <- 0
    # Statistical groups defined by a climatic event
    for ( i in seq_along(ABC_sampling[,4])){
      if (ABC_sampling[i,2]*Gen_time < Climatic_time){
        ABC_sampling[i,4] <- 0
      }
      if (ABC_sampling[i,2]*Gen_time > Climatic_time){
        ABC_sampling[i,4] <- 1
      }
    }
    writeLines(paste(length(ABC_sampling[which(ABC_sampling[,1] > 0)]), "sample groups"), fileConn)
    write.table(ABC_sampling[which(ABC_sampling[,1] > 0),], fileConn, row.names=F, col.names=F)
    ####################################################################################################
    ########                       4th block, used to define growth rate                        ########
    ####################################################################################################
    writeLines("// Growth rates", fileConn)
    writeLines(as.character(Growth_rate), fileConn)
    ####################################################################################################
    ########         5th block, used to define the number of migration matrices                 ########
    ####################################################################################################
    writeLines("// Number of migration matrices", fileConn)
    writeLines(as.character(Migration), fileConn)
    ####################################################################################################
    ########   6th block, used to define the parameters for the demographic process to model    ########
    ####################################################################################################
    writeLines("// Historical event: date, src, sink, %mig, new_pop_size, new_growth_rate, new_migration_matrix", fileConn)
    if(unique(ABC_table$Event_type) == "fast"){
      writeLines("1 event", fileConn)
      writeLines(paste(Event_time, Event_src, Event_sink, Event_mig, Event_new_Psize, Event_new_r, Event_new_mig, sep=" "), fileConn)
    }
    ####################################################################################################
    ########        7th block, used to define the mutation rate for the full population         ########
    ####################################################################################################
    writeLines("// mutation rate, average mutation number of mutations per generation per nucleotide, times the number of nucleotides", fileConn)
    writeLines(as.character(Mutation_rate), fileConn)
    ####################################################################################################
    ########          8th block, used to define the length of the sequences to simulate         ########
    ####################################################################################################
    writeLines("// Number of loci, sequnce length", fileConn)
    writeLines(as.character(nchar(ABC_DB$Sequence[ABC_DB$Species == species_vector[s]])[1]), fileConn)
    ####################################################################################################
    ########              9th block, used to define the type of DNA marker to simulate          ########
    ####################################################################################################
    writeLines("// Data type : either DNA, RFLP, or MICROSAT", fileConn)
    writeLines(paste("DNA", as.character(0.33333), sep=" "),  fileConn)
    ####################################################################################################
    ########                      10th block, used to define the mutation rate model            ########
    ####################################################################################################
    writeLines("// Mutation rates: gamma parameters, theta and k", fileConn)
    if(!is.na(Gamma)){
      writeLines(paste(Gamma, "10", sep=" "), fileConn)
    }
    if(is.na(Gamma)){
      writeLines(paste("0", "0", sep=" "), fileConn)
    }
    ####################################################################################################
    ########                    11th block, used to define the abstract priors                  ########
    ####################################################################################################
    writeLines("// Abstract priors", fileConn)
    writeLines(as.character(Abstract_prior), fileConn)
    close(fileConn)
  }
}
ABC_infile_single("species_constant.txt", "DATABASE.txt")
