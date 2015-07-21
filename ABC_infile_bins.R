rm(list=ls())
#====== Create input files using species database =====
ABC_infile_bins <- function(ABC_table_species, Database){
#### Only one population will be defined
#### Define the variables for the function
ABC_table <- read.delim(ABC_table_species ,header=T, sep="\t", stringsAsFactors=FALSE, na.strings="NA")
species_vector <- as.vector(ABC_table$Species)
for (s in seq_along(species_vector)){
    for (i in seq_along(colnames(ABC_table))){
      assign(colnames(ABC_table)[i], ABC_table[s,i])
    }
# Create a conection (e.i. an empty file) to write the text
  fileConn <-file(paste(Species, Npop, Event_type, paste(Event_time, ".par", sep=""), sep="_"), "w")
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
# Create an empty matrix to save the values for the sampling groups
  ABC_sampling <- matrix(NA, nrow=10, ncol=4)
# Send the sequence sampling to the first column
  temp_hist <- hist(ABC_DB$Median_Age[which(ABC_DB$Species == species_vector[s] & nchar(ABC_DB$Sequence) > 2)], breaks=seq(0,50000,5000),plot=F)
  ABC_sampling[,1] <- temp_hist$counts
# Send the age of the bins divided by the generation time to the second column
  ABC_sampling[,2] <- temp_hist$breaks[-11]/Gen_time
# Only one population 0
  ABC_sampling[,3] <- rep(0, length(ABC_sampling[,3]))
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
  writeLines(Growth_rate, fileConn)
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
    writeLines(paste(Event_time, Event_src, Event_sink, Event_mig, Event_new_Psize, Event_new_r, sep=" "), fileConn)
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
########                      10th block, used to define the mutation rate                  ########
####################################################################################################
  writeLines("// Mutation rates: gamma parameters, theta and k", fileConn)
  if(unique(ABC_table$Mutation_model) == "HKY"){
    writeLines("0.4 10", fileConn)
  }
####################################################################################################
########                    11th block, used to define the abstract priors                  ########
####################################################################################################
writeLines("// Abstract priors", fileConn)
writeLines(ABC_table$Abstract_prior, fileConn)
close(fileConn)

}
}
ABC_infile_bins("species_fast.txt", "DATABASE.txt")
