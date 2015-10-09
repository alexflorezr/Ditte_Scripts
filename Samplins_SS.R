rm(list=ls())
setwd("/Users/afr/Desktop/Ch1_Figures/Figure_one")
Full_DB <- read.delim("DATABASE.txt", header=T, stringsAsFactors=F)
Full_DB$Longitude <- as.numeric(Full_DB$Longitude)
Full_DB$Latitude <- as.numeric(Full_DB$Latitude)
full_vector <- as.vector(NULL)
for (i in seq_along(Full_DB$Latitude)){
  if (is.na(Full_DB$Longitude[i]) || is.na(Full_DB$Latitude[i])){
    full_vector <- c(full_vector,i)
  }
}
# remove the records without longitude and or latitude
Full_DB_LL <- Full_DB[-full_vector,]
# remove the records older than 50000 years
DATABASE <- Full_DB_LL[-which(Full_DB_LL$Median_Age > 50000),]
############################################################################################################################
############################################################################################################################
require(ape)
require(pegas)
require(adegenet)
require(strataG)
require(sidier)
############################################################################################################################
############################################################################################################################
setwd("~/Desktop/MSc_Ditte/")
Sp_fast_slow <- read.delim("Velocity_4_species", stringsAsFactors=F, header=T)
Sp_fast_slow <- Sp_fast_slow[1,]
str(Sp_fast_slow)
#Sampling_aDNA <- function(Sp_fast_slow, database, ind_color){
Sp_bef_aft <- as.data.frame(matrix(nrow=dim(Sp_fast_slow)[1], ncol=7))
colnames(Sp_bef_aft) <- c("Species", "Bef", "Aft","Nuc_div_Total", "Phist", "Fst", "Nei_DA")
Event_name <- "Fast and slow events"
hist_rec_seq <- TRUE
hist_add_SS <- TRUE
hist_2gether <- FALSE
subsampling <- FALSE
plot_sub <- FALSE
pop_stats <- FALSE
adjacent_nuc <- TRUE
buffer <- 7000
setwd("~/Desktop/Ch1_Figures/Figure_one/Same_length_seqs/")
plot.new()
###########################################################################################
### Plotting all histograms together #######################################################################################################
if (hist_2gether == TRUE){
    nfrow <- as.vector(n2mfrow(length(Sp_fast_slow[,1])))
    par(mfrow=nfrow , mar=c(2,3,3,1), yaxs="i",xaxs="i",lwd=0.5)
}else{
    par(mar=c(4,3,3,1), yaxs="i",xaxs="i",lwd=0.5) 
}
### ENDS HERE #########################################################################################################################
################################################################################################################################
Nuc_bef_aft_raw <- as.data.frame(matrix(nrow=length(Sp_fast_slow$Species)*2, ncol= 5))
colnames(Nuc_bef_aft_raw) <- c("Species", "Nuc_bef", "Nuc_aft", "Nuc_diff", "Event_type")
for (sp in seq_along(Sp_fast_slow$Species)){
  ### STARTS HERE ##############################  
  if(is.element(Sp_fast_slow$Species[sp], unique(DATABASE$Species))){
    ## STARTS: histogram for the fossil record and sequences ##
    temp_sp_rec_db <- Full_DB[which(Full_DB$Species == Sp_fast_slow$Species[sp] & Full_DB$Mean_Age <= 50000),]
    temp_hist_rec <- hist(temp_sp_rec_db$Mean_Age, breaks=seq(0, 50000, 2000), plot=F)
    temp_hist_rec$counts[which(temp_hist_rec$counts == 0)] <- NA
    temp_max <- max(na.omit(temp_hist_rec$counts))
    temp_genus <- strsplit(Sp_fast_slow$Species[sp], split="_")[[1]][1]
    temp_specific <- strsplit(Sp_fast_slow$Species[sp], split="_")[[1]][2]
    main_name <- paste(temp_genus, temp_specific, sep=" ")
    sp_short_name <- paste(strsplit(temp_genus, split="")[[1]][1], strsplit(temp_specific, split="")[[1]][1], sep="")
    temp_sp_seq_db <- temp_sp_rec_db[nchar(temp_sp_rec_db$Sequence) > 1,]
    temp_hist_seq <- hist(temp_sp_seq_db$Mean_Age,breaks=seq(0, 71000, 2000), plot=F )
    temp_hist_seq$counts[which(temp_hist_seq$counts == 0)] <- NA
    sp_color <- "#00BFFF"
    #Events <- as.numeric(Sp_fast_slow[sp,c(2,3)])
    Events <- c(Sp_fast_slow$Fast[sp],Sp_fast_slow$Slow[sp])
    temp_sp <-  paste(strsplit(temp_genus, split="")[[1]][1], strsplit(temp_specific, split="")[[1]][1], sep="")
### STARTS HERE #########################################################################################################################
### Plotting histograms #######################################################################################################
    if( hist_rec_seq == TRUE){
      hist(temp_sp_rec_db$Mean_Age, breaks=seq(0, 50000, 2000), main=NULL, xaxt="n", yaxt="n", ylab=NULL, xlab=NULL, col="#00BFFF50", border=NA)
      axis(side=1)
      axis(side=2)
      mtext(side=3, main_name, line=1)
      hist(temp_sp_seq_db$Mean_Age,breaks=seq(0, 50000, 2000), add=T, col=paste(sp_color, 90, sep=""), border=NA)
      #labels=as.character(temp_hist_rec$counts)
      #if (ind_color == T) {sp_color <= paste(Sp_fast_slow$Color_color[sp], 90, sep="")}
      text(x=temp_hist_seq$mids, y=0.5, labels=as.character(temp_hist_seq$counts), cex=0.7, adj=c(0,0))
      x_fast <- c(Sp_fast_slow$Fast_end[sp],Sp_fast_slow$Fast_end[sp], Sp_fast_slow$Fast_start[sp],Sp_fast_slow$Fast_start[sp])
      x_fast_buffer <- x_fast+c(-buffer, -buffer, buffer, buffer)
      y_fast <- c(0,temp_max, temp_max, 0 )
      x_slow <- c(Sp_fast_slow$Slow_end[sp],Sp_fast_slow$Slow_end[sp], Sp_fast_slow$Slow_start[sp],Sp_fast_slow$Slow_start[sp])
      x_slow_buffer <- x_slow+c(-buffer, -buffer, buffer, buffer)
      y_slow <- c(0,temp_max, temp_max, 0 )
      polygon(x=x_fast_buffer, y = y_fast, col="#FF7F0099", border="#FF7F0099")
      polygon(x=x_fast, y = y_fast, col="#FF7F00", border="#FF7F00")
      polygon(x=x_slow_buffer, y = y_slow, col= "#EEC59199", border= "#EEC59199")
      polygon(x=x_slow, y = y_slow, col= "#EEC591", border= "#EEC591")
      
    }
### ENDS HERE #########################################################################################################################
############################################################################################################################
    temp_fasta_stats <- read.dna(paste(sp_short_name, "_fasta_new.fasta", sep=""), format = "fasta", as.character = FALSE, as.matrix=NULL)
    temp_data_stats <- as.data.frame(matrix(nrow=dim(temp_fasta_stats)[1], ncol=3))
#### STARTS HERE #########################################################################################################################
#### Estimating population summary statistics ########################################################################################################################
    if (pop_stats == TRUE){
      for (name in 1:dim(temp_data_stats)[1]){
        temp_data_stats[name,1] <- labels(temp_fasta_stats)[name]
        temp_data_stats[name,2] <- as.numeric(strsplit(labels(temp_fasta_stats), spli="_")[[name]][4])
        if( as.numeric(strsplit(labels(temp_fasta_stats), spli="_")[[name]][4]) > Events){
          temp_data_stats[name,3] <- "Pleistocene"
        }else{ temp_data_stats[name,3] <- "Holocene"
        }
      }
      colnames(temp_data_stats) <- c("name", "age", "era")
      haps <-FindHaplo(align=temp_fasta_stats, saveFile=F)
      as.data.frame(haps)
      colnames(haps)<-c("UniqueInd", "haplotype")
      full_data <- cbind(temp_data_stats, haps)
      if(length(unique(temp_data_stats[,3])) >= 2){
        if (sum(full_data[,1] == full_data[,4]) == length(full_data[,1])){
          N_haps<-length(unique(full_data$haplotype))###add number of haplotypes as one of statistics to be output
          sp_haps <-GetHaplo(align=temp_fasta_stats, saveFile=T, outname=paste(sp_short_name,"_Pops.fasta", sep=""), format="fasta", seqsNames="Inf.Hap") #haps are now a DNAbin
          sp_ghaps<-read.fasta(paste(sp_short_name,"_Pops.fasta", sep="")) #imports haps as gtypes
          sp_gtype<-gtypes(gen.data=data.frame(full_data$UniqueInd,full_data$era,full_data$haplotype),id.col=1,strata.col=2,locus.col=3,dna.seq=sp_ghaps)
          temp_nuc_div_total <- nuc.div(temp_fasta_stats, variance = F, pairwise.deletion = FALSE)
          temp_phist <- stat.phist(sp_gtype)$estimate
          temp_FST <- stat.fst(sp_gtype)$estimate
          temp_Nei_DA <- nei.Da(sp_gtype)$Da
          #Sp_stats[sp,] <- c(sp_short_name, temp_nuc_div_total, temp_phist, temp_FST, temp_Nei_DA)
        }else{
          print("Haplotype table and data table do not match the names")
        }
      }else{
        print(paste("no data for ", main_name, sep=""))
      }
    }
### ENDS HERE #########################################################################################################################
############################################################################################################################
    for (event in seq_along(Events)){
      temp_fasta <- read.FASTA(paste(tolower(temp_sp), "_fasta_new.fasta", sep=""))
      temp_vec <- c("temp_db_seqs_bef","temp_db_seqs_aft")
      start <- 1
      end <- dim(temp_db_seqs_bef)[1]
      if(adjacent_nuc == T){
        if (event == 1){
          Events_start_end <- c(Sp_fast_slow$Fast_start[sp],Sp_fast_slow$Fast_end[sp], Sp_fast_slow$Slow_start[sp],Sp_fast_slow$Slow_end[sp])
          temp_db_seqs_bef <- temp_sp_seq_db[(temp_sp_seq_db$Mean_Age < Events_start_end[1]+buffer) & (temp_sp_seq_db$Mean_Age > Events_start_end[1]),]
          temp_db_seqs_aft <- temp_sp_seq_db[(temp_sp_seq_db$Mean_Age > Events_start_end[2]-buffer)& (temp_sp_seq_db$Mean_Age < Events_start_end[2]),]
        }
        if (event == 2){
          Events_start_end <- c(Sp_fast_slow$Fast_start[sp],Sp_fast_slow$Fast_end[sp], Sp_fast_slow$Slow_start[sp],Sp_fast_slow$Slow_end[sp])
          temp_db_seqs_bef <- temp_sp_seq_db[(temp_sp_seq_db$Mean_Age < Events_start_end[3]+buffer) & (temp_sp_seq_db$Mean_Age > Events_start_end[3]),]
          temp_db_seqs_aft <- temp_sp_seq_db[(temp_sp_seq_db$Mean_Age > Events_start_end[4]-buffer)& (temp_sp_seq_db$Mean_Age < Events_start_end[4]),]
        }
      }
      for (bef_aft in seq_along(temp_vec)){
        temp_count_seq <- dim(get(temp_vec[bef_aft]))[1]
        temp_acc <- get(temp_vec[bef_aft])$Accession_GB
        temp_age <- get(temp_vec[bef_aft])$Mean_Age
        temp_err <- get(temp_vec[bef_aft])$Cal_Sigma
        temp_rep <- get(temp_vec[bef_aft])$Repeat
        temp_Coor <- cbind(get(temp_vec[bef_aft])$Latitude, get(temp_vec[bef_aft])$Longitude)
        Temp_seq_name_fasta <- paste(temp_sp,temp_acc, temp_rep, temp_age,temp_err, sep="_") 
        temp_seqs  <- temp_fasta[as.vector(na.exclude(match(Temp_seq_name_fasta, labels(temp_fasta))))]
        ## pop gen SS ##
        if (length(temp_seqs) > 0){
          temp_nuc_div <- nuc.div(temp_seqs)
          temp_num_hap <- haplotype(temp_seqs)
          temp_segSites <- seg.sites(temp_seqs)
          temp_theta  <- theta.s(length(temp_segSites), length(temp_seqs))
          temp_pairdiff <- mean(as.vector(dist.dna(temp_seqs, model = "raw")))
          
        } else{
          temp_nuc_div <- -999
          temp_num_hap <- as.matrix(-999, -999)
          temp_segSites <- -999
          temp_theta  <- -999
          temp_pairdiff <- -999
        }
        if (temp_vec[bef_aft] == "temp_db_seqs_bef"){
          temp_seqs_bef <- temp_seqs
          SS_bef <- c(temp_count_seq, temp_pairdiff, temp_nuc_div, round(dim(temp_num_hap)[1]), round(length(temp_segSites)), temp_theta)
          temp_Coor_bef <- cbind(get(temp_vec[bef_aft])$Latitude, get(temp_vec[bef_aft])$Longitude)
          }
        if (temp_vec[bef_aft] == "temp_db_seqs_aft"){
          temp_seqs_aft <- temp_seqs
          SS_aft <- c(temp_count_seq, temp_pairdiff, temp_nuc_div, round(dim(temp_num_hap)[1]), round(length(temp_segSites)), temp_theta)
          temp_Coor_aft <- cbind(get(temp_vec[bef_aft])$Latitude, get(temp_vec[bef_aft])$Longitude)
          }
      }
### STARTS HERE #########################################################################################################################
### Adding Summay statistics to the histograms #######################################################################################################
    if( hist_add_SS == TRUE){
      for (ss_bef in seq_along(SS_bef)){
        text(round(SS_bef[ss_bef], digits=3), x=Events[event]+2000, y=(temp_max/2.5)+sum(rep(temp_max/12, times=ss_bef)), col="#363636", adj=c(0,0), cex=0.7)
      }
      for (ss_aft in seq_along(SS_aft)){
        text(round(SS_aft[ss_aft], digits=3), x=Events[event]-2000, y=(temp_max/2.5)+sum(rep(temp_max/12, times=ss_aft)), col="#363636", adj=c(1,0), cex=0.7)
      }
      SS_name_vector <- c("n","dist", "nuc div", "hap", "seg", "theta")
      for (name in seq_along(SS_name_vector)){
        if((min(Events) > 10000) & (!is.na(min(Events)))){
        text(SS_name_vector[name], x=5000, y=(temp_max/2.5)+sum(rep(temp_max/12, times=name)), col="#363636", adj=c(1,0), cex=0.7)
        }
        if((min(Events) <= 10000) | is.na(min(Events))){
        text(SS_name_vector[name], x=40000, y=(temp_max/2.5)+sum(rep(temp_max/12, times=name)), col="#363636", adj=c(0,0), cex=0.7)
        }
      }
    }
    if (event == 1 ){
      Nuc_bef_aft_raw[(sp*2)-1,] <- c(Sp_bef_aft$Species[sp],SS_bef[3], SS_aft[3],  (SS_bef[3] - SS_aft[3]), "Fast")
    }
    if (event == 2 ){
      Nuc_bef_aft_raw[(sp*2),] <- c(Sp_bef_aft$Species[sp],SS_bef[3], SS_aft[3],  (SS_bef[3] - SS_aft[3]), "Slow")
    }
    }
    
### ENDS HERE #########################################################################################################################
############################################################################################################################
### STARTS HERE #########################################################################################################################
### Subsampling sequences to estimate nucleotide diversity with equal sampling #########################################################################################################################
    if (subsampling == TRUE){
      nucdiv_sub <- as.numeric(vector(length=10000))
      if (length(temp_seqs_bef) > length(temp_seqs_aft)){ 
        for (sub in 1:10000){
          sub_seqs <- sample((1:length(temp_seqs_bef)), size=length(temp_seqs_aft), replace=T)
          nucdiv_sub[sub] <- nuc.div(temp_seqs_bef[sub_seqs])
        }
      }
      if (length(temp_seqs_bef) < length(temp_seqs_aft)){
        for (sub in 1:10000){
          sub_seqs <- sample((1:length(temp_seqs_aft)), size=length(temp_seqs_bef), replace=T)
          nucdiv_sub[sub] <- nuc.div(temp_seqs_aft[sub_seqs])
        }
      }
      if (length(temp_seqs_bef) == length(temp_seqs_aft)){
        temp_seqs_bef <- temp_seqs_bef[-sample((1:length(temp_seqs_bef)), 1)]
        for (sub in 1:10000){
          sub_seqs <- sample((1:length(temp_seqs_aft)), size=length(temp_seqs_bef), replace=T)
          nucdiv_sub[sub] <- nuc.div(temp_seqs_aft[sub_seqs])
        }
      }  
    }
### ENDS HERE #########################################################################################################################
############################################################################################################################
    
    if (colnames(Sp_fast_slow)[1+event] == "Fast"){
      col_event<- "#FFD39B"
    }
    if (colnames(Sp_fast_slow)[1+event] == "Slow"){
      col_event<- "#FFD39B"
    }
### STARTS HERE #########################################################################################################################
### Plotting histograms for the nucleotide diversity based on the subsampling #########################################################################################################################
    if (plot_sub == TRUE){
      if (event == 1){
        temp_hist_sub <- hist(nucdiv_sub, plot=F)
        temp_max_y_sub <- max(temp_hist_sub$counts)
        
        temp_max_x_sub <- max(temp_hist_sub$mids)
        temp_max_x_sub <- 0.035
        hist(nucdiv_sub, col=ifelse(length(temp_seqs_bef) > length(temp_seqs_aft),"#8B735580", "#FF8C0090" ),border=ifelse(length(temp_seqs_bef) > length(temp_seqs_aft),"#8B735580", "#FF8C0090" ), main=paste(main_name, " (", Event_name, ")", sep=""), xlab="Nucleotide diversity")
        abline(v=c(mean(nucdiv_sub), mean(nucdiv_sub)+sd(nucdiv_sub), mean(nucdiv_sub)-sd(nucdiv_sub)), lwd=c(3,3,3), lty=c(2,3,3), col=ifelse(length(temp_seqs_bef) > length(temp_seqs_aft),"#8B735580", "#FF8C0090" ))
        abline(v=SS_bef[2], lwd=6, col="#8B7355")
        text(paste("Before", " (", length(temp_seqs_bef),")", sep="") ,x=temp_max_x_sub*0.96, y=temp_max_y_sub*0.8, col="#8B7355", adj=c(0,0.5))
        #segments(x0=temp_max_x_sub*0.95, x1=temp_max_x_sub*0.91, y0=temp_max_y_sub*0.8, y1=temp_max_y_sub*0.8, col="#8B7355", lwd=6)
        abline(v=SS_aft[2], lwd=6, col="#FF8C00")
        text(paste("After", " (", length(temp_seqs_aft),")", sep=""),x=temp_max_x_sub*0.96, y=temp_max_y_sub*0.75, col="#FF8C00", adj=c(0,0.5))
        #segments(x0=temp_max_x_sub*0.95, x1=temp_max_x_sub*0.91, y0=temp_max_y_sub*0.75, y1=temp_max_y_sub*0.75, col="#FF8C00", lwd=6)
      }
      #if (event == 2){
      #hist(nucdiv_sub, col=paste(col_event, 50, sep=""), border=NA, xlim=c(0,0.07), main=main_name)
      #abline(v=c(mean(nucdiv_sub), mean(nucdiv_sub)+sd(nucdiv_sub), mean(nucdiv_sub)-sd(nucdiv_sub)), lwd=c(3,2,2), lty=c(2,3,3), col=col_event)
      #abline(v=SS_bef[2], lwd=4, col=col_event)
      #abline(v=SS_aft[2], lwd=4, col="yellow") 
      
    }
### ENDS HERE #########################################################################################################################
############################################################################################################################
    
  }
  #Sp_bef_aft[sp,] <- c(Sp_fast_slow$Species[sp], round(SS_bef[2], digits=4), round(SS_aft[2], digits=4), c(temp_nuc_div_total, temp_phist, temp_FST, temp_Nei_DA)) 
}


### STARTS HERE #########################################################################################################################
############################################################################################################################
plot_nuc_bef_aft  <- Nuc_bef_aft_raw
plot_nuc_bef_aft <- plot_nuc_bef_aft[-which(plot_nuc_bef_aft$Nuc_bef == -999),]
boxplot(as.numeric(plot_nuc_bef_aft$Nuc_diff) ~ as.factor(plot_nuc_bef_aft$Event_type), plot=F)
plot(as.numeric(plot_nuc_bef_aft$Nuc_diff), col=ifelse(plot_nuc_bef_aft$Event_type == "Fast", "red", "green"), pch=16, cex=2)
for (spe in seq_along(plot_nuc_bef_aft$Species)){
  text(x=spe, y=as.numeric(plot_nuc_bef_aft$Nuc_diff)[spe], plot_nuc_bef_aft$Species[spe], srt=90, adj=c(1,0), cex=0.8)
}
  







temp_plot_by <- as.data.frame(matrix(nrow=length(Sp_bef_aft[,1]),ncol=8))
temp_plot_by[,c(1,2,3)] <- Sp_bef_aft[,c(1,2,3)]
temp_plot_by[,c(5,6,7,8)] <- Sp_bef_aft[,c(4,5,6,7)]
temp_plot_by[,4] <- as.numeric(Sp_bef_aft[,3]) - as.numeric(Sp_bef_aft[,2])
colnames(temp_plot_by) <- c(colnames(Sp_bef_aft)[c(1,2,3)], "Diff", colnames(Sp_bef_aft)[c(4,5,6,7)])
# sort by Nuc div difference
Plot_bef_aft_by_diff <- temp_plot_by[order(as.numeric(temp_plot_by$Diff), decreasing=F),]
###### for a simple bar plot #####
sp_3_colors <- read.delim("Sp_3_colors.txt", sep="\t", header=T, stringsAsFactors = F)
plot(x=seq(1,dim(Plot_bef_aft_by_diff)[1], by=1), ylim=c(-0.01,0.082), xlim=c(0,dim(Plot_bef_aft_by_diff)[1]+2), frame = F, xaxt="n", xlab=NA, ylab=NA, yaxt="n")
Title <- "Change in nucleotide diversity in the Holarctic"
mtext(Title, line=1, cex=1.4, col="#3D3D3D")
width <- 0.85
for (nuc in seq_along(Plot_bef_aft_by_diff[,1])){
  if (as.numeric(Plot_bef_aft_by_diff[nuc,3]) > 0){
    if (as.numeric(Plot_bef_aft_by_diff[nuc,2]) > as.numeric(Plot_bef_aft_by_diff[nuc,3])){
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(0,Plot_bef_aft_by_diff[nuc,3], Plot_bef_aft_by_diff[nuc,3],0), col="#FF8C00", border=NA)
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(Plot_bef_aft_by_diff[nuc,3],Plot_bef_aft_by_diff[nuc,2],Plot_bef_aft_by_diff[nuc,2],Plot_bef_aft_by_diff[nuc,3]), col="#8B5A2B99", border=NA)
      sq_col <- sp_3_colors$Color[which(Plot_bef_aft_by_diff[nuc,1] == sp_3_colors[,1])]
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(-0.0035, -0.0003, -0.0003, -0.0035), col = paste(sq_col, 99, sep=""), lwd=2, border=NA)
    }
    if (as.numeric(Plot_bef_aft_by_diff[nuc,2]) < as.numeric(Plot_bef_aft_by_diff[nuc,3])){
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(Plot_bef_aft_by_diff[nuc,3],Plot_bef_aft_by_diff[nuc,2],Plot_bef_aft_by_diff[nuc,2],Plot_bef_aft_by_diff[nuc,3]), col="#FF8C00", border=NA)
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(0,Plot_bef_aft_by_diff[nuc,2], Plot_bef_aft_by_diff[nuc,2],0), col="#8B5A2B99", border=NA)
      sq_col <- sp_3_colors$Color[which(Plot_bef_aft_by_diff[nuc,1] == sp_3_colors[,1])]
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(-0.0035, -0.0003, -0.0003, -0.0035), col = paste(sq_col, 90, sep=""), lwd=2, border=NA)
      
    }
  }
  if (as.numeric(Plot_bef_aft_by_diff[nuc,3]) < 0){
    polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(0,Plot_bef_aft_by_diff[nuc,2], Plot_bef_aft_by_diff[nuc,2],0), col="#8B5A2B99", border=NA)
    #text(x=nuc+width/2, y=0 ,"T", col="#FF8C00", cex=1.5, adj=c(0.5,0))
    sq_col <- sp_3_colors$Color[which(Plot_bef_aft_by_diff[nuc,1] == sp_3_colors[,1])]
    polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(-0.0035, -0.0003, -0.0003, -0.0035), col = paste(sq_col, 90, sep=""), lwd=2, border=NA)
    
  }
  
  sp_letters <- paste(strsplit(strsplit(Plot_bef_aft_by_diff[nuc,1], split=c("_"))[[1]][1], split="")[[1]][1], strsplit(strsplit(Plot_bef_aft_by_diff[nuc,1], split=c("_"))[[1]][2], split="")[[1]][1], sep="")
  text(y=-0.005, x=nuc+(width/3), sp_letters, adj=c(0,1), cex=0.8)
}
polygon(x=c(2,2,3,3), y=c(0.07,0.073,0.073,0.07), col="#FF8C00",border=NA)
text("Holocene", x=3.1, y=0.07, adj=c(0, 0), cex=1.2, col="#3D3D3D")
polygon(x=c(2,2,3,3), y=c(0.066,0.069,0.069,0.066), col="#8B5A2B99",border=NA)
text("Late Pleistocene", x=3.1, y=0.066, adj=c(0, 0), cex=1.2, col="#3D3D3D")
for (line in c(0, 0.01, 0.02, 0.04, 0.08)){
  segments(x0=0.5,x1= nuc+width+phist+0.1, y0=line, y1=line, col="#3D3D3D", lwd=2, lty=3)
}
text(x=0.7, y=c(0,0.01, 0.02, 0.04, 0.08),  c(0,0.01, 0.02, 0.04, 0.08), cex=0.7, adj=c(0.5,-0.5))
mtext(side=2, "Nucleotide diversity", col="#3D3D3D", line=-1, cex=1.5)


###### for a plot using stats #####
Title <- "Change in nucleotide diversity and temporal population structure (FST)"
plot(x=seq(1,dim(Plot_bef_aft_by_diff)[1], by=1), ylim=c(-0.01,0.082), xlim=c(0,dim(Plot_bef_aft_by_diff)[1]+2), frame = F, xaxt="n", xlab=NA, ylab=NA, yaxt="n")
mtext(Title, line=1, cex=1.4)
width <- 0.4
for (nuc in seq_along(Plot_bef_aft_by_diff[,1])){
  # Using phist
  #phist <- as.numeric(Plot_bef_aft_by_diff[nuc,5])*6.5
  # Using FST
  phist <- as.numeric(Plot_bef_aft_by_diff[nuc,6])/2
  polygon(x=c(nuc+phist,nuc+phist ,nuc+width+phist,nuc+width+phist), y=c(0,Plot_bef_aft_by_diff[nuc,3], Plot_bef_aft_by_diff[nuc,3],0), col="#FF8C00", border=NA)
  polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(0,Plot_bef_aft_by_diff[nuc,2],Plot_bef_aft_by_diff[nuc,2],0), col="#8B5A2B99", border=NA)
  midpoint <- nuc+((nuc+width+phist-nuc)/2)
  polygon(x=c(midpoint, nuc, nuc+width+phist), y=c(-0.002, -0.0003, -0.0003), border="#3D3D3D", lwd=2)
  sp_letters <- paste(strsplit(strsplit(Plot_bef_aft_by_diff[nuc,1], split=c("_"))[[1]][1], split="")[[1]][1], strsplit(strsplit(Plot_bef_aft_by_diff[nuc,1], split=c("_"))[[1]][2], split="")[[1]][1], sep="")
  text(y=-0.005, x=midpoint, sp_letters, adj=c(0.5,0), cex=0.8)
}
polygon(x=c(2,2,3,3), y=c(0.07,0.073,0.073,0.07), col="#FF8C00",border=NA)
text("Holocene", x=3.1, y=0.07, adj=c(0, 0), cex=1.2, col="#3D3D3D")
polygon(x=c(2,2,3,3), y=c(0.066,0.069,0.069,0.066), col="#8B5A2B99",border=NA)
text("Late Pleistocene", x=3.1, y=0.066, adj=c(0, 0), cex=1.2, col="#3D3D3D")
for (line in c(0.01, 0.02, 0.04, 0.08)){
  segments(x0=0.5,x1= nuc+width+phist+0.1, y0=line, y1=line, col="#3D3D3D", lwd=2, lty=3)
}
text(x=0.7, y=c(0,0.01, 0.02, 0.04, 0.08),  c(0,0.01, 0.02, 0.04, 0.08), cex=0.7, adj=c(0.5,-0.5))
