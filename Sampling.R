
setwd("/Users/afr/Desktop/MSc_Ditte/Ditte_Data/")
par(mfrow=c(5,4), mar=c(2,3,3,1))
Full_DB_LL$Species
Sp_colors <- read.delim("Sp_color", stringsAsFactors=F, header=T)
for (sp in seq_along(Sp_colors$Species)){
  if(is.element(Sp_colors$Species[sp], unique(Full_DB_LL$Species))){
    temp_sp_db <- Full_DB[which(Full_DB$Species == Sp_colors$Species[sp] & Full_DB$Mean_Age <= 50000),]
    temp_hist_rec<- hist(temp_sp_db$Mean_Age, breaks=seq(0, 50000, 2000), plot=F)
    temp_hist_rec$counts[which(temp_hist_rec$counts == 0)] <- NA
    hist(temp_sp_db$Mean_Age, breaks=seq(0, 50000, 2000), main=NULL, xaxt="n", yaxt="n", ylab=NULL, xlab=NULL, col="#C1CDCD", border="#C1CDCD90")
    #labels=as.character(temp_hist_rec$counts) , 
    axis(side=2)
    mtext(side=3, Sp_colors$Species[sp], line=1) 
    temp_hist_seq <- hist(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1],breaks=seq(0, 50000, 2000), plot=F )
    temp_hist_seq$counts[which(temp_hist_seq$counts == 0)] <- NA
    hist(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1],breaks=seq(0, 50000, 2000), add=T, col=Sp_colors$Color_color[sp], border=Sp_colors$Color_color[sp])
    #text(x=temp_hist_seq$mids, y=0.5, labels=as.character(temp_hist_seq$counts), col="white", cex=0.7)
    abline(v=Sp_slow_fast$Slow[sp], col="#87CEFA90", lwd=10)
    Temp_count_seq_before_slow <- sum(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1] > Sp_slow_fast$Slow[sp])
    Temp_count_seq_after_slow <- sum(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1] < Sp_slow_fast$Slow[sp])
    text(Temp_count_seq_before_slow, x=Sp_slow_fast$Slow[sp]+2000, y=20, col="#87CEFA")
    text(Temp_count_seq_after_slow, x=Sp_slow_fast$Slow[sp]-2000, y=20, col="#87CEFA")
    abline(v=Sp_slow_fast$Fast[sp], col="#FF6A6A90", lwd=10)
    Temp_count_seq_before_fast <- sum(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1] > Sp_slow_fast$Fast[sp])
    Temp_count_seq_after_fast <- sum(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1] < Sp_slow_fast$Fast[sp])
    text(Temp_count_seq_before_fast, x=Sp_slow_fast$Fast[sp]+2000, y=20, col="#FF6A6A")
    text(Temp_count_seq_after_fast, x=Sp_slow_fast$Fast[sp]-2000, y=20, col="#FF6A6A")
    
    
  }
}
setwd("/Users/afr/Desktop/MSc_Ditte/Ditte_Data/")
par(mfrow=c(4,4), mar=c(2,3,3,1))
Sp_slow_fast <- read.delim("Species_slow_fast", stringsAsFactors=F, header=T)
for (sp in seq_along(Sp_slow_fast$Species)){
  if(is.element(Sp_slow_fast$Species[sp], unique(Full_DB_LL$Species))){
    temp_sp_db <- Full_DB[which(Full_DB$Species == Sp_slow_fast$Species[sp] & Full_DB$Mean_Age <= 50000),]
    temp_hist_rec <- hist(temp_sp_db$Mean_Age, breaks=seq(0, 50000, 2000), plot=F)
    temp_hist_rec$counts[which(temp_hist_rec$counts == 0)] <- NA
    temp_max <- max(na.omit(temp_hist_rec$counts))
    hist(temp_sp_db$Mean_Age, breaks=seq(0, 50000, 2000), main=NULL, xaxt="n", yaxt="n", ylab=NULL, xlab=NULL, col="#C1CDCD", border="#C1CDCD90")
    #labels=as.character(temp_hist_rec$counts) , 
    axis(side=2)
    mtext(side=3, Sp_slow_fast$Species[sp], line=1) 
    temp_hist_seq <- hist(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1],breaks=seq(0, 71000, 2000), plot=F )
    temp_hist_seq$counts[which(temp_hist_seq$counts == 0)] <- NA
    hist(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1],breaks=seq(0, 50000, 2000), add=T, col="#838B8B", border="#838B8B")
    #text(x=temp_hist_seq$mids, y=0.5, labels=as.character(temp_hist_seq$counts), col="white", cex=0.7)
    Event_1 <- 21000
    abline(v=Event_1 , col="#87CEFA90", lwd=5)
    Temp_count_seq_before_EV1 <- sum(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1] > Event_1 )
    Temp_count_seq_after_EV1 <- sum(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1] < Event_1 )
    text(Temp_count_seq_before_EV1, x=Event_1 +2500, y=temp_max*0.7, col="#87CEFA")
    text(Temp_count_seq_after_EV1, x=Event_1 -2500, y=temp_max*0.7, col="#87CEFA")
    Event_2 <- 9000
    abline(v=Event_2, col="#FF6A6A90", lwd=5)
    Temp_count_seq_before_EV2 <- sum(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1] > Event_2)
    Temp_count_seq_after_EV2 <- sum(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1] < Event_2)
    text(Temp_count_seq_before_EV2, x=Event_2+2500, y=temp_max*0.7, col="#FF6A6A")
    text(Temp_count_seq_after_EV2, x=Event_2-2500, y=temp_max*0.7, col="#FF6A6A")
    
    
  }
}










temp_sp_db <- Full_DB[which(Full_DB$Species == s),]
temp_hist_rec <- hist(temp_sp_db$Mean_Age, labels=T, breaks=seq(0, 71000, 2000), plot=F)
temp_hist_rec$counts[which(temp_hist_rec$counts == 0)] <- NA
hist(temp_sp_db$Mean_Age, labels=as.character(temp_hist_rec$counts) , breaks=seq(0, 71000, 2000), xlab="Time", main=NULL, xaxs="i", yaxs="i")
mtext(side=3, s, line=1)
temp_hist_seq <- hist(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1], labels=T,breaks=seq(0, 71000, 2000), plot=F )
temp_hist_seq$counts[which(temp_hist_seq$counts == 0)] <- NA
hist(temp_sp_db$Mean_Age[nchar(temp_sp_db$Sequence) > 1],breaks=seq(0, 71000, 2000), add=T, col="#838B8B")
text(x=temp_hist_seq$mids, y=0.5, labels=as.character(temp_hist_seq$counts), col="white", cex=0.7)

