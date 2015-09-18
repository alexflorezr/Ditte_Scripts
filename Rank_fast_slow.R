rm(list=ls())
library(raster)
library(rgdal)
# Upload the climate velocity files
setwd("/Users/afr/Desktop/Regression/Cli_vel/")
cli_vel_Tmp <- stack("climate_change_velocity_perYear_tmp_rescaled_truncated_0_000005.grd")
cli_vel_Prc <- stack("climate_change_velocity_perYear_prec_rescaled_truncated_0_000005.grd")
cli_vel_TnP <- stack("climate_change_velocity_perYear_tmp_prec_rescaled_truncated_0_000005.grd")
# Upload the biomes files
setwd("/Users/afr/Desktop/Regression/Biomes/Koppen_full_cor/")
biome_kpp <- stack(sort(dir(),decreasing=T))
setwd("/Users/afr/Desktop/Regression/Biomes/Bio4_CO2/")
biome_bio <- stack(sort(dir(),decreasing=T))
# objects describing the variables
## Vector for the time bins from the present to the past
bins_0to50 <- c(seq(0,21000, by=1000), seq(22000, 50000, by=2000))
bins_50to0 <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
layers <- seq(25, 61, by=1)
## Vector for the time bins middle point from the present to the past
bins_mid_point <- c(seq(-500, 21500, by=1000), seq(23000, 49000, by=2000))
# main loop for extracting values for every species
## Define the species to use in the loop
setwd("/Users/afr/Desktop/Regression")
Single_sp <- read.delim("Single_sp_regression", header=T, stringsAsFactors=F)
#########################################################################################################################
#########################################################################################################################
##### Read the database
## TODO ## Convert this part of the script in a function
#########################################################################################################################
#########################################################################################################################
Full_DB <- read.delim(file.choose(), header=T, stringsAsFactors=F)
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
Full_DB_LL <- Full_DB_LL[-which(Full_DB_LL$Median_Age > 50000),]
#########################################################################################################################
#########################################################################################################################
##### Extract the climatic values for every fossil record in the database
## TODO ## Convert this part of the script in a function
#########################################################################################################################
#########################################################################################################################
species_set_clim <- as.data.frame(matrix(nrow=0, ncol=12))
species_set_BSP <- as.data.frame(matrix(nrow=0, ncol=6))
for (s in seq_along(Single_sp$Species)){
  # Database to work with, contains all the species and data
  temp_DB_climate <- Full_DB_LL[which(Full_DB_LL$Species==Single_sp$Species[s]),]
  # Empty data frame to save the data for every species
  temp_points <- as.data.frame(matrix(nrow=nrow(temp_DB_climate), ncol=13))
  colnames(temp_points) <- c("Species","Longitude", "Latitude", "Mean_date_record", "Time_bin", "Layer", "Rec_type", "Cell", "Cli_vel_tmp", "Cli_vel_prc", "Cli_vel_tnp", "Biome_kpp","Biome_bio")
  counter1 <- 1
  for (bin in seq_along(bins_50to0)[-37]){
    records <- which(temp_DB_climate$Mean_Age <= bins_50to0[bin] & temp_DB_climate$Mean_Age > bins_50to0[bin+1])
    points <- as.matrix(cbind.data.frame(as.numeric(temp_DB_climate$Longitude[records]), as.numeric(temp_DB_climate$Latitude[records])))
    colnames(points) <- c("longitude", "Latitude")
    if(length(records > 0)){
      temp_points[(counter1:(counter1 + length(records)-1)),1] <- Single_sp$Species[s]
      temp_points[(counter1:(counter1 + length(records)-1)),2:7] <- cbind.data.frame(as.numeric(temp_DB_climate$Longitude[records]), as.numeric(temp_DB_climate$Latitude[records]),as.numeric(temp_DB_climate$Mean_Age[records]), rep(bins_50to0[bin], times=length(records)), rep(layers[bin], times=length(records)),ifelse(nchar(temp_DB_climate$Sequence[records]) > 1, "Seq", "Fossil") )
      temp_points[(counter1:(counter1 + length(records)-1)),c(which(colnames(temp_points) == "Cell"), which(colnames(temp_points) == "Cli_vel_tmp")) ] <- extract(cli_vel_Tmp, layer=layers[bin], nl=1, y=points, cellnumbers=T)
      temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == "Cli_vel_prc")]<- extract(cli_vel_Prc, layer=layers[bin], nl=1, y=points)
      temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == "Cli_vel_tnp")]<- extract(cli_vel_TnP, layer=layers[bin], nl=1, y=points)
      #temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == "Biome_kpp")]<- extract(biome_kpp, layer=(layers[bin])+1, nl=1, y=points)
      #temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == "Biome_bio")]<- extract(biome_bio, layer=(layers[bin])+1, nl=1, y=points)
      counter1 <- counter1+length(records)
    }
  } 
  assign(paste(Single_sp$Sp[s], "points", sep="_"), temp_points)
  if (Species_set_clim ==TRUE){
    species_set_clim <- rbind(species_set_clim, get(paste(Single_sp$Sp[s], "points", sep="_")))
  }
}

####### geometric mean estimates for species per bin #####
### working !!!! ####
### using geometric mean ###
g_mean_species <- matrix(nrow=length(bins_50to0), ncol=length(Single_sp$Species)+1)
rownames(g_mean_species) <- bins_50to0
colnames(g_mean_species) <- c(Single_sp$Species, "n")
g_mean_species[,18] <- as.vector(rep(0, length(bins_50to0)))
for (species in seq_along(colnames(g_mean_species))[-18]){
  for (bin in seq_along(rownames(g_mean_species))){
    g_mean_species[bin,species] <- exp(mean(log(na.omit(species_set_clim$Cli_vel_prc[species_set_clim$Species == colnames(g_mean_species)[species] & species_set_clim$Time_bin == rownames(g_mean_species)[bin]]))))
    g_mean_species[bin,length(Single_sp$Species)+1] <- g_mean_species[bin,length(Single_sp$Species)+1] + length(na.omit(species_set_clim$Cli_vel_prc[species_set_clim$Species == colnames(g_mean_species)[species] & species_set_clim$Time_bin == rownames(g_mean_species)[bin]]))
  }
}
sum(g_mean_species[which(as.numeric(rownames(g_mean_species)) > 21000,18])
### Plot the geometric mean per bin per species, using different colors for herbivorous, carnivorous, and rodents####
### working !!!! ####
for (sp in seq_along(colnames(g_mean_species))){
  if (sp == 1){
    plot(as.numeric(rownames(g_mean_species)), g_mean_species[,sp], pch=15, cex=1.3, col=paste(Single_sp$Sp_color[which(Single_sp$Species == colnames(g_mean_species)[sp])], "99", sep=""), frame=F, xlab="Time (years BP)", ylab=("Climate velocity (km/year)"), ylim=c(-0.05,0.5), xaxt="n")
    axis(side=1, line=-1)
    text(x=as.numeric(rownames(g_mean_species)), y=-0.03, labels=g_mean_species[,18], srt=90, cex=0.5, adj=c(0,0.5))
  } 
  else{
    points(as.numeric(rownames(g_mean_species)), g_mean_species[,sp], pch=15,cex=1.3, col=paste(Single_sp$Sp_color[which(Single_sp$Species == colnames(g_mean_species)[sp])], "99", sep=""))
  }
}
points(rownames(g_mean_region)[37:1], g_mean_region[37:1,1], col="black", pch=17)
abline(h=0.25, col="#000080", lwd=3)
text("Polar region",  col="#000080", x=50000, y=0.26, cex=0.85, adj=c(1,0))
abline(h=0.45, col="#1874CD", lwd=3)
text("Cold region",  col="#1874CD", x=50000, y=0.46, cex=0.85, adj=c(1,0))
abline(h=0.20, col="#40E0D0", lwd=3)
text("Temperate region",  col="#40E0D0", x=50000, y=0.21, cex=0.85, adj=c(1,0))
### plot variation in regional climate velocity ####
### Number of cells with climate velocity values per time bin ####
install.packages("rgl")
install.packages("plot3D")
library(rgl)
library(plot3D)
cell_bin <- matrix(nrow=length(bins_50to0), ncol=2)
clim_vel_bin <- matrix(nrow=length(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])]), ncol=0)
rownames(cell_bin) <- bins_50to0
colnames(cell_bin) <- c("Bin", "Cells")
for (point in seq_along(bins_50to0)){
  cell_bin[point,] <- (c(bins_50to0[point],length(na.omit(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])]))))
  clim_vel_bin <- cbind.data.frame(clim_vel_bin, cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])])
}
colnames(clim_vel_bin) <- bins_50to0
par(mar=c(7,7,7,7))
barplot(cell_bin[37:1,2], border=NA, col="#CDAA7D", ylim=c(0,12000),xaxs="i", yaxs="i",names.arg=as.character(cell_bin[37:1,1]), las=2)
mtext("Number of cells with climate velocity values per time bin", side=3, line=2, cex=1.5)

### histogram for the regional values ####
bins_col_func <- colorRampPalette(c("#8B4500", "#8B6508", "#FF7F00","#FFA500","#698B22", "#90EE90", "#B4EEB4"))
bins_col <- bins_col_func(37)
hist_breaks <- seq(log(0.000000001),log(2000), 0.1)
par(lwd=0.5)
for (point in seq_along(bins_50to0)){
  if (point <= 1){
  hist(log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])] + 0.00000001), breaks=hist_breaks, xaxt="n",yaxt="n",col=paste(bins_col[point], 90, sep=""), border="white", lwd=0.1, ylim=c(0,0.5),xlim=log(c(0.00001, 1000)), freq=F, yaxs="i", xaxs="i", main=NULL, xlab=NULL, ylab=NULL)
  temp_d <- density(log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])] + 0.00000001), na.rm=T)
  polygon(temp_d, col="red", border="blue")
  axis(side=2)
  axis(side=1, at=log(c(0.000001,0.00001, 0.0001, 0.001,0.01,0.1,1,10, 100, 1000)), labels=c(0.000001,0.00001, 0.0001, 0.001,0.01,0.1,1,10, 100, 1000), line=-0.9)
  mtext(side=3, "Climate velocity through time for the Holarctic region", line=2, cex=2)
  mtext(side=1, "Climate velocity (km/yr)", line=2)
  
  }
  else {
  hist(log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])] + 0.00000001), breaks=hist_breaks, xaxt="n",col=paste(bins_col[point], 90, sep=""), border="white", lwd=0.1, ylim=c(0,0.5),xlim=log(c(0.00001, 1000)), freq=F, yaxs="i", xaxs="i", main=NULL, xlab=NULL, add=T)
  }
  }
legend(x=log(50), y=0.5, fill=bins_col, legend=bins_50to0, cex=0.45, border=F, bty="n", x.intersp=0.1)
segments(x0=log(0.25), x1=log(0.25), y0=0,y1=0.5, lwd=2, lty=2, col="#00008B")
text(x=log(0.18),y=0.35,"Polar region", srt=90, adj=c(0,1),col="#00008B", cex=0.7)
segments(x0=log(0.26), x1=log(0.26), y0=0,y1=0.5, lwd=2, lty=2, col="#4F94CD")
text(x=log(0.28),y=0.35,"Temperate region", srt=90,adj=c(0,1), col="#4F94CD", cex=0.7)
segments(x0=log(0.4), x1=log(0.4), y0=0,y1=0.5, lwd=2, lty=2, col="#00CDCD")
text(x=log(0.45),y=0.35,"Cold region", srt=90,adj=c(0,1), col="#00CDCD", cex=0.7)



### Density comparison for the regional values ####
install.packages("sm")
library(sm)
for_density <- na.omit(cli_vel_tmp_hol)
for_density[,1] <- log(for_density[,1]+ 0.00000001)
sm.density.compare(for_density[,1], for_density[,2], col=bins_col, lty=1)
### Density for the regional values ####
bins_col_brown <- colorRampPalette("#8B4513")
bins_col_orange  <- colorRampPalette("#FF7F00")
bins_col_yellow  <- colorRampPalette("#EEC900")
bins_col_green <- colorRampPalette("#008B00")
bins_col_grey <- colorRampPalette("#8EE5EE")
bins_col <- c(bins_col_brown(19), bins_col_orange(8), bins_col_yellow(4), bins_col_green(2), bins_col_grey(4))
hist_breaks <- seq(log(0.000000001),log(2000), 0.1)
par(lwd=0.5)
plot.new()
par(mar=c(7,7,7,7))
for (point in seq_along(bins_50to0)){
  if (point <= 1){
    temp_d <- density(log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])] + 0.00000001), na.rm=T)
    plot(temp_d, col=paste(bins_col[point], 50, sep=""), xlim=log(c(0.00001, 100)), ylim=c(-0.05,0.5), main=NA, frame=F, xaxt="n", xaxs="i", yaxs="i", xlab="")
    polygon(temp_d, col=paste(bins_col[point], 50, sep="") , border=bins_col[point], lwd=2)
    x3=mean(log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])] + 0.00000001), na.rm=T)
    segments(x0=x3, x1=x3, y0=0, y1=-0,047, col=bins_col[point], lwd=3)
    axis(side=1, at=log(c(0.000001,0.00001, 0.0001, 0.001,0.01,0.1,1,10, 100, 1000)), labels=c(0.000001,0.00001, 0.0001, 0.001,0.01,0.1,1,10, 100, 1000), line=0.001)
    mtext(side=1, "Climate velocity (km/yr)", line=3)
    axis(side=2)
    mtext(side=3, "Climate velocity through time for the Holarctic region", line=3, cex=2)
    
    
  }
  else {
    temp_d <- density(log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])] + 0.00000001), na.rm=T)
    polygon(temp_d, col=paste(bins_col[point], 50, sep="") , border=bins_col[point] , lwd=2)
    x3=mean(log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])] + 0.00000001), na.rm=T)
    segments(x0=x3, x1=x3, y0=0, y1=-0.047, col=bins_col[point], lwd=3, lty=1)
    }
}
legend(x=log(20), y=0.5, fill=bins_col, legend=bins_50to0, cex=0.35, border=F, bty="n", x.intersp=0.1)
segments(x0=log(0.25), x1=log(0.25), y1=0, y0=-0.047, lwd=3, lty=3, col="#00008B")
#text(x=log(10),y=-0.005,"Polar region", adj=c(0,1),col="#00008B", cex=0.7)
segments(x0=log(0.26), x1=log(0.26), y1=0, y0=-0.047, lwd=3, lty=3, col="#4F94CD")
#text(x=log(10),y=-0.017,"Temperate region",adj=c(0,1), col="#4F94CD", cex=0.7)
segments(x0=log(0.4), x1=log(0.4), y1=0, y0=-0.047, lwd=3,lty=3, col="#00CDCD")
#text(x=log(10),y=-0.030,"Cold region",adj=c(0,1), col="#00CDCD", cex=0.7)
abline(h=0, col="white", lwd=3)
abline(h=-0.048, col="white", lwd=2)
### t-test for the regional values
t_test_region <- 
tt<-t.test(log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[1])] + 0.00000001), log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[1])] + 0.00000001))
kk <- aov(formula= cli_vel_tmp_hol$V2 ~ log(cli_vel_tmp_hol$V1 + 0.00000001))

for_t_test <- na.omit(cli_vel_tmp_hol)
for_t_test[,1] <- log(for_t_test$V1 + 0.00000001)
colnames(for_t_test) <- c("Clim_vel", "Time")
lm_region <- lm(for_t_test$Clim_vel ~ as.factor(for_t_test$Time))
anova_region <- anova(lm_region)
aov_region <- aov(for_t_test$Clim_vel ~ as.factor(for_t_test$Time))
tukey_region <- TukeyHSD(x=aov_region, 'as.factor(for_t_test$Time)', conf.level=0.95)
install.packages("agricolae")
library(agricolae)
HSD_region  <- HSD.test(aov_region, 'as.factor(for_t_test$Time)', group=T, alpha=)
HSD_region_mat <- as.data.frame(HSD_region$groups, stringsAsFactors=F)
HSD_region_mat[,1] <- as.numeric(levels(HSD_region_mat[,1]))
HSD_region_mat <- HSD_region_mat[order(HSD_region_mat[,1]),]



plot(HSD_region_mat[,3], HSD_region_mat[,1], pch=19)
tukey_region_df <- as.data.frame(tukey_region$'as.factor(for_t_test$Time)')

ajuste <- lm(chocolate$Sabor ~ chocolate$Tipo + chocolate$Provador)
summary(ajuste)
anova(ajuste)
a1 <- aov(chocolate$Sabor ~ chocolate$Tipo + chocolate$Provador)
posthoc <- TukeyHSD(x=a1, 'chocolate$Tipo', conf.level=0.95)

plot(aov_region)

### histogram for the regional values ENDS here ####
bins_col_func <- colorRampPalette(c("#8B4500", "#8B6508", "#FF7F00","#FFA500","#698B22", "#90EE90", "#B4EEB4"))
bins_col <- bins_col_func(28)
hist_breaks <- seq(log(0.000000001),log(2000), 1)
hist_matrix <- matrix(ncol=37, nrow=length(hist_breaks)-1)
for (point in seq_along(bins_50to0)){
    for_hist_matrix <- hist(log(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == bins_50to0[point])] + 0.00000001), breaks=hist_breaks, xaxt="n",yaxt="n",col=paste(bins_col[point], 90, sep=""), border="white", lwd=0.1, ylim=c(0,0.5),xlim=log(c(0.00001, 1000)), freq=F, yaxs="i", xaxs="i", main=NULL, xlab=NULL, ylab=NULL)
    hist_matrix[,point] <- for_hist_matrix$density         
  }
x <- hist_breaks[-length(hist_breaks)]
y <- bins_50to0
z <- hist_matrix
col_3D <- bins_col
col_3D <- terrain.colors(28)
col_var <- x
for (t in seq_along(x)[-1]){
  col_3D  <- rbind(col_3D , terrain.colors(28))
}
for (t in 1:36){
  col_var <- cbind(col_var , x)
}
hist3D(x=x, y=y, z=z, colvar=y, col=as.matrix(col_3D),phi=0, theta=40, border="white", lwd=0.5, xlab="Climate velocity", ylab="Time", labels=x, zlab="", colkey=as.list(bins_col))
########### Ends here ########

dim(col_3D)
factor(z, levels=seq(19, 99, by=10))
dim(col_var)
dim(z)






### using median ###
g_median_species <- matrix(nrow=37, ncol=17)
rownames(g_median_species) <- bins_50to0
colnames(g_median_species) <- Single_sp$Species
for (species in seq_along(colnames(g_median_species))){
  for (bin in seq_along(rownames(g_median_species))){
    g_median_species[bin,species] <- median(na.omit(species_set_clim$Cli_vel_prc[species_set_clim$Species == colnames(g_median_species)[species] & species_set_clim$Time_bin == rownames(g_median_species)[bin]]))
  }
}

### Plot the median per bin per species, using different colors for herbivorous, carnivorous, and rodents####
### working !!!! ####
for (sp in seq_along(colnames(g_median_species))){
  if (sp == 1){
    plot(as.numeric(rownames(g_median_species)), g_median_species[,sp], pch=19, cex=1.3, col=paste(Single_sp$Sp_color[which(Single_sp$Species == colnames(g_median_species)[sp])], "99", sep=""), frame=F, xlab="Time (years BP)", ylab=("Climate velocity (km/year)"), ylim=c(0, 0.5))
  } 
  else{
    points(as.numeric(rownames(g_median_species)), g_median_species[,sp], pch=19,cex=1.3, col=paste(Single_sp$Sp_color[which(Single_sp$Species == colnames(g_mean_species)[sp])], "99", sep=""))
  }
}
### add points estimating the climate velocity for the full region per bin
setwd("/Users/afr/Desktop/Regression/Cli_vel/")
cli_vel_Tmp <- stack("climate_change_velocity_perYear_tmp_rescaled_truncated_0_000005.grd")
e <- extent(-180,180,30, 90)
cropped <- crop(cli_vel_Tmp, e)
cli_vel_tmp_hol <- as.data.frame(matrix(ncol=3, nrow=0))
for (layer in 25:61){
  cli_vel_tmp_hol<- rbind(cli_vel_tmp_hol, cbind(extract(cropped[[layer]], e), rep(bins_50to0[layer-24], times=21600)))
}
### mean per bin
g_mean_region <- matrix(nrow=37, ncol=1)
rownames(g_mean_region) <- bins_50to0
for (bin in seq_along(rownames(g_mean_region))){
  g_mean_region[bin,1] <- exp(mean(log(as.vector(na.omit(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == rownames(g_mean_region)[bin])] +0.000000000000001)))))
}

kk <- log(as.vector(na.omit(cli_vel_tmp_hol$V1[which(cli_vel_tmp_hol$V2 == rownames(g_mean_region)[bin])])))
plot(g_mean_region[,1], rownames(g_mean_region))
### median per bin
g_median_region <- matrix(nrow=37, ncol=1)
rownames(g_median_region) <- bins_50to0
for (bin in seq_along(rownames(g_median_region))){
  g_median_region[bin,1] <- median(na.omit(cli_vel_tmp_hol$V1[cli_vel_tmp_hol$V2 == rownames(g_median_region)[bin]]))
}
points(rownames(g_median_region)[37:1], g_median_region[37:1,1], col="yellow", pch=17)














##### Estimate and plot climate change based on greenlandic database
clim <- read.delim("~/Desktop/PhD/Thesis/1stChapter/clim_greenland", header=T, stringsAsFactors=F)
#### Extract the first 50000 years >>>>  Starts here####
clim_temp <- clim[clim$Years <= 50000,]
clim_50 <- clim_temp[c(0,seq(1,2002,2)),]
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i",xaxt="n", yaxt="n", ylab="",xlab="",frame.plot=F,xlim=c(0, 51000),ylim=c(-46, -10) )
for (bin in seq_along(bins_50to0)){
  plot(clim_50$Years == bins_50to0)
}
temp_match <- match(bins_50to0, clim_50$Years)
bins_green <- cbind.data.frame(clim_50$Years[temp_match], clim_50$Del_18O[temp_match])
diff_green <- matrix(ncol=2, nrow=length(bins_green[,1])-1)
diff_green[,1] <- bins_green[,1][-37]
for (bin in seq_along(diff_green[,1])){
  diff_green[bin,2] <- bins_green[bin,2] - bins_green[bin+1,2]
}
plot(diff_green, type="l")
#### Extract the first 50000 years <<<<  Ends here ####



#### Extract the first 50000 years >>>>  Starts here####
clim_temp <- clim[clim$Years <= 50000,]
clim_50 <- clim_temp[c(0,seq(1,2002,2)),]
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i",xaxt="n", yaxt="n", ylab="",xlab="",frame.plot=F,xlim=c(0, 51000),ylim=c(-46, -10) )
for (bin in seq_along(bins_50to0)){
  plot(clim_50$Years == bins_50to0)
}
temp_match <- match(bins_50to0, clim_50$Years)
bins_green <- cbind.data.frame(clim_50$Years[temp_match], clim_50$Del_18O[temp_match])
diff_green <- matrix(ncol=2, nrow=length(bins_green[,1])-1)
diff_green[,1] <- bins_green[,1][-37]
for (bin in seq_along(diff_green[,1])){
  diff_green[bin,2] <- bins_green[bin,2] - bins_green[bin+1,2]
}
plot(diff_green, type="l")
#### Extract the first 50000 years <<<<  Ends here ####
clim_temp <- clim[clim$Years <= 51000,]
clim_50 <- clim_temp[c(0,seq(1,2042,2)),]
green_mean <- matrix(ncol=2, nrow=length(bins_green[,1])-1)
green_mean[,1] <- bins_green[,1][-37]
for (bin in seq_along(green_mean[,1])){
  green_mean[bin,2] <- median(clim_50$Del_18O[clim_50$Years <= green_mean[bin,1]+1000 & clim_50$Years > green_mean[bin,1]-1000 ])
}
diff_green_mean <- matrix(ncol=2, nrow=length(bins_green[,1])-1)
diff_green_mean[,1] <- green_mean[,1][-37]
for (bin in seq_along(diff_green_mean[,1])){
  diff_green_mean[bin,2] <- green_mean[bin,2] - green_mean[bin+1,2]
}
plot(abs(diff_green_mean), type="l")













cor_clim_spp <- matrix(ncol=17, nrow=17)
colnames(cor_clim_spp) <- Single_sp$Species
rownames(cor_clim_spp) <- Single_sp$Species
for (r_spp in seq_along(rownames(cor_clim_spp))){
  for (c_spp in seq_along(rownames(cor_clim_spp))){
    cor_clim_spp[r_spp,c_spp] <-  cor(na.omit(cbind.data.frame(g_mean_species[,r_spp], g_mean_species[,c_spp])))[1,2]
}
par(mar=c(7,7,3,3))
image(cor_clim_spp^2, col=gray.colors(10), axes=F)
axis(1, at=(seq(0,1,by=1/16)), labels=colnames(cor_clim_spp), las=2)
axis(2, at=(seq(0,1,by=1/16)), labels=rownames(cor_clim_spp), las=2)




glm(g_mean_species[,sp] ~ as.numeric(rownames(g_mean_species)))
######### TODO, plot an image for the geometric mean. Not use the ranking but other way to select (color) the slow and fast bins
Herb_spp <- c("Mammuthus_primigenius",  "Coelodonta_antiquitatis", "Saiga_tatarica", "Bison_sp", "Ovibos_moschatus", "Rangifer_tarandus", "Cervus_elaphus", "Equus_caballus")
Carn_spp <- c("Alopex_lagopus", "Canis_lupus" , "Crocuta_crocuta","Panthera_leo", "Ursus_arctos", "Ursus_spelaeus")
Rodn_spp <- c("Castor_fiber", "Dicrostonyx_torquatus", "Microtus_gregalis")
####### geometric mean #####


image(g_mean_species[37:1,], axes=F, col=terrain.colors(33))
axis(1, at=(seq(0,1,by=1/36)), labels=rownames(mean_species[37:1,]), las=2)
axis(2, at=(seq(0,1,by=1/16)), labels=colnames(mean_species), las=2)
for (x in 1:ncol(g_mean_species[37:1,])){
  counter <- 0
  #for (y in 1:nrow(mean_species[37:1,])){
    #extemes <- c(max(mean_species[37:1,][,x], na.rm=T), min(mean_species[37:1,][,x], na.rm=T))
    text(seq(0,1,by=1/36), (seq(0,1,by=1/16))[x] , rank(g_mean_species[37:1,x], na.last="keep"), cex=0.5)
    #text((x-1)/36, ((y-1)/16)-0.01, mean_species[37:1,][x,y], cex=0.5)
    counter <- counter + 1
  }




rank(mean_species[,1])
hist(mean_species[,1])

hist(species_set_clim$Cli_vel_prc[species_set_clim$Species == colnames(mean_species)[species] & species_set_clim$Time_bin == rownames(mean_species)[bin]])
bp_temp <- boxplot(species_set_clim$Cli_vel_prc[species_set_clim$Species == colnames(mean_species)[species] & species_set_clim$Time_bin == rownames(mean_species)[bin]])
points(as.vector(bp_temp$out), col="red")


length(bp_temp$out)
length()

#########################################################################################################################
# Climate velocity
#########################################################################################################################
Full_DB <- read.delim(file.choose(), header=T, stringsAsFactors=F)
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
Full_DB_LL <- Full_DB_LL[-which(Full_DB_LL$Median_Age > 50000),]
# add two empty colums to assign color and type of point in the map
For_map <- matrix(NA,nrow=length(Full_DB_LL$Latitude), ncol=2)
colnames(For_map) <- c("Map_color", "Map_type")
Full_DB_map <- cbind(Full_DB_LL, For_map)





bins<-c(seq(0, 21000, by=1000), seq(22000, 50000, by=2000))
bins_image <-c(seq(-500, 21500, by=1000), seq(23000, 49000, by=2000))
bins_length<- c(rep(1000, ))
library(raster)
# Upload the climate velocity files
setwd("/Users/afr/Desktop/Regression/Cli_vel/")
dir()
velocity_map_year <- stack("climate_change_velocity_perYear_prec_rescaled_truncated_0_000005.grd")
# Define the directory where the files
setwd("/Users/afr/Desktop/Evolution_ppt/Corr_result/")
image_climate <- matrix(NA,ncol=42, nrow=21)
colnames(image_climate) <- c("Q0", "Q25", "Q50", "Q75", "Q100" , seq(0, 21000, by=1000), seq(22000, 50000, by=2000))
rownames(image_climate) <- Single_sp
image_climate_rank <- image_climate
image_climate_n <- image_climate
#Single_sp <- unique(Full_DB_LL$Species)
Single_sp <- c("Alopex_lagopus","Bison_sp","Canis_lupus","Castor_fiber","Cervus_elaphus","Coelodonta_antiquitatis","Crocuta_crocuta","Dicrostonyx_torquatus","Equus_caballus","Mammuthus_primigenius","Microtus_gregalis","Ovibos_moschatus","Panthera_leo","Rangifer_tarandus","Saiga_tatarica","Ursus_arctos","Ursus_spelaeus")
for (s in seq_along(Single_sp)){
  temp_DB_climate <- Full_DB_LL[which(Full_DB_LL$Species==Single_sp[s]),]
  temp_points <- as.data.frame(matrix(nrow=nrow(temp_DB_climate), ncol=7))
  colnames(temp_points) <- c("Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "cell","Velocity")
  matrix_time <- matrix("numeric", nrow=38, ncol=2 )
  matrix_time[,1] <- c(seq(50000, 22000, by=-2000), seq(21000, -1000, by=-1000))
  matrix_time[,2] <- seq(25, 62, by=1)
  for(i in seq_along(temp_DB_climate$Median_Age)){
    for (j in 1:(dim(matrix_time)[1]-1)){
      if(temp_DB_climate$Median_Age[i]<= as.numeric(matrix_time[j]) & temp_DB_climate$Median_Age[i]> as.numeric(matrix_time[j+1])){
        temp_points[i,c(1,2,3,4,5)] <- c(temp_DB_climate$Longitude[i], temp_DB_climate$Latitude[i],temp_DB_climate$Median_Age[i], matrix_time[j,1], matrix_time[j,2] )
      }
    }
  }
  for(k in seq_along(temp_points$Longitude)){
    temp_points[k,c(6,7)]<- extract(velocity_map_year, layer=as.numeric(temp_points[k,5]), nl=1, y=matrix(as.numeric(temp_points[k,c(1,2)]), nrow=1, ncol=2), cellnumbers=T)
  }
  points_plot <- temp_points[!is.na(temp_points$Velocity),]
  bp_all_points <- boxplot(points_plot$Velocity, plot=F)
  #boxplot(points_plot$Velocity)
  bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), plot=F)
  tick <- as.vector(as.numeric(bp$names))
  # Check point
  sort(as.numeric(points_plot$Time_bin))
  # Save the velocity values for every species
  ### setwd("/Users/afr/Desktop/kk_temp/")
  ### write.table(x=temp_points, file=paste(tolower(Single_sp[s]), "_vel_bin.txt", sep=""), sep="\t",row.names=F)
  species_sp <- paste(strsplit(strsplit(tolower(Single_sp[s]), split = "_")[[1]][1], split = "")[[1]][1], strsplit(strsplit(tolower(Single_sp[s]), split = "_")[[1]][2], split = "")[[1]][1], sep="")
  bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), main=Single_sp[s], border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50))
  tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
  axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
  axis(side=4, at=seq(0, 1, by=0.2), cex=0.5)
  #lines(as.numeric(bp$names)/1000, bp$stats[1,], col="red", lwd=2)
  #lines(as.numeric(bp$names)/1000, bp$stats[2,], col="pink", lwd=2)
  #lines(as.numeric(bp$names)/1000, bp$stats[3,], col="grey", lwd=2)
  #lines(as.numeric(bp$names)/1000, bp$stats[4,], col="lightblue", lwd=2)
  #lines(as.numeric(bp$names)/1000, bp$stats[5,], col="green", lwd=2)
  image_climate[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- bp$stats[3,]
  image_climate[s, c(1,2,3,4,5)] <-as.vector(bp_all_points$stats)
  image_climate_rank[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- rank(bp$stats[3,])
  image_climate_n[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- bp$n                                                                                    
}

ln<- layout(as.matrix(cbind(c(1,2), c(4,3))), heights=c(20,80),widths=(c(80,20)))
par(mar=c(1,13,4,2))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylab="",frame.plot=F,ylim=c(-46, -34), main="Climate velocity per bin (color relative to ranking in the full database)")
axis(2, at=(c(-46, -34)), labels=c(-46, -34))
mtext("Temperature", side=2, line=2)
par(mar=c(5,13,0,2))
image(t(image_climate_rank[,6:42]), axes=F)
axis(1, at=(seq(0,1,by=1/36)), labels=colnames(image_climate_rank[,6:42]), las=2)
axis(2, at=(seq(0,1,by=1/20)), labels=rownames(image_climate_rank), las=2)
axis(4, at=(seq(0,1,by=1/20)), labels=rowsum(image_climate_rank, ), las=2)
for (x in 1:ncol(image_climate_rank)){
  counter <- 0
  for (y in 1:nrow(image_climate_rank)){
    extemes<- c(max(t(image_climate_rank[,6:42])[,y], na.rm=T), min(t(image_climate_rank[,6:42])[,y], na.rm=T))
    text((x-1)/36, ((y-1)/20)+0.01, t(image_climate_rank[,6:42])[x,y], cex=ifelse(is.na(sum(match(t(image_climate_rank[,6:42])[x,y],extemes))),0.3, 0.9))
    text((x-1)/36, ((y-1)/20)-0.01, t(image_climate_n[,6:42])[x,y], cex=0.5)
    counter <- counter + 1
    print(sum(match(t(image_climate_rank)[x,y],extemes)))
  }
}
par(mar=c(0,0,0,0))
plot.new()
legend.gradient(cbind(c(0,0.5,0.5,0), c(1,1,0.5,0.5)), cols=paste(colores, 99, sep=""), limits=c("Slowest", "Fastest"), title="Ranking")

plot(image_climate[1,])
boxplot(mean_species[,])
  