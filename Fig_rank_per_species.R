# Reset
rm(list=ls())

# Files needed
setwd("/Users/afr/Desktop/Regression/Cli_vel/")
cli_vel_Tmp <- stack("climate_change_velocity_perYear_tmp_rescaled_truncated_0_000005.grd")
cli_vel_Prc <- stack("climate_change_velocity_perYear_prec_rescaled_truncated_0_000005.grd")
setwd("/Users/afr/Desktop/Figures_Ditte")
Single_sp <- read.delim("Sp_information.txt", header=T, stringsAsFactors=F)

# Import the DATABASE
setwd("/Users/afr/Desktop/Figures_Ditte")
Full_DB <- read.delim("DATABASE_for_ditte.txt", header=T, stringsAsFactors=F)
Full_DB$Longitude <- as.numeric(Full_DB$Longitude)
Full_DB$Latitude <- as.numeric(Full_DB$Latitude)
full_vector <- as.vector(NULL)
for (i in seq_along(Full_DB$Latitude)){
  if (is.na(Full_DB$Longitude[i]) || is.na(Full_DB$Latitude[i])){
    full_vector <- c(full_vector,i)
  }
}
DATABASE <- Full_DB[-full_vector,]
DATABASE <- DATABASE[-which(DATABASE$Mean_Age >= 50000),]

# objects needed 
bins_0to50 <- c(seq(0,21000, by=1000), seq(22000, 50000, by=2000))
bins_50to0 <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
bins_mid_point <- c(seq(500, 21500, by=1000), seq(23000, 49000, by=2000))
layers <- seq(25, 61, by=1)
setwd("/Users/afr/Desktop/Evolution_ppt/Corr_result/")
image_cli_vel_tmp <- matrix(NA,ncol=41, nrow=length(Single_sp$Species))
colnames(image_cli_vel_tmp) <- c("Q0", "Q25", "Q50", "Q75", "Q100" , seq(1000, 21000, by=1000), seq(22000, 50000, by=2000))
rownames(image_cli_vel_tmp) <- Single_sp$Species
image_climate_rank <- image_cli_vel_tmp[,6:41]
colores <- colorRampPalette(c("#E5F5F9", "#CCECE6", "#99D8C9", "#66C2A4", "#41AE76", "#238B45", "#005824","#00441B"))(6)
image_quartiles <- matrix(NA,ncol=36, nrow=length(Single_sp$Species))

# Script body 
Box_bin <- FALSE
Rank_species <- TRUE
for (s in seq_along(Single_sp$Species)){
  # Database to work with, contains all the species and data
  temp_DB_climate <- DATABASE[which(DATABASE$Species==Single_sp$Species[s]),]
  # Empty data frame to save the data for every species
  temp_colnames <- c("Species","Longitude", "Latitude", "Mean_date_record", "Time_bin", "Layer", "Rec_type", "Cell", "Cli_vel_tmp", "Cli_vel_prc")
  temp_points <- as.data.frame(matrix(nrow=nrow(temp_DB_climate), ncol=length(temp_colnames)))
  colnames(temp_points) <- temp_colnames
  counter1 <- 1
  for (bin in seq_along(bins_50to0)[-37]){
    records <- which(temp_DB_climate$Mean_Age < bins_50to0[bin] & temp_DB_climate$Mean_Age >= bins_50to0[bin+1])
    points <- as.matrix(cbind.data.frame(as.numeric(temp_DB_climate$Longitude[records]), as.numeric(temp_DB_climate$Latitude[records])))
    colnames(points) <- c("longitude", "Latitude")
    if(length(records > 0)){
      temp_points[(counter1:(counter1 + length(records)-1)),1] <- Single_sp$Species[s]
      temp_points[(counter1:(counter1 + length(records)-1)),2:7] <- cbind.data.frame(as.numeric(temp_DB_climate$Longitude[records]), as.numeric(temp_DB_climate$Latitude[records]),as.numeric(temp_DB_climate$Mean_Age[records]), rep(bins_50to0[bin], times=length(records)), rep(layers[bin], times=length(records)),ifelse(nchar(temp_DB_climate$Sequence[records]) > 1, "Seq", "Fossil") )
      temp_points[(counter1:(counter1 + length(records)-1)),c(which(colnames(temp_points) == "Cell"), which(colnames(temp_points) == "Cli_vel_tmp")) ] <- extract(cli_vel_Tmp, layer=layers[bin], nl=1, y=points, cellnumbers=T)
      temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == "Cli_vel_prc")]<- extract(cli_vel_Prc, layer=layers[bin], nl=1, y=points)
      counter1 <- counter1+length(records)
    }
  }
  points_plot <- temp_points[!is.na(temp_points$Cli_vel_tmp),]
  bp_all_points <- boxplot(points_plot$Cli_vel_tmp, plot=F)
  #boxplot(points_plot$Velocity)
  bp <- boxplot(points_plot$Cli_vel_tmp ~ as.numeric(points_plot$Time_bin), plot=F)
  tick <- as.vector(as.numeric(bp$names))
  # Check point
  sort(as.numeric(points_plot$Time_bin))
  # Save the velocity values for every species
  ### setwd("/Users/afr/Desktop/kk_temp/")
  ### write.table(x=temp_points, file=paste(tolower(Single_sp[s]), "_vel_bin.txt", sep=""), sep="\t",row.names=F)
  species_sp <- tolower(Single_sp$Sp[s])
  if (Box_bin == TRUE){
  maxs <- c(max(bp$stats[5,]), max(bp$out))
  bp <- boxplot(points_plot$Cli_vel_tmp ~ as.numeric(points_plot$Time_bin), main=Single_sp$Species[s], border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50), ylim=c(0, max(maxs)), xaxs="i", yaxs="i")
  tick_labels <- bins_0to50 
  axis(side=1, at=tick_labels/1000, labels=tick_labels/1000, cex=0.5)
  axis(side=2, cex=0.5)
  }
  image_cli_vel_tmp[s,as.vector(na.omit(match(bp$names, colnames(image_cli_vel_tmp))))] <- bp$stats[3,]
  image_cli_vel_tmp[s, c(1,2,3,4,5)] <-as.vector(bp_all_points$stats)
  image_climate_rank[s,as.vector(na.omit(match(bp$names, colnames(image_climate_rank))))] <- rank(bp$stats[3,])
}
if(sort_sp == TRUE){
  image_cli_vel_tmp <- image_cli_vel_tmp[order(Single_sp$Generation_time), ]
  image_climate_rank <- image_climate_rank[order(Single_sp$Generation_time), ]
}
# Script fot plotting
if(Rank_species == TRUE){
ima_row <- length(Single_sp$Species)
ima_col <- length(bins_mid_point)
X0 <- matrix(ncol=ima_col, nrow=ima_row)
for (x in 1:ima_row){
  X0[x,] <- bins_mid_point
}
Y0<- matrix(ncol=ima_col, nrow=ima_row)
for (y in 1:ima_col){
  Y0[,y] <- seq(1,ima_row, by=1)
}

for (q in 1:dim(image_cli_vel_tmp[,6:41])[1]){
  for (r in 1:dim(image_cli_vel_tmp[,6:41])[2]){
    if (is.na(image_cli_vel_tmp[q,r+5])){
      image_quartiles[q,r] <- NA
    }
    else if (image_cli_vel_tmp[q,r+5] < image_cli_vel_tmp[q,2]){
      image_quartiles[q,r] <- 25
    }
    else if (image_cli_vel_tmp[q,r+5] >= image_cli_vel_tmp[q,2] & image_cli_vel_tmp[q,r+5] < image_cli_vel_tmp[q,3]){
      image_quartiles[q,r] <- 50
    }
    else if (image_cli_vel_tmp[q,r+5] >= image_cli_vel_tmp[q,3] & image_cli_vel_tmp[q,r+5] < image_cli_vel_tmp[q,4]){
      image_quartiles[q,r] <- 75
    }
    else if (image_cli_vel_tmp[q,r+5] >= image_cli_vel_tmp[q,4]){
      image_quartiles[q,r] <- 100
    }
  }
}
ln<- layout(matrix(c(1,2), ncol=2), heights=1,widths=(c(80,20)))
layout.show(ln)
par(mar=c(5,5,7,0))
poly.image(X0, Y0, z=image_quartiles,col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", yaxt="n",xaxt="n", xlim=c(0, 50000), border="white", frame=F)
axis(1, at = bins_0to50, labels = bins_0to50/1000)
axis(2, at=Y0[,1], labels=rownames(image_cli_vel_tmp), las=2)
mtext("ranking per species", side=3,line=1, cex=2)
# add the ranking
for (x in 1:ncol(image_climate_rank)){
  counter <- 0
  for (y in 1:nrow(image_climate_rank)){
    extremes <- c(max(t(image_climate_rank)[,y], na.rm=T), min(t(image_climate_rank)[,y], na.rm=T))
    text(bins_mid_point[x], y, t(image_climate_rank)[x,y], cex=ifelse(is.na(sum(match(t(image_climate_rank)[x,y],extremes))),0.35, 0.5))
    counter <- counter + 1
    print(sum(match(t(image_climate_rank)[x,y],extremes)))
  }
}
# add the legend
par(mar=c(5,3,7,3))
plot(1,1, col="white", frame.plot=F, axes=F, xlim=c(0,2), ylab=NA,xlab=NA)
legend.gradient(cbind(c(0,1,1,0), c(1.375,1.375,0.593,0.593)), cols=paste(colores, 99, sep=""),limits=c("",""), title="")
text(1.80,1.344, labels="Fast")
text(1.80,0.625, labels="Slow", srt=90)
}

# plotting for all the database values
ln<- layout(as.matrix(rbind(c(1,2),c(3,4))), heights=c(70,30),widths=(c(80,20)))
layout.show(ln)

#if(Rank_all == TRUE){
  ima_row <- length(Single_sp$Species)
  ima_col <- length(bins_mid_point)
  X0 <- matrix(ncol=ima_col, nrow=ima_row)
  for (x in 1:ima_row){
    X0[x,] <- bins_mid_point
  }
  Y0<- matrix(ncol=ima_col, nrow=ima_row)
  for (y in 1:ima_col){
    Y0[,y] <- seq(1,ima_row, by=1)
  }
  par(mar=c(0,5,5,0))
  poly.image(X0, Y0, z=image_cli_vel_tmp[, 6:41],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", yaxt="n",xaxt="n", xlim=c(0, 50000), border="white", frame=F)
  axis(1, at = bins_0to50, labels = NA)
  axis(2, at=Y0[,1], labels=rownames(image_cli_vel_tmp), las=2)
  mtext("ranking full database", side=3,line=1, cex=2)
  # add the legend
  par(mar=c(1,0,5,3))
  plot(1,1, col="white", frame.plot=F, axes=F, xlim=c(0,2), ylab=NA,xlab=NA)
  legend.gradient(cbind(c(0,1,1,0), c(1.375,1.375,0.593,0.593)), cols=paste(colores, 99, sep=""),limits=c("",""), title="")
  text(1.80,1.344, labels="Fast")
  text(1.80,0.625, labels="Slow", srt=90)


clim <- read.delim("~/Desktop/PhD/Thesis/1stChapter/clim_greenland", header=T, stringsAsFactors=F)
#### Extract the first 50000 years
clim_temp <- clim[clim$Years <= 50000,]
clim_50 <- clim_temp[c(0,seq(1,2000,2)),]
par(mar=c(3,5,0,0))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i",xaxt="n", yaxt="n", ylab="",xlab="",frame.plot=F,xlim=c(0, 50000),ylim=c(-46, -32))
axis(1, at=bins, labels=bins/1000)
axis(2, at=(c(-46, -34)), labels=c(-46, -34), cex=0.7, las=2)
mtext("Temperature", side=2, line=2,at=-40)





#### Plot image using all the values NO LAYOUT, add=T#####
ln<- layout(matrix(c(1,2), ncol=2), heights=1,widths=(c(80,20)))
layout.show(ln)
#plot.new()
par(mar=c(5,14,5,6))
mtext("ranking relative to all the records", side=3,line=0)
#par(mar=c(0,0,0,0))
#poly.image(add=T,X1, Y1, z=image_climate[,6:28],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(-500, 22000), border="white")
#par(mar=c(0,0,0,0))
#poly.image(add=T,X2, Y2, z=image_climate[,29:42],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(22000, 50000), border="white")
axis(1, at=bins,pos=-32.5, labels=NA)
#image(t(image_climate[,6:42]), col=paste(colores, 90, sep=""), axes=F,)
#axis(1, at=bins_image, labels=colnames(image_climate[,6:42]), las=2)
axis(2, at=Y0[,1], labels=rownames(image_climate), las=2)
par(mar=c(4.5,7,5,2))
#plot(1,1)
plot(1,1, col="white", frame.plot=F, axes=F, xlim=c(0,2), ylab=NA,xlab=NA)
legend.gradient(cbind(c(1.07,1.62,1.62,1.07), c(1.375,1.375,0.593,0.593)), cols=paste(colores, 99, sep=""),limits=c("",""), title="")
segments(1.07,0.593,1.62,0.593)
segments(1.07,0.593,1.07,1.375)
segments(1.62,0.593,1.62,1.375)
segments(1.62,1.375,1.07,1.375)
text(1.80,1.344, labels="Fast", srt=90)
text(1.80,0.625, labels="Slow", srt=90)
#text("Slowest", x=(1.44-0.55), y=0.7)
#text("Fastest", x=(1.44-0.55), y=1.2)


