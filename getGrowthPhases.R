################################################################
# Name: getGrowthPhases.R
# To-do: Identifies different growth phases from fermentation data
# Date: 14-07-2017
# Author: Jens Christian Nielsen, jens.c.nielsen@chalmers.se | jcfnielsen@gmail.com
################################################################

###############
# Import data #
###############
# Specify input data file
time_file <- "time.csv"
accumulated_CO2_file <- "accumulated_CO2.csv" 
exhaust_CO2_file <- "exhaust_CO2.csv"

# read input data
time_mat <- as.matrix( read.csv(file = time_file, header = T) )
accu_CO2_mat <- as.matrix( read.csv(file = accumulated_CO2_file, header = T) )
exhust_CO2_mat <- as.matrix( read.csv(file = exhaust_CO2_file, header = T) )
#time_mat <- time_mat[,49:54]

accu_CO2_mat[accu_CO2_mat<0] = 0
exhust_CO2_mat[is.na(exhust_CO2_mat)] <- 0.04


# Save plots as pdf
pdf("growth_phases.pdf",width = 12,height = 8)
par(mfrow=c(2,3))
par(mar = c(5,5,2,5))

######################
# Specify parameters #
######################
rsq_cutoff <- 0.99 # r-squared cut-off
points = seq(50,200,2) # Number of points for the linear fit
max_percent_change <- 0.1

####################################
# Calculate growth characteristics #
####################################
# Initialization
data <- matrix(NA,nrow = dim(time_mat)[2],ncol=6)
colnames(data) <- c("mu_max","R-squared","datapoints","exp_end","max_CO2","data")
rownames(data) <- colnames(time_mat)

for( fermentation in colnames(time_mat) ) {
  # Get data for each fermentation
  toUse <- which( !is.na(time_mat[,fermentation]) & !is.na(accu_CO2_mat[,fermentation]) )
  accu_CO2 <- log(accu_CO2_mat[,fermentation][toUse])
  exhaust_CO2 <-  exhust_CO2_mat[,fermentation][toUse]
  accu_CO2[is.infinite(accu_CO2)] = unique(accu_CO2[order(accu_CO2,decreasing = F)])[2]
  time <- time_mat[,fermentation][toUse]
  
  
  # re-initialization for each run
  max_slope = 0.001
  end = 0
  
  for(point in points) {
    
    # loop over all possible slopes
    for(i in 1:(length(accu_CO2)-point)) {

      min_point <- i
      max_point <- i + point
      
      # Get linear fit
      fit <- lm( accu_CO2[min_point:max_point] ~ time[min_point:max_point] )
      intercept <- coef(fit)[1]
      slope <- coef(fit)[2]
      rsq <- summary(fit)$r.squared
      
      
      if (is.na(rsq)) {next}
      percent_change <- (max_slope-slope)/max_slope
      
      # save slope as mslope if:
      # - slope is new max or no more than 10 % samller than max
      # - r_squared (rsq) is more than the cut-off
      # - 
      #if (rsq > rsq_cutoff && time[min_point] > 10 && percent_change < max_percent_change && max_point > end) {
      if (rsq > rsq_cutoff && percent_change < max_percent_change && max_point > end) {
        best_fit <- fit
        if (slope > max_slope) {max_slope <- slope}
        max_intercept <- intercept
        max_rsq <- rsq
        start <- min_point
        end <- max_point
        data_points_used <- max_point-min_point
      }
    }
  }
  
  # Plot best fit  
  plot(time,accu_CO2,xlab="time (h)", ylab = expression("log(accumulated CO"[2]*")"),cex.lab=1.2,cex.axis=1.2,cex.main=1.2,main=fermentation,col="black",pch=19)
  points(time[start:end],accu_CO2[start:end],pch=19,col="blue")
  abline(best_fit,lwd=2,col="red")
  par(new=T)
  plot(time,exhaust_CO2,col="grey",axes=F, xlab=NA, ylab=NA,pch=19)
  axis(side = 4,cex.axis=1.2)
  mtext(side = 4, line = 3, expression("CO"[2]*" exhaust (%)"),cex.lab=0.7)
  
  index_max_co2 <- which(exhaust_CO2==max(exhaust_CO2))[1] 
  rect(time[end],-1,time[index_max_co2],10,col = rgb(0.5,0.5,0.5,1/4),border = NA)

  legend("bottomright", bty="n", cex=1.2,
         legend=paste0("R-squared: ",format(max_rsq, digits=4),"\n","Max CO2 rate: ",format(max_slope,digits = 4),"\n", "No. data points: ", data_points_used,"\n","End exp phase: ",format(time[end],digits=0), "\n","CO2 max: ",format(time[index_max_co2],digits=0),"\n\n\n\n")
         )
  
  # Save summary data of each fermentation
  data[fermentation,"mu_max"] <- format(max_slope,digits = 4)
  data[fermentation,"R-squared"] <- format(max_rsq, digits=4)
  data[fermentation,"datapoints"] <- data_points_used
  data[fermentation,"exp_end"] <- as.numeric(format(time[end],digits=0))
  data[fermentation,"max_CO2"] <- as.numeric(format(time[index_max_co2],digits=0))
  data[fermentation,"data"] <- as.numeric(time[index_max_co2]-time[end])
  colnames(data) <- c("mu_max","R-squared","datapoints","exp_end","max_CO2","data")
}
dev.off()
