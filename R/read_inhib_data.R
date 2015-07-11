##' Reads in data frames, corrects bad data, and returns a master list of the data
##' 
##' @param printDataCheck prints plots of the data, to show what bad data points were thrown out.
##' @return a list containing each data set
##' @export
##' 

read_inhib_data <- function(printDataCheck=FALSE) {
  ########
  # Function (sort of) to read in all the data from the substrate inhibition project
  ########
  
  ###################
  # Arg-AMC, TN river
  ###################
  # Read raw data
  #rawFl <- read.csv("2013-03-07Arg-AMCcompetitioncurve.csv")
  #X2013_03_07Arg_AMCcompetitioncurve <- rawFl
  #save(X2013_03_07Arg_AMCcompetitioncurve, file="data/X2013_03_07Arg_AMCcompetitioncurve.Rdata")
  load("data/X2013_03_07Arg_AMCcompetitioncurve.Rdata")
  rawFl <- X2013_03_07Arg_AMCcompetitioncurve
  rawFl$pNaSub <- as.character(rawFl$pNaSub)
  rawFl <- na.omit(rawFl) # There were a few spots where I didn't take data
  attr(rawFl$time, "units") <- "days"
  rawFl$incTime <- (rawFl$time-min(rawFl$time))*24 
  attr(rawFl$incTime, "units") <- "hours"
  
  # Check for bad data points
  argTNdataCheckPre <- ggplot(rawFl, aes(x=incTime, y=fl)) +
    geom_point() + geom_line() +
    facet_grid(pNaSub~conc.uM) +
    scale_x_continuous("incubation time, hrs") +
    ggtitle("unedited arg-AMC TN")
  if(printDataCheck) {
    print(argTNdataCheckPre)
  }
  
  # Remove the two bad killed control points. All the other points look good, qualitatively
  rawFl <- rawFl[-which(rawFl$pNaSub=="killed" & rawFl$conc.uM==0 & rawFl$fl > 2000), ]
  rawFl <- rawFl[-which(rawFl$pNaSub=="killed" & rawFl$conc.uM==280 & rawFl$fl < 3000), ]
  rawFl <- rawFl[-which(rawFl$pNaSub==10 & rawFl$conc.uM==320 & rawFl$fl==3798), ]
  
  # Confirm that bad data points are gone
  argTNdataCheckPost <- ggplot(rawFl, aes(x=incTime, y=fl)) +
    geom_point() + geom_line() +
    facet_grid(pNaSub~conc.uM) +
    scale_x_continuous("incubation time, hrs") +
    ggtitle("edited arg-AMC TN")
  if(printDataCheck) {
    print(argTNdataCheckPost)
  }
  
  
  # Rename column names appropriately (and get rid of useless columns) & rename
  rawFl$AMC.subs <- "arg-AMC"
  argTNrawFl <- rename(rawFl, replace=c("pNaSub" = "pNA.subs", "conc.uM" = "conc.pNA"))
  argTNrawFl$pNA.subs <- as.factor(argTNrawFl$pNA.subs)
  
  # rename raw data, to ultimately put it into a list
  rm(rawFl)
  
  #############
  # Leu-AMC in TN river
  #############
  
  # Read raw data
  #rawFl <- read.csv("2013-02-10 leu-AMC competitive inhibition.csv")
  #X2013_02_10_leu_AMC_competitive_inhibition <- rawFl
  #save(X2013_02_10_leu_AMC_competitive_inhibition, file="data/X2013_02_10_leu_AMC_competitive_inhibition.Rdata")
  load("data/X2013_02_10_leu_AMC_competitive_inhibition.Rdata")
  rawFl <- X2013_02_10_leu_AMC_competitive_inhibition
  
  rawFl <- na.omit(rawFl) # There were a few spots where I didn't take data
  attr(rawFl$time, "units") <- "days"
  rawFl$incTime <- (rawFl$time-min(rawFl$time))*24 
  attr(rawFl$incTime, "units") <- "hours"
  
  # Check for bad data points
  leuTNdataCheckPre <- ggplot(rawFl, aes(x=incTime, y=fl, colour=liveOrDead)) +
    geom_point() + geom_line() +
    facet_grid(pNaSub~conc.uM) +
    scale_x_continuous("incubation time, hrs") +
    ggtitle("unedited leu-AMC TN River")
  if(printDataCheck) {
    print(leuTNdataCheckPre)
  }
  
  # Remove bad data points:
  rawFl <- rawFl[-which(rawFl$pNaSub == 1 & rawFl$conc.uM == 40 & rawFl$liveOrDead=="live" & rawFl$fl < 2000), ]
  rawFl <- rawFl[-which(rawFl$pNaSub == 7 & rawFl$conc.uM == 240 & rawFl$liveOrDead=="live" &  rawFl$time < 0.57), ]
  
  # Confirm that bad data points are gone
  leuTNdataCheckPost <- ggplot(rawFl, aes(x=incTime, y=fl, colour=liveOrDead)) +
    geom_point() + geom_line() +
    facet_grid(pNaSub~conc.uM) +
    ggtitle("edited leu-AMC TN River")
  if(printDataCheck) {
    print(leuTNdataCheckPost)
  }
  
  # Rename column names appropriately (and get rid of useless columns)
  rawFl$AMC.subs <- "leu-AMC"
  leuTNrawFl <- rename(rawFl, replace=c("pNaSub" = "pNA.subs", "conc.uM" = "conc.pNA"))
  
  # 'recast' (in the colloquial sense) rawFl so that 'killed' is a level of $pNA.subs
  leuTNkilledControl <- leuTNrawFl[leuTNrawFl$liveOrDead=="killed", ]
  leuTNkilledControl$pNA.subs="killed"
  leuTNrawFl <- rbind(leuTNrawFl[leuTNrawFl$liveOrDead=="live", ], leuTNkilledControl)
  leuTNrawFl$pNA.subs <- as.factor(leuTNrawFl$pNA.subs)
  
  # rename raw data, to ultimately put it into a list
  rm(rawFl)
  
  #############
  # Arg-AMC in Bogue Sound
  #############
  
  # Read in raw data
  #argRaw <- read.csv("2013-03-25 arg competitive inhibition.csv")
  #X2013_03_25_arg_competitive_inhibition <- argRaw
  #save(X2013_03_25_arg_competitive_inhibition, file="data/X2013_03_25_arg_competitive_inhibition.Rdata")
  load("data/X2013_03_25_arg_competitive_inhibition.Rdata")
  argRaw <- X2013_03_25_arg_competitive_inhibition
  
  # Change times to relative incubation time
  argRaw$Rtimes <- ymd_hms(paste("2013-03-25 ", as.character(argRaw$time), ":00", sep=""))
  argRaw$incTime <- as.numeric(argRaw$Rtimes - min(argRaw$Rtimes))/3600
  attr(argRaw$incTime, "units") <- "hours"
  
  # Plot the raw data
  argNCdataCheckPost <- ggplot(argRaw, aes(x=incTime, y=fl)) +
    geom_point() + geom_line() +
    facet_grid(pNA.subs~conc.pNA) +
    scale_x_continuous("incubation time, hrs") +
    ggtitle("unedited arg-AMC Bogue Sound")
  if(printDataCheck) {
    print(argNCdataCheckPost)
  }
  
  # Note that there don't seem to be any bad data points for arg-AMC in Bogue Sound
  
  
  # rename raw data, to ultimately put it into a list
  argNCrawFl <- argRaw
  rm(argRaw)
  
  #############
  # Leu-AMC in Bogue Sound
  #############
  
  # Read in raw data
  #leuRaw <- read.csv("2013-03-25 leu competitive inhibition.csv")
  #X2013_03_25_leu_competitive_inhibition <- leuRaw
  #save(X2013_03_25_leu_competitive_inhibition, file="data/X2013_03_25_leu_competitive_inhibition.Rdata")
  load(file="data/X2013_03_25_leu_competitive_inhibition.Rdata")
  leuRaw <- X2013_03_25_leu_competitive_inhibition
  
  # Change times to relative incubation time
  leuRaw$Rtimes <- ymd_hms(paste("2013-03-25 ", as.character(leuRaw$time), ":00", sep=""))
  leuRaw$incTime <- as.numeric(leuRaw$Rtimes - min(leuRaw$Rtimes))/3600
  attr(leuRaw$incTime, "units") <- "hours"
  
  # Plot the raw data
  leuNCdataCheckPre <- ggplot(leuRaw, aes(x=incTime, y=fl)) +
    geom_point() + geom_line() +
    facet_grid(pNA.subs~conc.pNA) +
    scale_x_continuous("incubation time, hrs") +
    ggtitle("unedited leu-AMC Bogue Sound")
  if(printDataCheck) {
    print(leuNCdataCheckPre)
  }
  
  # Remove bad data points 
  # pNA 1, conc 80 and 120 are all bad due to pipetting error
  leuEdit <- leuRaw[-which(leuRaw$conc.pNA==120 & leuRaw$pNA.subs == 1), ]
  leuEdit <- leuEdit[-which(leuEdit$conc.pNA==80 & leuEdit$pNA.subs == 1), ]
  leuEdit <- leuEdit[-which(leuEdit$conc.pNA==280 & leuEdit$pNA.subs == 6 & leuEdit$fl < 1000), ]
  leuEdit <- leuEdit[-which(leuEdit$conc.pNA==320 & leuEdit$pNA.subs == 2 & leuEdit$fl < 1000), ]
  leuEdit <- leuEdit[-which(leuEdit$conc.pNA==80 & leuEdit$pNA.subs == 9 & leuEdit$fl < 1500), ]  
  leuEdit <- leuEdit[-which(leuEdit$conc.pNA==120 & leuEdit$pNA.subs == 9 & leuEdit$fl < 1500), ]  
  leuEdit <- leuEdit[-which(leuEdit$conc.pNA==160 & leuEdit$pNA.subs == 9 & leuEdit$fl < 1500), ]  
  
  leuNCdataCheckPost <- ggplot(leuEdit, aes(x=incTime, y=fl)) +
    geom_point() + geom_line() +
    facet_grid(pNA.subs~conc.pNA) +
    scale_x_continuous("incubation time, hrs") +
    ggtitle("edited leu-AMC Bogue Sound")
  if(printDataCheck) {
    print(leuNCdataCheckPost)
  }
  
  # rename raw data, to ultimately put it into a list
  leuNCrawFl <- leuEdit
  rm(leuRaw)
  rm(leuEdit)
  
  #############
  # Pro-AMC in Bogue Sound
  #############

  # Read in raw data
  #proRaw <- read.csv("2013-03-25 pro competitive inhibition IMS.csv")
  #X2013_03_25_pro_competitive_inhibition_IMS <- proRaw
  #save(X2013_03_25_pro_competitive_inhibition_IMS, file="data/X2013_03_25_pro_competitive_inhibition_IMS.Rdata")
  load(file="data/X2013_03_25_pro_competitive_inhibition_IMS.Rdata")
  proRaw <- X2013_03_25_pro_competitive_inhibition_IMS
  
  # Change times to relative incubation time
  proRaw$Rtimes <- ymd_hms(paste("2013-03-25 ", as.character(proRaw$time), ":00", sep=""))
  proRaw$incTime <- as.numeric(proRaw$Rtimes - min(proRaw$Rtimes))/3600
  attr(proRaw$incTime, "units") <- "hours"
 
  # Plot the raw data
  proNCdataCheckPre <- ggplot(proRaw, aes(x=incTime, y=fl)) +
    geom_point() + geom_line() +
    facet_grid(pNA.subs~conc.pNA) +
    scale_x_continuous("incubation time, hrs") + 
    ggtitle("unedited pro-AMC Bogue Sound")
  if(printDataCheck) {
    print(proNCdataCheckPre)
  }
  
  # Excise bad data - only one bad data point
  proEdit <- proRaw
  proEdit <- proRaw[-which(proRaw$conc.pNA==360 & proRaw$pNA.subs=="killed" & proRaw$fl==1795), ]
  
  proNCdataCheckPost <- ggplot(proEdit, aes(x=incTime, y=fl)) +
    geom_point() + geom_line() +
    facet_grid(pNA.subs~conc.pNA) +
    scale_x_continuous("incubation time, hrs") +
    ggtitle("edited pro-AMC Bogue Sound")
  if(printDataCheck) {
    print(proNCdataCheckPost)
  }
  
  proNCrawFl <- proEdit
  rm(proRaw)
  rm(proEdit)
  
  # Create a list to hold all the raw data (could use a data frame I suppose)
  readInhibData <- list(argTNrawFl=argTNrawFl,
                        leuTNrawFl=leuTNrawFl,
                        argNCrawFl=argNCrawFl,
                        leuNCrawFl=leuNCrawFl,
                        proNCrawFl=proNCrawFl)
  
  #readInhibData <- llply(readInhibData, corr_pNA_conc)
  readInhibData
  
}