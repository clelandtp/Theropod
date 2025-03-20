
setwd("output")

.libPaths(paste("C:/Users/",Sys.getenv("USERNAME"),"/AppData/Local/R/win-library/4.3",sep=""))

library(tidyverse)
library(mzR)

MHfile <- list.files(pattern = "-MH.csv$")
print(MHfile)

MHfile2 <- gsub("\\.csv$","", list.files(pattern="\\-MH.csv$"))

MH <- read_csv(MHfile)
MH <- MH %>% 
  select(mz, intensity)
MH3 <- MH 
MH3['mz'] <- MH3$mz+MH3$mz*(20/1e6)
MH2 <- data.matrix(MH3)
MHlist <- list(MH2)

seqNum <- 1                   
acquisitionNum <- 1           
msLevel <- 1                  
polarity  <- 1                
peaksCount  <- 8390              
totIonCurrent <- sum(MH3$intensity)            
retentionTime <- 0.616            
basePeakMZ <-  MH3$mz[which.max(MH3$intensity)]              
basePeakIntensity <- max(MH3$intensity)        
collisionEnergy <- as.numeric(NA)          
ionisationEnergy  <- 0        
lowMZ <- min(MH3$mz)                    
highMZ  <- max(MH3$mz)                  
precursorScanNum <- as.numeric(NA)          
precursorMZ  <- as.numeric(NA)             
precursorCharge   <- as.numeric(NA)        
precursorIntensity   <- as.numeric(NA)     
mergedScan  <- as.numeric(NA)              
mergedResultScanNum <- as.numeric(NA)       
mergedResultStartScanNum <- as.numeric(NA)  
mergedResultEndScanNum <- as.numeric(NA)   
injectionTime <- 0.6            
filterString  <- "FTMS + p ESI sid=15.00 Full ms [500.00-2000.00]"            
spectrumId  <- "controllerType=0 controllerNumber=1 scan=1"              
centroided  <- FALSE              
ionMobilityDriftTime <- as.numeric(NA)     
isolationWindowTargetMZ  <- as.numeric(NA) 
isolationWindowLowerOffset <- as.numeric(NA)
isolationWindowUpperOffset <- as.numeric(NA)
scanWindowLowerLimit <- min(MH3$mz)     
scanWindowUpperLimit  <- max(MH3$mz)           

MHhdr <- data.frame(seqNum, acquisitionNum, msLevel, polarity, peaksCount, 
                    totIonCurrent, retentionTime, basePeakMZ, basePeakIntensity, 
                    collisionEnergy, ionisationEnergy, lowMZ, highMZ, precursorScanNum, 
                    precursorMZ, precursorCharge, precursorIntensity, mergedScan, mergedResultScanNum, 
                    mergedResultStartScanNum, mergedResultEndScanNum, injectionTime, filterString, 
                    spectrumId, centroided, ionMobilityDriftTime, isolationWindowTargetMZ, 
                    isolationWindowLowerOffset, isolationWindowUpperOffset, scanWindowLowerLimit, scanWindowUpperLimit
)


writeMSData(MHlist, paste0(MHfile2,".mzML"), MHhdr, backend = 'pwiz', outformat = 'mzml', rtime_seconds = TRUE)

Mfile <- list.files(pattern = "-M.csv$")
print(Mfile)

Mfile2 <- gsub("\\.csv$","", list.files(pattern="\\-M.csv$"))

M <- read_csv(Mfile)
M <- M %>% 
  select(mz, intensity)
M3 <- M 
M3['mz'] <- M3$mz+M3$mz*(20/1e6)
M2 <- data.matrix(M3)
Mlist <- list(M2)

seqNum <- 1                   
acquisitionNum <- 1           
msLevel <- 1                  
polarity  <- 1                
peaksCount  <- 8390              
totIonCurrent <- sum(M3$intensity)            
retentionTime <- 0.616            
basePeakMZ <-  MH3$mz[which.max(M3$intensity)]              
basePeakIntensity <- max(M3$intensity)        
collisionEnergy <- as.numeric(NA)          
ionisationEnergy  <- 0        
lowMZ <- min(MH3$mz)                    
highMZ  <- max(MH3$mz)                  
precursorScanNum <- as.numeric(NA)          
precursorMZ  <- as.numeric(NA)             
precursorCharge   <- as.numeric(NA)        
precursorIntensity   <- as.numeric(NA)     
mergedScan  <- as.numeric(NA)              
mergedResultScanNum <- as.numeric(NA)       
mergedResultStartScanNum <- as.numeric(NA)  
mergedResultEndScanNum <- as.numeric(NA)   
injectionTime <- 0.6            
filterString  <- "FTMS + p ESI Full ms [500.00-2000.00]"            
spectrumId  <- "controllerType=0 controllerNumber=1 scan=1"              
centroided  <- FALSE              
ionMobilityDriftTime <- as.numeric(NA)     
isolationWindowTargetMZ  <- as.numeric(NA) 
isolationWindowLowerOffset <- as.numeric(NA)
isolationWindowUpperOffset <- as.numeric(NA)
scanWindowLowerLimit <- min(M3$mz)     
scanWindowUpperLimit  <- max(M3$mz)              


Mhdr <- data.frame(seqNum, acquisitionNum, msLevel, polarity, peaksCount, 
                    totIonCurrent, retentionTime, basePeakMZ, basePeakIntensity, 
                    collisionEnergy, ionisationEnergy, lowMZ, highMZ, precursorScanNum, 
                    precursorMZ, precursorCharge, precursorIntensity, mergedScan, mergedResultScanNum, 
                    mergedResultStartScanNum, mergedResultEndScanNum, injectionTime, filterString, 
                    spectrumId, centroided, ionMobilityDriftTime, isolationWindowTargetMZ, 
                    isolationWindowLowerOffset, isolationWindowUpperOffset, scanWindowLowerLimit, scanWindowUpperLimit
)


writeMSData(Mlist, paste0(Mfile2,".mzML"), Mhdr, backend = 'pwiz', outformat = 'mzml', rtime_seconds = TRUE)



            