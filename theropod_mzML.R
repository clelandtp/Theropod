setwd("W:/Proteomics/Cleland/singleionprotmod092424/20250115_eno_90nM_0-3ulmin_0-5MIT_500ns_10000scans_240k/output")



library(tidyverse)
library(mzR)

MHfile <- list.files(pattern = "-MH.csv$")

MHfile2 <- gsub("\\.csv$","", list.files(pattern="\\-MH.csv$"))

MH <- read_csv(MHfile, id = "file_name")
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
    
lowMZ <- min(MH3$mz)                    
highMZ  <- max(MH3$mz)                  
           
spectrumId  <- "controllerType=0 controllerNumber=1 scan=1"              
centroided  <- FALSE              


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

Mfile2 <- gsub("\\.csv$","", list.files(pattern="\\-M.csv$"))

M <- read_csv(Mfile, id = "file_name")
M <- M %>% 
  select(mz, intensity)
M3 <- M 
#M3['mz'] <- M3$mz*(20.2/1e6)+M3$mz
M2 <- data.matrix(M3)
Mlist <- list(M2)

seqNum <- 1                   
acquisitionNum <- 1           
msLevel <- 1                  
polarity  <- 1                
      
lowMZ <- min(M3$mz)                    
highMZ  <- max(M3$mz)                  
           
spectrumId  <- "controllerType=0 controllerNumber=1 scan=1"              
centroided  <- FALSE              


Mhdr <- data.frame(seqNum, acquisitionNum, msLevel, polarity, peaksCount, 
                    totIonCurrent, retentionTime, basePeakMZ, basePeakIntensity, 
                    collisionEnergy, ionisationEnergy, lowMZ, highMZ, precursorScanNum, 
                    precursorMZ, precursorCharge, precursorIntensity, mergedScan, mergedResultScanNum, 
                    mergedResultStartScanNum, mergedResultEndScanNum, injectionTime, filterString, 
                    spectrumId, centroided, ionMobilityDriftTime, isolationWindowTargetMZ, 
                    isolationWindowLowerOffset, isolationWindowUpperOffset, scanWindowLowerLimit, scanWindowUpperLimit
)


writeMSData(Mlist, paste0(Mfile2,".mzML"), Mhdr, backend = 'pwiz', outformat = 'mzml', rtime_seconds = TRUE)



            