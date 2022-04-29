# interpolations among lake commmunity time series (diatoms, pollen and agropastoralism indicators)

library(mgcv)
library(tidyverse)

# read in XRF Llaviucu data (Majoi=pollen record)
llaviucu_xrf <- read.csv("data/XRFTransformed3000yrs.csv") %>% rename(age=age_calBP) 
Ti_xrf <- na.omit(llaviucu_xrf[, c("age", "Ti")])

## read in PrC datasets: diatoms, pollen and agropastoralism
diatomPrC <- readRDS("outputs/PrC-cores-diatoms.rds") %>%
  filter(Lake=="llaviucu") %>%
  #mutate(Age=upper_age*(-1)+1950) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 
length(diatomPrC$Age)

## Pollen PrC
pollenPrC <- read.csv("outputs/pollen-PrC.csv") %>%
  mutate(Age=upper_age) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 
#filter(Age > 0)
pollenPrC <- pollenPrC[order(pollenPrC$Age),] #order time
median(diff(pollenPrC$upper_age))
length(pollenPrC$upper_age)

## Agropastoralism PrC
agropastPrC <- read.csv("outputs/agropastoralism-PrC.csv") %>%
  mutate(Age=upper_age) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 

agropastPrC <- agropastPrC[order(agropastPrC$Age),] #order time
median(diff(agropastPrC$upper_age))

# Read XRF Llaviucu Tobi's data and merge with diatom PrC
llaviucu_diatoms_xrf <- read.csv("data/llaviucu_xrf.csv") %>%
  mutate(Depth = depth) %>%
  left_join(diatomPrC, by="Depth") %>% #this line is the diatoms PrC 
  mutate(elapsedTime = abs(upper_age - lower_age)) %>%
  mutate(rootPrC=sqrt(PrC)) %>%
  mutate(Fe_Mn = Fe/Mn) %>%
  mutate(Si_Ti = Si/Ti) %>%
  mutate(K_Ti = K/Ti) %>%
  filter(Age > 0)

#Arrange diatom PrC in time order
diatomPrC <- diatomPrC[order(diatomPrC$Age),] #order time
diatomPrC <- na.omit(diatomPrC[, c("Age", "PrC", "elapsedTime")])
diatomPrC$Age <- round(diatomPrC$Age, 1)
mean(diff(diatomPrC$Age))

### Interpolation to the finer resolution: diatoms core
# linear interpolation diatom ~ 2009 Ti (pollen record)
Ti_diat_prc_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, diatomPrC$Age)$y)
Ti_diat_prc_i$age <- diatomPrC$Age
colnames(Ti_diat_prc_i) <- c("Ti", "age")
Ti_diat_prc_i$diat_prc <- diatomPrC$PrC
diff(Ti_diat_prc_i$age)
Ti_diat_prc_i$elapsedTime <- diatomPrC$elapsedTime
plot(Ti_diat_prc_i[,1], Ti_diat_prc_i[,3])

# linear interpolation 2009 pollen ~ 2009 Ti
Ti_pollen_prc_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, pollenPrC$Age)$y)
Ti_pollen_prc_i$age <- pollenPrC$Age
colnames(Ti_pollen_prc_i) <- c("Ti", "age")
Ti_pollen_prc_i$pollen_prc <- pollenPrC$PrC
Ti_pollen_prc_i$elapsedTime <- pollenPrC$elapsedTime
plot(Ti_pollen_prc_i[,1], Ti_pollen_prc_i[,3])

# linear interpolation agropastoralism ~ Ti
Ti_agropast_prc_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, agropastPrC$Age)$y)
Ti_agropast_prc_i$age <- agropastPrC$Age
colnames(Ti_agropast_prc_i) <- c("Ti", "age")
Ti_agropast_prc_i$agrospast_prc <- agropastPrC$PrC

plot(Ti_agropast_prc_i[,1], Ti_agropast_prc_i[,3])

# Merge all datasets finer dataset (diatoms)
interpolatedData <- as.data.frame(cbind(Ti_diat_prc_i$age, Ti_diat_prc_i$Ti,
                                        Ti_pollen_prc_i$elapsedTime, #elapsed time in equal time resolution than Ti
                                        Ti_diat_prc_i$elapsedTime, #elapsed time in equal time resolution than diatoms
                                        Ti_diat_prc_i$diat_prc,
                                        Ti_pollen_prc_i$pollen_prc, Ti_agropast_prc_i$agrospast_prc))

colnames(interpolatedData) <- c("Age", "Ti", "elapsedTime_ti", "elapsedTime_diat", "diatPrC", "pollenPrC", "agropastPrC")
median(diff(interpolatedData$Age))

# save data
write.csv(interpolatedData, "outputs/principalcurves_ti_interp_finer.csv")
####

### Interpolation to the coarser dataset: pollen
# linear interpolation diatom to pollen PrCs
diatom_to_pollen_prc_i <- as.data.frame(approx(diatomPrC$Age, diatomPrC$PrC, pollenPrC$Age)$y)
diatom_to_pollen_prc_i$age <- pollenPrC$Age
colnames(diatom_to_pollen_prc_i) <- c("PrC", "age")
diatom_to_pollen_prc_i$pollen_prc <- pollenPrC$PrC
head(diatom_to_pollen_prc_i)
diff(diatom_to_pollen_prc_i$age)
diatom_to_pollen_prc_i$elapsedTime <- pollenPrC$elapsedTime
plot(diatom_to_pollen_prc_i[,1], diatom_to_pollen_prc_i[,3])

# linear interpolation diatom to agropastoralism PrCs
diatom_to_agropast_prc_i <- as.data.frame(approx(diatomPrC$Age, diatomPrC$PrC, agropastPrC$Age)$y)
diatom_to_agropast_prc_i$age <- agropastPrC$Age
colnames(diatom_to_agropast_prc_i) <- c("PrC", "age")
diatom_to_agropast_prc_i$agropast_prc <- agropastPrC$PrC
head(diatom_to_agropast_prc_i)
diff(diatom_to_agropast_prc_i$age)
diatom_to_agropast_prc_i$elapsedTime <- agropastPrC$elapsedTime
plot(diatom_to_agropast_prc_i[,1], diatom_to_agropast_prc_i[,3])


# linear interpolation 2009 pollen ~ 2009 Ti--this is the same for agropast
Ti_pollen_prc_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, pollenPrC$Age)$y)
Ti_pollen_prc_i$age <- pollenPrC$Age
colnames(Ti_pollen_prc_i) <- c("Ti", "age")
Ti_pollen_prc_i$pollen_prc <- pollenPrC$PrC
Ti_pollen_prc_i$elapsedTime <- pollenPrC$elapsedTime
plot(Ti_pollen_prc_i[,1], Ti_pollen_prc_i[,3])


# Merge all datasets finer dataset (diatoms)
interpolatedData <- as.data.frame(cbind(diatom_to_pollen_prc_i$age, diatom_to_pollen_prc_i$PrC,
                                        diatom_to_pollen_prc_i$elapsedTime, #elapsed time in equal time resolution than Ti
                                        diatom_to_pollen_prc_i$pollen_prc, 
                                        diatom_to_agropast_prc_i$agropast_prc,
                                        Ti_pollen_prc_i$Ti))

colnames(interpolatedData) <- c("Age", "diatPrC", "elapsedTime_ti", "pollenPrC", "agropastPrC", "Ti")
median(diff(interpolatedData$Age))

# save data
write.csv(interpolatedData, "outputs/principalcurves_ti_interp_coarser.csv")
###

### Old stuff

  # linear interpolation diatom to pollen PrCs
  diatom_to_pollen_prc_i <- as.data.frame(approx(diatomPrC$Age, diatomPrC$PrC, pollenPrC$Age)$y)
  diatom_to_pollen_prc_i$age <- pollenPrC$Age
  colnames(diatom_to_pollen_prc_i) <- c("PrC", "age")
  diatom_to_pollen_prc_i$pollen_prc <- pollenPrC$PrC
  head(diatom_to_pollen_prc_i)
  diff(diatom_to_pollen_prc_i$age)
  diatom_to_pollen_prc_i$elapsedTime <- pollenPrC$elapsedTime
  plot(diatom_to_pollen_prc_i[,1], diatom_to_pollen_prc_i[,3])

  # linear interpolation diatom to agropastoralism PrCs
  diatom_to_agropast_prc_i <- as.data.frame(approx(diatomPrC$Age, diatomPrC$PrC, agropastPrC$Age)$y)
  diatom_to_agropast_prc_i$age <- agropastPrC$Age
  colnames(diatom_to_agropast_prc_i) <- c("PrC", "age")
  diatom_to_agropast_prc_i$agropast_prc <- agropastPrC$PrC
  head(diatom_to_agropast_prc_i)
  diff(diatom_to_agropast_prc_i$age)
  diatom_to_agropast_prc_i$elapsedTime <- agropastPrC$elapsedTime
  plot(diatom_to_agropast_prc_i[,1], diatom_to_agropast_prc_i[,3])
  

  #linear interpolation diatom Ti ~ pollen
  TiInt<-approx(llaviucu_diatoms_xrf$Age, llaviucu_diatoms_xrf$Ti,  method="constant", xout=0:3000)
  pollenPrC<-pollenPrC[pollenPrC$Age>0,]
  llaviucu_TiInt<- as.data.frame(TiInt$y[pollenPrC$Age]) # This is the result of when Ti visits are present
  llaviucu_TiInt$age <- pollenPrC$Age
  colnames(llaviucu_TiInt) <- c("Ti", "age")
  
  #linear interpolation diatom PrC ~ pollen
  DiatPrcInt<-approx(llaviucu_diatoms_xrf$Age, llaviucu_diatoms_xrf$PrC,  method="constant", xout=0:3000)
  pollenPrC<-pollenPrC[pollenPrC$Age>0,]
  llaviucu_diatInt<- as.data.frame(DiatPrcInt$y[pollenPrC$Age]) # This is the result of when diatom PrC visits are present
  llaviucu_diatInt$age <- pollenPrC$Age
  colnames(llaviucu_diatInt) <- c("PrC", "age")

  #linear interpolation pollen PrC ~ pollen
  PollenPrcInt<-approx(pollenPrC$Age, pollenPrC$PrC,  method="constant", xout=0:3000)
  pollenPrC<-pollenPrC[pollenPrC$Age>0,]
  llaviucu_pollenInt<- as.data.frame(PollenPrcInt$y[pollenPrC$Age]) # This is the result of when diatom visits are present
  llaviucu_pollenInt$age <- pollenPrC$Age
  colnames(llaviucu_pollenInt) <- c("PrC", "age")
  


  #linear interpolation diatom PrC ~ pollen
  AgropastPrcInt<-approx(agropastPrC$Age, agropastPrC$PrC,  method="constant", xout=0:3000)
  pollenPrC<-pollenPrC[pollenPrC$Age>0,]
  llaviucu_agropastInt<- as.data.frame(AgropastPrcInt$y[pollenPrC$Age]) # This is the result of when diatom visits are present
  llaviucu_agropastInt$age <- pollenPrC$Age
  colnames(llaviucu_agropastInt) <- c("PrC", "age")


# Merge all datasets coarser dataset (pollen)
# interpolatedData_coarser <- as.data.frame(cbind(Ti_pollen_prc_i$age, #Ti_pollen_prc_i$Ti,
#                                         Ti_pollen_prc_i$elapsedTime, #elapsed time in equal time resolution than Ti
#                                         #Ti_diat_prc_i$elapsedTime, #elapsed time in equal time resolution than diatoms
#                                         Ti_diat_prc_i$diat_prc,
#                                         Ti_pollen_prc_i$pollen_prc, Ti_agropast_prc_i$agrospast_prc))

# Merge all datasets coarser dataset (pollen)
interpolatedData2 <- as.data.frame(cbind(llaviucu_TiInt$age, llaviucu_TiInt$Ti,
                                        pollenPrC$elapsedTime,
                                        llaviucu_diatInt$PrC,
                                        llaviucu_pollenInt$PrC, llaviucu_agropastInt$PrC))

colnames(interpolatedData2) <- c("Age", "Ti", "elapsedTime_ti", "diatPrC", "pollenPrC", "agropastPrC")
median(diff(interpolatedData2$Age))

# save data
write.csv(interpolatedData2, "outputs/principalcurves_ti_interp_coarser.csv")

### LOESS interpolation (Blas Benito's code)
source("scripts/functions_custom.R")
#

# pollen.i <- pollenPrC[, c("Age", "PrC")]
# mean(diff(pollen.i$Age))
# 
# char.interpolated<-interpolateDatasets(
#   datasets.list=list(pollen.i=pollen.i),
#   age.column.name="Age",
#   interpolation.time.step=60 #ka
#   )




