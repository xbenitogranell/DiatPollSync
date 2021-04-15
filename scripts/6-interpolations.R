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

# linear interpolation diatom ~ Ti (pollen record)
Ti_diat_prc_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, diatomPrC$Age)$y)
Ti_diat_prc_i$age <- diatomPrC$Age
colnames(Ti_diat_prc_i) <- c("Ti", "age")
Ti_diat_prc_i$diat_prc <- diatomPrC$PrC
diff(Ti_diat_prc_i$age)
Ti_diat_prc_i$elapsedTime <- diatomPrC$elapsedTime

plot(Ti_diat_prc_i[,1], Ti_diat_prc_i[,3])

## Pollen PrC
pollenPrC <- read.csv("outputs/pollen-PrC.csv") %>%
  mutate(Age=upper_age) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 
  #filter(Age > 0)
pollenPrC <- pollenPrC[order(pollenPrC$Age),] #order time

# linear interpolation pollen ~ Ti
Ti_pollen_prc_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, pollenPrC$Age)$y)
Ti_pollen_prc_i$age <- pollenPrC$Age
colnames(Ti_pollen_prc_i) <- c("Ti", "age")
Ti_pollen_prc_i$pollen_prc <- pollenPrC$PrC
Ti_pollen_prc_i$elapsedTime <- pollenPrC$elapsedTime

plot(Ti_pollen_prc_i[,1], Ti_pollen_prc_i[,3])

  #linear interpolation diatom Ti ~ pollen
  TiInt<-approx(llaviucu_diatoms_xrf$Age, llaviucu_diatoms_xrf$Ti,  method="constant", xout=0:3000)
  pollenPrC<-pollenPrC[pollenPrC$Age>0,]
  llaviucu_TiInt<- as.data.frame(TiInt$y[pollenPrC$Age]) # This is the result of when diatom visits are present
  llaviucu_TiInt$age <- pollenPrC$Age
  colnames(llaviucu_TiInt) <- c("Ti", "age")
  

## Agropastoralism PrC
agropastPrC <- read.csv("outputs/agropastoralism-PrC.csv") %>%
  mutate(Age=upper_age) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 

agropastPrC <- agropastPrC[order(agropastPrC$Age),] #order time

# linear interpolation agropastoralism ~ Ti
Ti_agropast_prc_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, agropastPrC$Age)$y)
Ti_agropast_prc_i$age <- agropastPrC$Age
colnames(Ti_agropast_prc_i) <- c("Ti", "age")
Ti_agropast_prc_i$agrospast_prc <- agropastPrC$PrC

plot(Ti_agropast_prc_i[,1], Ti_agropast_prc_i[,3])

# Merge all datasets
interpolatedData <- as.data.frame(cbind(Ti_diat_prc_i$age, Ti_diat_prc_i$Ti,
                          Ti_pollen_prc_i$elapsedTime, #elapsed time in equal time resolution than Ti
                          Ti_diat_prc_i$elapsedTime, #elapsed time in equal time resolution than diatoms
                          Ti_diat_prc_i$diat_prc,
                          Ti_pollen_prc_i$pollen_prc, Ti_agropast_prc_i$agrospast_prc))

colnames(interpolatedData) <- c("Age", "Ti", "elapsedTime_ti", "elapsedTime_diat", "diatPrC", "pollenPrC", "agropastPrC")

# save data
write.csv(interpolatedData, "outputs/principalcurves_ti_interp_v2.csv")


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




