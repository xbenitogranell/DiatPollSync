
library(tidyverse)

llaviucu_xrf <- read.csv("data/XRFTransformed3000yrs.csv")
llaviucu_xrf <- rename(llaviucu_xrf, c("age"="age_calBP")) 
Ti_xrf <- na.omit(llaviucu_xrf[, c("age", "Ti")])

## read in PrC datasets: diatoms, pollen and agropastoralism
diatomPrC <- readRDS("outputs/PrC-cores-diatoms.rds") %>%
  filter(Lake=="llaviucu") %>%
  mutate(Age=upper_age*(-1)+1950) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) %>%
  filter(Age > 0)
diatomPrC <- diatomPrC[order(diatomPrC$Age),] #order time
diatomPrC <- na.omit(diatomPrC[, c("Age", "PrC", "elapsedTime")])
diatomPrC$Age <- round(diatomPrC$Age, 1)
mean(diff(diatomPrC$Age))

pollenPrC <- read.csv("outputs/pollen-PrC.csv") %>%
  mutate(Age=upper_age*(-1)+1950) %>%
  filter(Age > 0)
pollenPrC <- pollenPrC[order(pollenPrC$Age),] #order time

agropastPrC <- read.csv("outputs/agropastoralism-PrC.csv") %>%
  mutate(Age=upper_age*(-1)+1950) %>%
  filter(Age > 0)
agropastPrC <- agropastPrC[order(agropastPrC$Age),] #order time


#linear interpolation
TiInt<-approx(Ti_xrf$age, Ti_xrf$Ti,  method="constant", xout=0:2000)
TiInt$y # This is the vector of Ti data (years 0 to 3000)
plot(TiInt$y) # This is the interpolated data
length(TiInt$y) # This is the interpolated data

plot.ts(diatom_Ti)
length(diatom_Ti[,1])

pollenInt<-approx(pollenPrC$Age, pollenPrC$PrC,  method="constant", xout=0:2000)
plot(pollenInt$y, type="l") # This is the interpolated data

agropastInt<-approx(agropastPrC$Age, agropastPrC$PrC,  method="constant", xout=0:2000)
plot(agropastInt$y, type="l") # This is the interpolated data

diatom_Ti <- as.data.frame(TiInt$y[diatomPrC$Age]) # This is the result of when diatom visits are present
diatom_pollen <- as.data.frame(pollenInt$y[diatomPrC$Age]) # This is the result of when diatom visits are present
diatom_agropast <- as.data.frame(agropastInt$y[diatomPrC$Age])

# Merge all datasets
interpolatedData <- cbind(diatomPrC,  diatom_Ti, diatom_pollen, diatom_agropast)
colnames(interpolatedData) <- c("Age", "diatomPrC", "elapsedTime" ,"Ti", "pollenPrC", "agropastPrC")

write.csv(interpolatedData, "outputs/interpolatedDatasets.csv")    
      
library(mgcv)
mod1 <- gam(diatomPrC ~ s(Age, k=10) + s(Ti) + s(pollenPrC) + s(agropastPrC),
            weights = elapsedTime / mean(elapsedTime),
            data = interpolatedData, method = "REML", select = TRUE, family = gaussian(link="identity"),
            na.action = na.omit) #places notes at the deciles of sample ages

plot(mod1, page=1, scale = 0)
summary(mod1)
gam.check(mod1)

                        
# ## LOESS interpolation
# source("scripts/functions_custom.R")
# 
# char.interpolated<-interpolateDatasets(
#   datasets.list=list(diatomPrC=diatomPrC), 
#   age.column.name="Age", 
#   interpolation.time.step=30 #ka
# )



