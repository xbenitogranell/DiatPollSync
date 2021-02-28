# interpolations and response-driver relationships using GAMs

library(mgcv)
library(tidyverse)

# read in XRF Llaviucu data
llaviucu_xrf <- read.csv("data/XRFTransformed3000yrs.csv") %>% rename(age=age_calBP) 
Ti_xrf <- na.omit(llaviucu_xrf[, c("age", "Ti")])


#linear interpolation
#TiInt<-approx(Ti_xrf$age, Ti_xrf$Ti,  method="constant", xout=0:3000)
#TiInt$y # This is the vector of Ti data (years 0 to 3000)
#plot(TiInt$y) # This is the interpolated data
#length(TiInt$y) # This is the interpolated data


# Read diatom-based rate of change
diat <- read.csv("outputs/pollen-agropast-diatom-deriv_24timeSteps.csv") %>%
  filter(proxy=="diatoms") %>% select(negAge,deriv_mean_med) %>%
  mutate(age=-negAge) %>% 
  mutate(age=round(age,0)) %>%
  data.frame()

# Read agropastoralism indicators rate of change
agropast <- read.csv("outputs/pollen-agropast-diatom-deriv_24timeSteps.csv") %>%
  filter(proxy=="agropastoralism") %>% select(negAge,deriv_mean_med) %>%
  mutate(age=-negAge) %>% 
  mutate(age=round(age,0)) %>%
  data.frame()

# Read pollen rate of change
pollen <- read.csv("outputs/pollen-agropast-diatom-deriv_24timeSteps.csv") %>%
  filter(proxy=="pollen") %>% select(negAge,deriv_mean_med) %>%
  mutate(age=-negAge) %>% 
  mutate(age=round(age,0)) %>%
  data.frame()

# Interpolate Ti into diatom rate-of-change temporal resolution 
Ti_diat_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, round(diat$age, 0))$y)
Ti_diat_i$age <- diat$age
colnames(Ti_diat_i) <- c("Ti", "age")
Ti_diat_i$diat_deriv <- diat$deriv_mean_med 

plot.ts(Ti_diat_i[,1])
plot(Ti_diat_i$Ti, Ti_diat_i$diat_deriv)
cor(Ti_diat_i[,1], Ti_diat_i[,3],  use = "complete")

# Interpolate Ti into agropastoralism rate-of-change temporal resolution 
Ti_agropast_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, round(agropast$age, 0))$y)
Ti_agropast_i$age <- agropast$age
colnames(Ti_agropast_i) <- c("Ti", "age")
Ti_agropast_i$agropast_deriv <- agropast$deriv_mean_med 

plot(Ti_agropast_i$Ti, Ti_agropast_i$agropast_deriv)
cor(Ti_agropast_i[,1], Ti_agropast_i[,3],  use = "complete")

# Interpolate Ti into pollen rate-of-change temporal resolution 
Ti_pollen_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, round(pollen$age, 0))$y)
Ti_pollen_i$age <- pollen$age
colnames(Ti_pollen_i) <- c("Ti", "age")
Ti_pollen_i$pollen_deriv <- pollen$deriv_mean_med 

plot(Ti_pollen_i$Ti, Ti_pollen_i$pollen_deriv)
cor(Ti_pollen_i[,1], Ti_pollen_i[,3],  use = "complete")


#Create a common dataframe with Ti and biological rate-of-change
Ti_bio <- data.frame(Ti_diat_i$age ,Ti_diat_i$Ti, Ti_diat_i$diat_deriv, Ti_agropast_i$agropast_deriv, Ti_pollen_i$pollen_deriv)
colnames(Ti_bio) <- c("age", "Ti", "diat", "agropast", "pollen")

Ti_bio <- Ti_bio %>% gather(key=taxa, value, -age, -Ti) %>% mutate(taxa_f=as.factor(taxa))

#plot response-driver relationships
ggplot(Ti_bio, aes(Ti, value, colour = taxa_f)) +
  geom_point() +
  geom_smooth()



# ## LOESS interpolation (Blas Benito's code)
# source("scripts/functions_custom.R")
# # 
# char.interpolated<-interpolateDatasets(
#  datasets.list=list(Ti_xrf=Ti), 
#  age.column.name="age", 
#  interpolation.time.step=0.024 #ka
#  )


## read in PrC datasets: diatoms, pollen and agropastoralism
diatomPrC <- readRDS("outputs/PrC-cores-diatoms.rds") %>%
  filter(Lake=="llaviucu") %>%
  #mutate(Age=upper_age*(-1)+1950) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 

diatomPrC <- diatomPrC[order(diatomPrC$Age),] #order time
diatomPrC <- na.omit(diatomPrC[, c("Age", "PrC", "elapsedTime")])
diatomPrC$Age <- round(diatomPrC$Age, 1)
mean(diff(diatomPrC$Age))

# linear interpolation
Ti_diat_prc_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, diatomPrC$Age)$y)
Ti_diat_prc_i$age <- diatomPrC$Age
colnames(Ti_diat_prc_i) <- c("Ti", "age")
Ti_diat_prc_i$diat_prc <- diatomPrC$PrC
Ti_diat_prc_i$elapsedTime <- diatomPrC$elapsedTime

plot(Ti_diat_prc_i[,1], Ti_diat_prc_i[,3])

## Pollen PrC
pollenPrC <- read.csv("outputs/pollen-PrC.csv") %>%
  mutate(Age=upper_age) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 
  #filter(Age > 0)
pollenPrC <- pollenPrC[order(pollenPrC$Age),] #order time

# linear interpolation
Ti_pollen_prc_i <- as.data.frame(approx(Ti_xrf$age, Ti_xrf$Ti, pollenPrC$Age)$y)
Ti_pollen_prc_i$age <- pollenPrC$Age
colnames(Ti_pollen_prc_i) <- c("Ti", "age")
Ti_pollen_prc_i$pollen_prc <- pollenPrC$PrC
Ti_pollen_prc_i$elapsedTime <- pollenPrC$elapsedTime

plot(Ti_pollen_prc_i[,1], Ti_pollen_prc_i[,3])

## Agropastoralism PrC
agropastPrC <- read.csv("outputs/agropastoralism-PrC.csv") %>%
  mutate(Age=upper_age) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 

agropastPrC <- agropastPrC[order(agropastPrC$Age),] #order time

# linear interpolation
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
write.csv(interpolatedData, "outputs/principalcurves_ti_interp.csv")



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
interpolatedDatasets <- read.csv("outputs/interpolatedDatasets.csv")
mod1 <- gam(diatomPrC ~ s(Age, k=10) + s(Ti) + s(pollenPrC) + s(agropastPrC),
            weights = elapsedTime / mean(elapsedTime),
            data = interpolatedDatasets, method = "REML", select = TRUE, family = gaussian(link="identity"),
            na.action = na.omit) #places notes at the deciles of sample ages

plot(mod1, page=1, scale = 0)
summary(mod1)
gam.check(mod1)

                        



