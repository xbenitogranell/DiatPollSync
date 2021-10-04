
##################
#### MUTI ########
##################

#install.packages("remotes")
#remotes::install_github("mdscheuerell/muti")
library(tidyverse)
library(muti)

derivAll <- read.csv("outputs/pollen-agropast-diatom-deriv_24timeSteps.csv", row.names = 1)
Ti <- read.csv("outputs/Ti-gam-derivatives_24timeSteps.csv", row.names=1) %>%
  select(5,6) 

## 

mean <- "deriv_mean_med"
sd <- "deriv_sd_med"

#extract datasets to compute muti among mean derivative values across all lake pairs
diat <- derivAll %>% filter(proxy=="diatoms") %>% select(deriv_mean_med,deriv_sd_med) 
pollen <- derivAll %>% filter(proxy=="pollen") %>% select(deriv_mean_med,deriv_sd_med) 
agropast <- derivAll %>% filter(proxy=="agropastoralism") %>% select(deriv_mean_med,deriv_sd_med) 
#Ti <- derivAll %>% filter(proxy=="Ti") %>% select(deriv_mean_med,deriv_sd_med) 

df <- cbind(diat,pollen,agropast,Ti)
colnames(df) <- c("diat_deriv_mean", "diat_sd_mean", "pollen_deriv_mean", "pollen_sd_mean", 
                  "agropast_deriv_mean", "agropast_sd_mean",
                  "ti_deriv_mean", "ti_sd_mean")

df_mean <- df[,c(1,3,5,7)]
df_sd <- df[,c(2,4,6,8)]


## Run the Mutual information analysis for mean derivatives
syncDeriv_mean <- list()
counter=0

for (i in 2:ncol(df_mean)) {
  counter = counter + 1
  x = df_mean[,1] #this is the diatom time series (response)
  y = df_mean[[i]] # these are pollen, agropastoralism and Ti
  syncDeriv_mean[[counter]] <- muti(x, y, sym = TRUE, normal=TRUE, mc=100, alpha=.05)
}

names(syncDeriv_mean) <- colnames(df_mean[,2:ncol(df_mean)])

## Run the Mutual information analysis for sd derivatives
syncSd_mean <- list()
counter=0

for (i in 2:ncol(df_sd)) {
  counter = counter + 1
  x = df_sd[,1] #this is the diatom time series (response)
  y = df_sd[[i]] # these are pollen and agropastoralism
  syncSd_mean[[counter]] <-muti(x, y, sym = TRUE, normal=TRUE, mc=100, alpha=.05)
}

names(syncSd_mean) <- colnames(df_sd[,2:ncol(df_sd)])

#extract deriv_mean lists
muti_deriv_mean <- plyr::ldply(syncDeriv_mean, data.frame)
colnames(muti_deriv_mean) <- c("var", "lag",  "MI", "MI_null")
muti_deriv_mean <- mutate(muti_deriv_mean,sig=ifelse(MI>MI_null, 1, 0)) #is MI significant?
muti_deriv_mean <- filter(muti_deriv_mean, lag<=0 & sig==1) #filter significant synchronous (lag=0) and delayed lags (negative to ask if time-delayed affects of pollen and agropastoralism affect diatoms?)

#extract sd_mean lists
muti_sd_mean <- plyr::ldply(syncSd_mean, data.frame)
colnames(muti_sd_mean) <- c("var", "lag",  "MI", "MI_null")
muti_sd_mean <- mutate(muti_sd_mean,sig=ifelse(MI>MI_null, 1, 0)) #is MI significant?
muti_sd_mean <- filter(muti_sd_mean, lag<=0 & sig==1) #filter significant synchronous (lag=0) and delayed lags (negative to ask if time-delayed affects of pollen and agropastoralism affect diatoms?)




pollen_agropast_diatom_deriv_24timeSteps <- read.csv("outputs/pollen-agropast-diatom-deriv_24timeSteps.csv", row.names=1)
head(pollen_agropast_diatom_deriv_24timeSteps)

Ti_gam_derivatives_24timeSteps <- read.csv("outputs/Ti-gam-derivatives_24timeSteps.csv", row.names=1) %>%
  select(5) %>%
  rename(titanium=deriv_mean_med)


pollen <- pollen_agropast_diatom_deriv_24timeSteps[,c(6,8)] %>%
  filter(proxy=="pollen") %>%
  mutate(deriv_pollen=deriv_mean_med) %>%
  select(-deriv_mean_med, -proxy) 

diatoms <- pollen_agropast_diatom_deriv_24timeSteps[,c(6,8)] %>%
  filter(proxy=="diatoms") %>%
  mutate(deriv_diat=deriv_mean_med) %>%
  select(-deriv_mean_med, -proxy)

agropast <- pollen_agropast_diatom_deriv_24timeSteps[,c(6,8)] %>%
  filter(proxy=="agropastoralism") %>%
  mutate(deriv_agropast=deriv_mean_med) %>%
  select(-deriv_mean_med, -proxy) %>%
  cbind(pollen, diatoms, Ti_gam_derivatives_24timeSteps)

df <- agropast
colnames(df) <- c("agropast", "pollen", "diatoms", "titanium")
df <- df[, c(3, 2, 1, 4)]

## muti analysis
library(muti)
syncCoupl <- list()
counter=0

for (i in 2:ncol(df)) {
  counter = counter + 1
  x = df$diatoms #this is the response 
  y = df[[i]]
  syncCoupl[[counter]] <- muti(x, y, sym = TRUE, normal=TRUE, mc=100, alpha=.05)
}

names(syncCoupl) <- colnames(df[,2:ncol(df)])

sync <- plyr::ldply(syncCoupl, data.frame) %>%
  mutate(sig=ifelse(MI_xy>MI_tv,1,0)) %>%
  filter(lag<=0)


