## Model temporal contributions of pollen and agrospatoralism PrC and Ti on diatom PrC

#Load libraries for functions used
library(tidyverse)
library(mgcv)

# Read in the interpolated dataset of Ti, pollen and agropastoralism data on diatom record
interpolatedData <- read.csv("outputs/principalcurves_ti_interp_finer.csv")
interpolatedData <- read.csv("outputs/principalcurves_ti_interp_coarser.csv")

# this chunk run a multivariate GAM weighting in the elapsed time of the pollen (Majoi) core
mod1 <- gam(diatPrC ~ s(Age, k=20, bs="ad") + s(Ti) + s(pollenPrC, k=10, bs="ad") + s(agropastPrC,k=10, bs="ad"),
            data = interpolatedData, method = "REML", 
            weights = elapsedTime_ti / mean(elapsedTime_ti),
            select = TRUE, family = gaussian(link="identity"),
            na.action = na.omit,
            knots = list(negAge=quantile(interpolatedData$Age, seq(0,1, length=10)))) #places notes at the deciles of sample ages

pacf(residuals(mod1)) # indicates AR1
plot(mod1, page=1, scale = 0)
gam.check(mod1)
summary(mod1)


# this chunk run the model weighting in the elapsed time of the diatom core (ULL! not applicable for the coarser dataset)
mod2 <- gam(diatPrC ~ s(Age, k=20) + s(Ti, k=15, bs="ad") + s(pollenPrC) + s(agropastPrC, k=10, bs="ad"),
            data = interpolatedData, method = "REML", 
            weights = elapsedTime_diat / mean(elapsedTime_diat),
            select = TRUE, family = gaussian(link="identity"),
            na.action = na.omit,
            knots = list(negAge=quantile(interpolatedData$Age, seq(0,1, length=10)))) 

pacf(residuals(mod2)) # indicates non AR1
plot(mod2, page=1, scale=0)
gam.check(mod2)
summary(mod2)



#accounting for temporal autocorrelation
mod1.car <- gamm(diatPrC ~ s(Age, k=20) + s(Ti) + s(pollenPrC) + s(agropastPrC),
                 data = interpolatedData, method = "REML", 
                 weights = elapsedTime_ti/ mean(elapsedTime_ti),
                 family = gaussian(link="identity"),
                 correlation = corCAR1(form = ~ Age),
                 knots = list(negAge=quantile(interpolatedData$Age, seq(0,1, length=10)))) #places notes at the deciles of sample ages


pacf(residuals(mod1.car$lme)) # indicates AR1

plot(mod1.car$gam, page=1, scale = 0)
summary(mod1.car$gam)
gam.check(mod1.car$gam)


#Compare different model fits using AIC
AIC_table <- AIC(mod1, mod1.car$lme)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

AIC_table

# Predict GAM
## data to predict at
#Check NAs and remove rows
row.has.na <- apply(interpolatedData, 1, function(x){any(is.na(x))})
sum(row.has.na)
interpolatedData <- interpolatedData[!row.has.na,]

## Predict
predGam <- cbind(interpolatedData, 
                 data.frame(predict.gam(mod1, interpolatedData, 
                                        type = "terms" , se.fit = TRUE)))

#plot
var <- predGam$fit.s.pollenPrC.
se.var <- predGam$se.fit.s.pollenPrC.

predGamPlt <- ggplot(predGam, aes(x = Age, y = var)) +
  geom_line() +
  geom_ribbon(aes(ymin = var + (2 * se.var), ymax = var - (2 * se.var), alpha = 0.1)) +
  #scale_x_reverse() +
  labs(y = "Diatom GAM PrC", x = "Cal yr BP", title = "")+
  theme(legend.position = "none")+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()+ theme(legend.position = "none")+
  ggtitle(expression("Agropastoralism" %->% "Diatoms"))

predGamPlt

ggsave("outputs/diatomPrC_GAM_covariates.png",
       plot = predGamPlt ,
       width = 10,
       height=8,
       units="in",
       dpi = 400)

####



# Read pollen ratios (Majoi's data)
llaviucu_ratios <- read.csv("data/llaviucu_pollen.csv")  %>%
  mutate(Age = age*(-1)+1950) %>%
  filter(Age > 0) %>%
  mutate(Sporormiella=log10(Sporormiella+1))%>%
  select(age, Age, Sporormiella, Zeamays, Grassland, Woodland, Charcoal.cc.)


# Read in pollen aand agropastoralism PrC
pollen_PrC <- read.csv("outputs/pollen-PrC.csv") 
agropast_PrC <- read.csv("outputs/agropastoralism-PrC.csv")

# Read in Diatom PrC
diatomPrC <- readRDS("outputs/PrC-cores-diatoms.rds") %>%
  filter(Lake=="llaviucu")

#subsetting Llaviucu data and removing NAs
char <- na.omit(llaviucu_ratios[, c("age","Charcoal.cc.")])

# check temporal resolution of diatom 
diatomPrC <- diatomPrC[order(diatomPrC$Age),]
diff(diatomPrC$Age)
dif <- c(NA, diff(diatomPrC$Age)) #NA
median(dif[-1])
plot.ts(dif)

# check temporal resolution of charcoal data 
char <- char[char$age >= 0,]
diff(char$age)
plot.ts(char$Charcoal.cc.)


# Read XRF Llaviucu Tobi's data and merge with diatom PrC
llaviucu_diatoms_xrf <- read.csv("data/llaviucu_xrf.csv") %>%
  mutate(Depth = depth) %>%
  left_join(diatomPrC, by="Depth") %>% #this line is the diatoms PrC 
  mutate(elapsedTime = abs(upper_age - lower_age)) %>%
  mutate(rootPrC=sqrt(PrC)) %>%
  mutate(Fe_Mn = Fe/Mn) %>%
  mutate(Si_Ti = Si/Ti) %>%
  mutate(K_Ti = K/Ti) %>%
  mutate(Age=upper_age*(-1)+1950) %>%
  filter(Age > 0)


## Linear interpolation
charcInt<-approx(llaviucu_ratios$age, llaviucu_ratios$Charcoal.cc.,  method="constant", xout=0:2000)
length(charcInt$y) # This is the interpolated data
plot.ts(charcInt$y) # This is the interpolated data

  #matching diatom ages with charcoal interpolated
  llaviucu_charc<- as.data.frame(charcInt$y[llaviucu_diatoms_xrf$Age]) # This is the result of when diatom visits are present
  plot.ts(llaviucu_charc)

## Linear interpolation
sporInt<-approx(llaviucu_ratios$age, llaviucu_ratios$Sporormiella,  method="constant", xout=0:2000)
plot(sporInt$y) # This is the interpolated data
length(sporInt$y) # This is the interpolated data

  # Extract Sporiormella data when each Diatom visit took place
  llaviucu_spor <- as.data.frame(sporInt$y[llaviucu_diatoms_xrf$Age]) # This is the interpolated data
  plot.ts(llaviucu_spor)
  
# Create dataframe with reponse (diatom PrC) and explanatory variables
sitesData <- as.data.frame(cbind(llaviucu_diatoms_xrf$Age, llaviucu_diatoms_xrf$Depth, 
                                 llaviucu_diatoms_xrf$negAge,
                                 llaviucu_diatoms_xrf$elapsedTime, llaviucu_diatoms_xrf$PrC, 
                                 llaviucu_diatoms_xrf$rootPrC,
                                 llaviucu_spor[,1], 
                                 llaviucu_charc[,1],
                                 llaviucu_diatoms_xrf$Fe_Mn, llaviucu_diatoms_xrf$Si_Ti,
                                 llaviucu_diatoms_xrf$K_Ti))
colnames(sitesData) <- c("Age", "depth", "negAge", "elapsedTime", "PrC", "rootPrC", "spor", "charcoal", "Fe_Mn", "Si_Ti", "K_Ti")


## Fit GAM using gam(), which assumes residuals are independent
## Llaviucu
mod1 <- gam(PrC ~ s(negAge, k=10) + s(spor) + s(charcoal) + s(K_Ti) + s(Fe_Mn),
            weights = elapsedTime / mean(elapsedTime),
            data = sitesData, method = "REML", select = TRUE, family = gaussian(link="identity"),
            knots = list(negAge=quantile(sitesData$Age, seq(0,1, length=10)))) #places notes at the deciles of sample ages

plot(mod1, page=1, scale = 0)
summary(mod1)
gam.check(mod1)


## add predictions from GAM model
data_pred <- sitesData %>%
  select(Age, negAge, spor, charcoal, Si_Ti, Fe_Mn, K_Ti)

## Predict
predGam <- cbind(sitesData, data.frame(predict.gam(mod1, data_pred, 
                                                   type = "terms" , se.fit = TRUE)))

#plot
var <- predGam$fit.s.K_Ti.
se.var <- predGam$se.fit.s.K_Ti.

predGamPlt <- ggplot(predGam, aes(x = Age, y = var)) +
  geom_line() +
  geom_ribbon(aes(ymin = var + (2 * se.var), ymax = var - (2 * se.var), alpha = 0.1)) +
  #scale_x_reverse() +
  labs(y = "Diatom GAM PrC", x = "Cal yr BP", title = "")+
  theme(legend.position = "none")+
  geom_hline(yintercept=0, linetype="dashed")+
  ggtitle("")
predGamPlt

