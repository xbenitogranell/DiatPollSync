#code for analyze Lake Llaviucu pollen data with PrC and GAM

# load libraries for functions used
library(tidyverse)
library(mgcv)
library(vegan)
library(analogue)
library(ggplot2)
library(cowplot)
library(gratia)


## Model temporal contributions of pollen PrC and Ti on diatom PrC
interpolatedData <- read.csv("outputs/principalcurves_ti_interp.csv")


# model temporal contributions of the covariates to diatom PrC
library(mgcv)
# this chunk run the model weighting in the elapsed time of the pollen (Majoi) core
mod1 <- gam(diatPrC ~ s(Age, k=10) + s(Ti) + s(pollenPrC) + s(agropastPrC, k=15),
            data = interpolatedData, method = "REML", 
            weights = elapsedTime / mean(elapsedTime),
            select = TRUE, family = gaussian(link="identity"),
            na.action = na.omit,
            knots = list(negAge=quantile(interpolatedData$Age, seq(0,1, length=10)))) #places notes at the deciles of sample ages

# this chunk run the model weighting in the elapsed time of the diatom core
mod1 <- gam(diatPrC ~ s(Age, k=10) + s(Ti) + s(pollenPrC) + s(agropastPrC, k=15),
            data = interpolatedData, method = "REML", 
            weights = elapsedTime_diat / mean(elapsedTime_diat),
            select = TRUE, family = gaussian(link="identity"),
            na.action = na.omit) 

pacf(residuals(mod1)) # indicates AR1

plot(mod1, page=1, scale = 0)
gam.check(mod1)
summary(mod1)

#accounting for temporal autocorrelation
mod1.car <- gamm(diatPrC ~ s(Age, k=10) + s(Ti) + s(pollenPrC) + s(agropastPrC, k=15),
                 data = interpolatedData, method = "REML", 
                 weights = elapsedTime/ mean(elapsedTime),
                 family = gaussian(link="identity"),
                 correlation = corCAR1(form = ~ Age),
                 knots = list(negAge=quantile(interpolatedData$Age, seq(0,1, length=10)))) #places notes at the deciles of sample ages


pacf(residuals(mod1.car$lme)) # indicates AR1

plot(mod1.car$gam, page=1, scale = 0)
summary(mod1.car$gam)
gam.check(mod1.car$gam)


#Compare different model fits using AIC
AIC_table <- AIC(mod1,mod1.car$lme)%>%
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
var <- predGam$fit.s.agropastPrC.
se.var <- predGam$se.fit.s.agropastPrC.

predGamPlt <- ggplot(predGam, aes(x = Age, y = var)) +
  geom_line() +
  geom_ribbon(aes(ymin = var + (2 * se.var), ymax = var - (2 * se.var), alpha = 0.1)) +
  #scale_x_reverse() +
  labs(y = "Diatom GAM PrC", x = "Cal yr BP", title = "")+
  theme(legend.position = "none")+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()+ theme(legend.position = "none")+
  ggtitle("")

predGamPlt

ggsave("outputs/diatomPrC_GAM_covariates.png",
       plot = predGamPlt ,
       width = 10,
       height=8,
       units="in",
       dpi = 400)


### fit GAM Gausian Location-Scale model
# read PrC pollen and diatom datasets
pollenPrC <- read.csv("outputs/pollen-PrC.csv", row.names = 1)
agropastolarismPrC <- read.csv( "outputs/agropastoralism-PrC.csv", row.names = 1)
diatomPrC <- readRDS("outputs/PrC-cores-diatoms.rds") %>%
  filter(Lake=="llaviucu")


#Transform dataframe to include elapsedtime ("years mud slice") and standarized PrC scores
diatoms <- transform(diatomPrC, Age = upper_age, negAge = - upper_age, elapsedTime = abs(upper_age - lower_age), 
                     rootPrC = sqrt(PrC), logPrC = log10(PrC+0.25))
diatoms <- arrange(diatoms,desc(negAge))

modDiatom <- gam(list(rootPrC ~ s(negAge, k = 30, bs = "ad"),
                      ~ elapsedTime + s(negAge, k=15)),
                 data = diatoms, method = "REML", family = gaulss(link = list("identity", "logb")))

summary(modDiatom)
gam.check(modDiatom)
plot(modDiatom, pages = 1, scale = 0, shade = TRUE)

#Transform dataframe to include elapsedtime ("years mud slice") and standarized PrC scores
pollen <- transform(pollenPrC, Age = upper_age, negAge = - upper_age, elapsedTime = abs(upper_age - lower_age), 
                    rootPrC = sqrt(PrC), logPrC = log10(PrC+0.25))
pollen <- arrange(pollen,desc(negAge))

modPollen <- gam(list(PrC ~ s(negAge, k = 30, bs = "ad"),
                      ~ elapsedTime + s(negAge, k=15)),
                 data = pollen, method = "REML", family = gaulss(link = list("identity", "logb")))

summary(modPollen)
gam.check(modPollen)
plot(modPollen, pages = 1, scale = 0, shade = TRUE)

#Transform dataframe to include elapsedtime ("years mud slice") and standarized PrC scores
agropastolarism <- transform(agropastolarismPrC, Age = upper_age, negAge = - upper_age, elapsedTime = abs(upper_age - lower_age), 
                    rootPrC = sqrt(PrC), logPrC = log10(PrC+0.25))
agropastolarism <- arrange(agropastolarism,desc(negAge))

modAgrosp <- gam(list(rootPrC ~ s(negAge, k = 30, bs = "ad"),
                      ~ elapsedTime + s(negAge, k=15)),
                 data = agropastolarism, method = "REML", family = gaulss(link = list("identity", "logb")))

summary(modAgrosp)
gam.check(modAgrosp)
plot(modAgrosp, pages = 1, scale = 0, shade = TRUE)


## Check temporal resolution of the time series
dif <- c(NA, diff(diatoms$upper_age)) #NA
median(dif[-1])
plot.ts(dif)

dif <- c(NA, diff(pollen$upper_age)) #NA
median(dif[-1])
plot.ts(dif)

dif <- c(NA, diff(agropastolarism$upper_age)) #NA
median(dif[-1])
plot.ts(dif)


## Predictions over range of both data sets
Nnew <- 500
elapsed <- 20
createNewData <- function(Age, n, elapsed) {
  out <- data.frame(Age = seq(min(Age), max(Age), length.out = n))
  out <- transform(out, negAge = - Age, elapsedTime = elapsed)
  out
}
newDiat <- with(diatoms, createNewData(Age, Nnew, elapsed))
newPoll <- with(pollen, createNewData(Age, Nnew, elapsed))
newAgropast <- with(agropastolarism, createNewData(Age, Nnew, elapsed))

predGLSS <- function(model, newdata) {
  N <- nrow(newdata)
  p <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  fam <- family(model)
  mulink <- fam[["linfo"]][[1L]][["linkinv"]]
  sigmalink <- fam[["linfo"]][[2L]][["linkinv"]]
  res <- data.frame(Age = rep(newdata$Age, 2),
                    negAge = rep(newdata$negAge, 2),
                    term = rep(c("mu","sigma"), each = N),
                    fitted = c(mulink(p$fit[,1]), 1 / sigmalink(p$fit[,2])),
                    upper  = c(mulink(p$fit[,1] + (2 * p$se.fit[,1])),
                               1 / sigmalink(p$fit[,2] + (2 * p$se.fit[,2]))),
                    lower  = c(mulink(p$fit[,1] - (2 * p$se.fit[,1])),
                               1 / sigmalink(p$fit[,2] - (2 * p$se.fit[,2]))))
  res
}

predDiat <- predGLSS(modDiatom, newDiat)
predPoll <- predGLSS(modPollen, newPoll)
predAgropast <- predGLSS(modAgrosp, newPoll)

## plot
predDiatPlt <- ggplot(predDiat, aes(x = Age, y = fitted, group = term)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.1) +
  geom_line() +
  facet_wrap(~ term, nrow = 2, labeller = label_parsed, scales = "free_y") +
  scale_x_reverse() +
  labs(y = "Fitted", x = "Age", title = "Diatoms") +
  theme_classic()
predDiatPlt

predPollPlt <- ggplot(predPoll, aes(x = Age, y = fitted, group = term)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.1) +
  geom_line() +
  facet_wrap(~ term, nrow = 2, labeller = label_parsed, scales = "free_y") +
  scale_x_reverse() +
  labs(y = "Fitted", x = "Age", title = "Pollen") +
  theme_classic() 
predPollPlt

predAgropastPlt <- ggplot(predAgropast, aes(x = Age, y = fitted, group = term)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "black", alpha = 0.1) +
  geom_line() +
  facet_wrap(~ term, nrow = 2, labeller = label_parsed, scales = "free_y") +
  scale_x_reverse() +
  labs(y = "Fitted", x = "Age", title = "Pollen (agrospastolarism indicators)") +
  theme_classic()
predAgropastPlt

# combine the three plots
All <- plot_grid(predDiatPlt, predPollPlt, predAgropastPlt,ncol = 3, align = "hv")
All

ggsave("outputs/llaviucu-diatom-pollen-GAULSS.png",
       plot = All,
       width=8,
       height=6,
       units="in",
       dpi = 400)




