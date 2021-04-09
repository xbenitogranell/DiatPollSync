#code for analyze Lake Llaviucu pollen data with PrC and GAM

# load libraries for functions used
library(tidyverse)
library(mgcv)
library(vegan)
library(analogue)
library(ggplot2)
library(cowplot)
library(gratia)



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




