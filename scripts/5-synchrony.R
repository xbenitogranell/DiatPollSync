# Synchrony models 

#############
## GLS models
#############
derivAll <- read.csv("outputs/pollen-agropast-diatom-deriv_24timeSteps.csv")

library(nlme)
library(tidyverse)

#extract datasets to get mean derivative values for each group
diat <- derivAll %>% filter(proxy=="diatoms") %>% select(deriv_mean_med,deriv_sd_med) 
agropast <- derivAll %>% filter(proxy=="agropastoralism") %>% select(deriv_mean_med,deriv_sd_med) 
pollen <- derivAll %>% filter(proxy=="pollen") %>% select(deriv_mean_med,deriv_sd_med) 

# agropastoralism ~ diatoms
df <- cbind(diat,agropast)
colnames(df) <- c("diat_deriv_mean", "diat_sd_mean", "agropast_deriv_mean", "agropast_sd_mean")

#fitting GLS model with autocorrelated residuals
diat.agropast.gls <- gls(diat_deriv_mean ~ agropast_deriv_mean, data=df)
summary(diat.agropast.gls)

#autocorrelation of the residuals
acf(residuals(diat.agropast.gls)) # evidence for correlated residuals
pacf(residuals(diat.agropast.gls)) # indicates AR1

#fitting GLS model with autocorrelated residuals
diat.agropast.gls <- gls(diat_deriv_mean ~ agropast_deriv_mean, 
                         data=df, correlation=corARMA(p=1))

#pseudo R2 (gls doesn't provide R2)
diat.agropast.gls.R2 <- cor(df$diat_deriv_mean, 
                          predict(diat.agropast.gls))^2 ##0.42

df$diat.agrost.predicted <- predict(diat.agropast.gls)

#plotting diatom vs agropastoralism
diatoms_agropast_plt <- ggplot(data=df, aes(x=agropast_deriv_mean, y=diat_deriv_mean)) + 
  geom_point(shape=21, fill="gray50", color="black", size=4, alpha=0.5) +
  geom_line(aes(x=agropast_deriv_mean, y=diat.agrost.predicted), size=2, color="red4", alpha=0.6) +
  ylab("Diatoms (st. dev of rate of change of log-abundance)") + 
  xlab("Agropastoralism-indicators (st. dev of rate of change of log-abundance)") +
  ggtitle("Diatoms vs. Agropastoralism") +
  theme(text=element_text(size=9), plot.title=element_text(size = 12))

diatoms_agropast_plt

# pollen ~ diatoms
df <- cbind(diat,pollen)
colnames(df) <- c("diat_deriv_mean", "diat_sd_mean", "pollen_deriv_mean", "pollen_sd_mean")

#fitting GLS model
diat.pollen.gls <- gls(diat_deriv_mean ~ pollen_deriv_mean, 
                         data=df)

#pseudo R2 (gls doesn't provide R2)
diat.pollen.gls.R2 <- cor(df$diat_deriv_mean, 
                            predict(diat.pollen.gls))^2 #0.02

#adding predicted and residual values to erica.char
df$diat.pollen.predicted <- predict(diat.pollen.gls)

#plotting diatom vs pollen
diatoms_pollen_plt <- ggplot(data=df, aes(x=pollen_deriv_mean, y=diat_deriv_mean)) + 
  geom_point(shape=21, fill="gray50", color="black", size=4, alpha=0.5) +
  geom_line(aes(x=pollen_deriv_mean, y=diat.pollen.predicted), size=2, color="red4", alpha=0.6) +
  ylab("Diatoms (st. dev of rate of change of log-abundance)") + 
  xlab("Pollen (st. dev of rate of change of log-abundance)") +
  ggtitle("Diatoms vs. Pollen") +
  theme(text=element_text(size=9), plot.title=element_text(size = 12))

diatoms_pollen_plt


# pollen ~ agropastoralism
df <- cbind(pollen,agropast)
colnames(df) <- c("pollen_deriv_mean", "pollen_sd_mean", "agropast_deriv_mean", "agropast_sd_mean")

#fitting GLS model
pollen.agropast.gls <- gls(pollen_deriv_mean ~ agropast_deriv_mean, 
                       data=df)

#pseudo R2 (gls doesn't provide R2)
pollen.agropast.gls.R2 <- cor(df$pollen_deriv_mean, 
                          predict(pollen.agropast.gls))^2 #0.004

#adding predicted and residual values to erica.char
df$pollen.agropast.predicted <- predict(pollen.agropast.gls)

#plotting pollen vs agropastoralism
pollen_agropastoralism_plt <- ggplot(data=df, aes(x=pollen_deriv_mean, y=agropast_deriv_mean)) + 
  geom_point(shape=21, fill="gray50", color="black", size=4, alpha=0.5) +
  geom_line(aes(x=agropast_deriv_mean, y=pollen.agropast.predicted), size=2, color="red4", alpha=0.6) +
  ylab("Agropastoralism-indicators (st. dev of rate of change of log-abundance)")+
  xlab("Pollen (st. dev of rate of change of log-abundance)") +
  ggtitle("Pollen vs. Agropastoralism") +
  theme(text=element_text(size=9), plot.title=element_text(size = 12))

pollen_agropastoralism_plt

# plot all three GLS
all_three_mean_deriv <- plot_grid(diatoms_agropast_plt, diatoms_pollen_plt, pollen_agropastoralism_plt,
                       ncol=3)

all_three_sd_deriv <- plot_grid(diatoms_agropast_plt, diatoms_pollen_plt, pollen_agropastoralism_plt,
                                  ncol=3)


################################
## Generate time-lagged datasets
################################

source("scripts/functions_custom.R")
#source("scripts/functions copy.R")

diat <- derivAll %>% filter(proxy=="diatoms") %>% select(negAge,deriv_mean_med) %>%
  mutate(age=negAge) %>% data.frame()
diat$age <- round(diat$age, 1)
diat$negAge <- NULL
#colnames(diat) <- c("ericaceae.par", "age")
colnames(diat) <- c("diatom_deriv", "age")
#diat <- diat[diat$age >= 0,]
round(diff(diat$age, 1))[1]


pollen <- derivAll %>% filter(proxy=="pollen") %>% select(negAge, deriv_mean_med) %>%
  mutate(age=negAge) %>% data.frame()
pollen$age <- round(pollen$age, 1)
pollen$negAge <- NULL
#colnames(pollen) <- c("charcoal.acc.rate", "age") 
colnames(pollen) <- c("pollen_deriv", "age")
#pollen <- pollen[pollen$age >= 0,]
round(diff(pollen$age, 1))[1] #check temporal resolution

agropastolarism <- derivAll %>% filter(proxy=="agropastoralism") %>% select(negAge, deriv_mean_med) %>%
  mutate(age=negAge) %>% data.frame()
agropastolarism$age <- round(agropastolarism$age, 1)
agropastolarism$negAge <- NULL
colnames(agropastolarism) <- c("pollen_deriv", "age") #note the change of the column name
#pollen <- pollen[pollen$age >= 0,]
round(diff(agropastolarism$age, 1))[1]

### Backward lags
# lags * temp resolution = 40 * 24 = 960
lags<-1:40

#backward dataset 
#to assess the effect of “past” pollen and agropastoralism values on diatoms
lag.data.backward <- backwardLags(
  lags=lags, 
  reference.data=diat, 
  data.to.lag=agropastolarism
)


#preparing plotting data (4 lags only)
temp.backward <- lag.data.backward
temp.backward$lag <- as.numeric(temp.backward$lag)
temp.backward <- temp.backward[temp.backward$lag%in% c(1, 3, 6, 10, 40),]

plot.past <- ggplot(data=temp.backward, aes(x=pollen_deriv, y=diatom_deriv, group=lag)) + 
  geom_point(shape=21, fill="gray50", color="black", size=2, alpha=0.5) +
  facet_wrap("lag", ncol=1) +
  xlab("Agropastoralism (avg rate of change)") +
  ylab("Diatoms (avg rate of change)") +
  ggtitle("Backward (Diatom ~ Pollen)") +
  theme(text=element_text(size=12),
        plot.title=element_text(size = 16), legend.position="none") +
  geom_smooth(method = lm, size=2, color="red4", se=FALSE, aes(alpha=0.5))

plot.past


#fitting a GLS model per lag on backward datasets
backward.results <- modelLagData(
  model.formula="diatom_deriv ~ pollen_deriv", 
  lagged.data=lag.data.backward
)

backward.results$value <- round(backward.results$value, 2)

# fitting a null model
backward.results.random <- modelRandomLagData(
  lagged.data=lag.data.backward, 
  model.formula="diatom_deriv ~ pollen_deriv", 
  iterations=1000
)

backward.results.random$value <- round(backward.results.random$value, 2)


# plot results
# plotModelOutput(backward.results=backward.results, 
#                 backward.results.random=backward.results.random, 
#                 width=9, 
#                 height=6, 
#                 title.size=16, 
#                 text.size=10)

                #filename="outputs/Figure2.pdf")

# plot model results with own code
#axes limits
max.lag = max(c(backward.results$lag))
max.coefficient = round(max(c(backward.results[backward.results$variable=="Coefficient", "value"], backward.results.random[backward.results.random$variable=="Coefficient", "upper"])) + 0.1, 1)
min.coefficient = round(min(c(backward.results[backward.results$variable=="Coefficient", "value"], backward.results.random[backward.results.random$variable=="Coefficient", "lower"])) - 0.1, 1)
max.R2 = round(max(c(backward.results[backward.results$variable=="R2", "value"])), 1)

library(viridis)
viridis.colors <- viridis(10, option="D")

backward.plot.coefficient <- ggplot(data=subset(backward.results, variable=="Coefficient"), aes(x=lag, y=value)) +
  geom_ribbon(data=subset(backward.results.random, variable=="Coefficient"), aes(ymin=lower, ymax=upper), alpha=0.3, fill="light grey") +
  geom_line(data=subset(backward.results.random, variable=="Coefficient"), aes(x=lag, y=value), alpha=0.6, color="light grey", size=1) +
  geom_hline(yintercept=0, color="black", linetype=2) +
  geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.3, fill=viridis.colors[2]) +
  geom_line(size=1.5, color=viridis.colors[1]) +
  ggtitle(expression("Agropastoralism" %->% "Diatoms")) +
  theme(legend.position="none") +
  xlab("") +
  ylab("Standardized coefficient") +
  scale_y_continuous(breaks=seq(min.coefficient, max.coefficient, by=0.8)) +
  scale_x_reverse()+
  theme(axis.text = element_text(size=12),
        plot.title = element_text(size=14),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x = element_blank()) +
  theme_classic()
  #coord_cartesian(ylim = c(min.coefficient, max.coefficient + 0.5))

backward.plot.coefficient


backward.plot.R2 <- ggplot(data=subset(backward.results, variable=="R2"), aes(x=lag, y=value, group=variable)) +
  geom_ribbon(data=subset(backward.results.random, variable=="R2"), aes(ymin=lower, ymax=upper), alpha=0.3, fill="light grey") +
  geom_line(data=subset(backward.results.random, variable=="R2"), aes(x=lag, y=value), alpha=0.6, color="light grey", size=1.5) +
  geom_line(size=1.5, color=viridis.colors[2]) +
  theme(legend.position="none") +
  xlab("Years (before Diatom samples)") +
  ylab("Pseudo R squared") +
  scale_y_continuous(breaks=seq(0, max.R2, by=0.1)) +
  scale_x_reverse()+
  theme(axis.text = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(size = 14),
        plot.margin = unit(c(0.2, 0.5, 0, 0), "cm")) +
  theme_classic() +
  coord_cartesian(ylim = c(0, max.R2 + 0.05))

backward.plot.R2

library(cowplot)
plot_composite <- plot_grid(backward.plot.coefficient, 
                            backward.plot.R2, ncol = 1, 
                            rel_heights = c(1, 1), align="v") + 
                            theme(plot.margin = unit(c(0.5, -1, 0.5, 0.5), "cm"))
plot_composite

# save plot
ggsave("outputs/agropastoralism_diatoms_asycnchronousModel.png",
       plot = plot_composite,
       width=8,
       height=6,
       units="in",
       dpi = 400)


##################
#### MUTI ########
##################

library(muti)

derivAll <- read.csv("outputs/diatom_pollen_derivatives.csv", row.names = 1)
derivAll <- read.csv("outputs/pollen-agropast-diatom-deriv_24timeSteps.csv", row.names = 1)

## 

mean <- "deriv_mean_med"
sd <- "deriv_sd_med"

#extract datasets to compute muti among mean derivative values across all lake pairs
diat <- derivAll %>% filter(proxy=="diatoms") %>% select(deriv_mean_med,deriv_sd_med) 
pollen <- derivAll %>% filter(proxy=="pollen") %>% select(deriv_mean_med,deriv_sd_med) 
agropast <- derivAll %>% filter(proxy=="agropastoralism") %>% select(deriv_mean_med,deriv_sd_med) 

df <- cbind(diat,pollen,agropast)
colnames(df) <- c("diat_deriv_mean", "diat_sd_mean", "pollen_deriv_mean", "pollen_sd_mean", "agropast_deriv_mean", "agropast_sd_mean")

df_mean <- df[,c(1,3,5)]
df_sd <- df[,c(2,4,6)]


## Run the Mutual information analysis for mean derivatives
syncDeriv_mean <- list()
counter=0

for (i in 2:ncol(df_mean)) {
  counter = counter + 1
  x = df_mean[,1] #this is the diatom time series (response)
  y = df_mean[[i]] # these are pollen and agropastoralism
  syncDeriv_mean[[counter]] <- muti(x, y, sym = TRUE, normal=TRUE, mc=100, alpha=.05)
}

names(syncDeriv_mean) <- colnames(df_mean[,2:ncol(df_mean)])

## Run the Mutual information analysis for sd derivatiives
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
muti_deriv_mean <- filter(muti_deriv_mean, lag<=0) #filter lags (negative to ask if time-delayed affects of pollen and agropastoralism changes on diatom changes?)

#extract sd_mean lists
muti_sd_mean <- plyr::ldply(syncSd_mean, data.frame)
colnames(muti_sd_mean) <- c("var", "lag",  "MI", "MI_null")
muti_sd_mean <- mutate(muti_sd_mean,sig=ifelse(MI>MI_null, 1, 0)) #is MI significant?
muti_sd_mean <- filter(muti_sd_mean, lag<=0) #filter lags (negative to ask if time-delayed affects of pollen and agropastoralism changes on diatom changes?)

