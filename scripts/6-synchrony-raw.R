
library(nlme)
library(tidyverse) 
library(analogue)
library(rioja)

source("scripts/functions_custom.R") #functions to lag datasets and model (Blas) and binning (Seddon)

#read in Principal curves data interpolated to diatom time visits
interpolatedData <- read.csv("outputs/principalcurves_ti_interp.csv")

# Binning samples by categories of ages
diff(interpolatedData$Age) #this is the temporal resolution among samples of diatom core

#removing samples with temporal resolution older than 60
#data.reduced <- interpolatedData %>% filter(elapsedTime_diat)

#make the age categories (60-years bins)
dataBinned <- binFunc(as.data.frame(interpolatedData), as.numeric(interpolatedData$Age), 24, 0, 3000) 

#remove NaNs
#dataBinned <- na.omit(dataBinned)
str(dataBinned)

dataBinned$Age <- rownames(dataBinned)

#diat <- dataBinned$diatPrC
#pollen <- dataBinned$pollenPrC
#agropastolarism <- dataBinned$agropastPrC


################################
## Generate time-lagged datasets
################################

diat <- dataBinned %>% select(diatPrC, Age)
#diat$age <- round(diat$Age, 1)
#diat$Age <- NULL
colnames(diat) <- c("diatom_deriv", "age")
diat$age <- as.numeric(diat$age)
diff(diat$age)
round(diff(diat$age, 1))[1]

pollen <- dataBinned %>% select(pollenPrC, Age)
#pollen$age <- round(pollen$Age, 1)
#pollen$Age <- NULL
colnames(pollen) <- c("pollen_deriv", "age")
pollen$age <- as.numeric(pollen$age)
diff(pollen$age)


agropastolarism <- dataBinned %>% select(agropastPrC, Age)
#agropastolarism$age <- round(agropastolarism$Age, 1)
#agropastolarism$Age <- NULL
colnames(agropastolarism) <- c("pollen_deriv", "age")
agropastolarism$age <- as.numeric(agropastolarism$age)
diff(agropastolarism$age)


### Backward lags
# lags * temp resolution = 40 * 24 = 960
lags<-1:40

#backward dataset 
#to assess the effect of “past” pollen or agropastoralism values on diatoms
lag.data.backward <- backwardLags(
  lags=lags, 
  reference.data=diat, 
  data.to.lag=pollen
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


# remove NaNs
lag.data.backward <- na.omit(lag.data.backward)

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
  ggtitle(expression("Pollen" %->% "Diatoms")) +
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

#Combine plots
library(cowplot)
plot_composite <- plot_grid(backward.plot.coefficient, 
                            backward.plot.R2, ncol = 1, 
                            rel_heights = c(1, 1), align="v") + 
  theme(plot.margin = unit(c(0.5, -1, 0.5, 0.5), "cm"))
plot_composite

# save plot
ggsave("outputs/pollen_diatoms_asycnchronousModel_rawdata.png",
       plot = plot_composite,
       width=8,
       height=6,
       units="in",
       dpi = 400)






##################
## POLLEN datasets 
##################

# read data (drop Lycopodium counts)
pollen <- read.csv("data/llaviucu_pollen_counts.csv") %>%
  select(everything(), -contains("Lycopodium"))

#select human disturbance pollen taxa
agropastolarism_indicators <- c("Zea", "Hedyosmum", "Phaseolus", "Ipomoea", "Rumex", 
                                "Alnus", "Cyperaceae", "Cecropia", "Asteracea", "sporormiella")
agropastolarism <- select(pollen, contains(agropastolarism_indicators))
agropastolarism <- data.frame(sapply(agropastolarism, function(x) as.numeric((x))))
agropastolarism[is.na(agropastolarism)] <- 0 #Replace NA (if any) by 0

#drop agropastoralism indicators from pollen dataset
pollen <- pollen %>%  select(everything(), -contains(agropastolarism_indicators))


## extract agedepth variables
agedepth <- pollen[,names(pollen) %in% c("depth", "upper_age", "lower_age")]
pollen <- pollen[,!(names(pollen) %in% c("upper_age", "lower_age", "depth"))]

pollen <- data.frame(sapply(pollen, function(x) as.numeric((x))))
pollen[is.na(pollen)] <- 0 #Replace NA (if any) by 0

## assign dataframe to analyze
llaviucu_pollen <- pollen
#llaviucu_pollen <- agropastolarism

##Calculate relative abundance
total <- apply(llaviucu_pollen, 1, sum)
llaviucu_pollen <- llaviucu_pollen / total * 100

#Select spp 
abund <- apply(llaviucu_pollen, 2, max)
n.occur <- apply(llaviucu_pollen>0, 2, sum)
llaviucu_pollen <- llaviucu_pollen[, abund>1 & n.occur >2] #more than 3% of RA and present in >2 samples


#dissimilarity rate of change
roc2 <- as.matrix(paldist(llaviucu_pollen, dist.method="sq.chord"))
roc2 <- roc2[row(roc2) == col(roc2) + 1] ## extract off-diagonal

intervals <- diff(agedepth$upper_age)
roc2 <- roc2 / intervals
#roc2 <- roc2*100

plot(roc2 ~ head(agedepth$upper_age, -1), type = "l", ylab = "Rate of Change", xlab = "Age")
points(roc2 ~ head(agedepth$upper_age, -1), type = "h")
roc2

