##loading libraries for functions used
library(analogue) #to join diatom datasets on their common spp
library(rioja) #to merge diatom datasets on their common spp
library(vegan) #to perform multivariate ordination analysis (dca, nmds, rda and varpart)
library(mgcv) #allow to perform GAM analyses
library(tidyverse) #allow to join dataframes by common column
library(ggplot2) #allow to make fancy graphs

theme_set(theme_bw())

#read diatom core datasets
mergedCores <- read.csv("data/mergedCores_counts4.csv") #this is a dataframe with absolute counts containing all the spp

agedepth <- mergedCores[, names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
diat <- mergedCores[, !names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
diat[is.na(diat)] <- 0

diatoms_save <- cbind(agedepth, diat)

changes <- read.csv("data/old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
#new1: ecological groups
#new2: harmonized taxonomic names

#this is to transform to tidy format, calculate % and subset more common species
new <- diatoms_save %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -lake) %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth, taxa, lake, upper_age, lower_age) %>%
  summarise(count = sum(count)) %>%
  filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  filter(!upper_age==0.0) %>% #this is to drop extraneous ages
  ungroup() %>%
  group_by(depth, lake) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  ungroup()

# #this is to calculate planktic:benthic ratios
# new <- diatoms_save %>% 
#   gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -lake) %>%
#   mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_1)) %>% #ecological grouping
#   group_by(depth, taxa, lake, upper_age, lower_age) %>%
#   summarise(count = sum(count)) %>%
#   filter(!count == "0" ) %>% #this is to remove empty samples (rows)
#   filter(!upper_age==0.0) %>% #this is to drop extraneous ages
#   ungroup() %>%
#   group_by(depth, lake) %>%
#   mutate(relative_abundance_percent = count / sum(count) * 100) %>%
#   mutate(plank=sum(count[taxa=="freshwater_planktic" | taxa=="tycoplanktonic"])) %>%
#   mutate(benthic=sum(count[taxa=="epiphytics"| taxa== "saline" | taxa=="benthic"])) %>%
#   mutate(P_B=plank/benthic) %>%
#   mutate(P_B2=(plank-benthic)/(plank+benthic)) %>% #[-1(benthic dominated) to 1(planktic dominated)]
#   ungroup() 

#make it wide
# core_counts_ratios <- new %>%
#   select(depth, lake, upper_age, lower_age, taxa, P_B2, relative_abundance_percent) %>%
#   spread(key = taxa, value = relative_abundance_percent) 
# 
# ## split cores by lakes and reassemble
# coresList <- split(core_counts_ratios, core_counts_ratios$lake)
# saveRDS(coresList, file = "diatomratios.rds")

# filter more abundant taxa; format need to be on long-wide format-->no spreaded 
core_common_taxa <- new %>%
  group_by(taxa) %>%
  summarise(max_rel_abund = max(relative_abundance_percent)) %>%
  filter(max_rel_abund >= 5) %>%
  arrange(max_rel_abund) %>%
  pull(taxa)

# select from initial table
core_counts_common <- new %>%
  filter(taxa %in% core_common_taxa) %>%
  mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
  arrange(taxa)

#make it wide
core_counts_wide <- core_counts_common %>%
  select(depth, lake, upper_age, lower_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent) 

## split cores by lakes and reassemble
coresList <- split(core_counts_wide, core_counts_wide$lake)


# function to calculate rate of change (squared chord distance standardized by age intervals)
SCD <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[ , -which(names(core) %in% c("depth","upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
  core[is.na(core)] <- 0
  core <- core[, colSums(core) > 0] #select only present species
  core <- tran(core, method="hellinger") #Hellinger transform relative abundance data
  scd <- distance(core, method = "SQchord")
  upper_age <- coresList[[i]]$upper_age
  lower_age <- coresList[[i]]$lower_age
  elapsedTime <- abs(upper_age - lower_age)
  age <- cbind(upper_age, lower_age, elapsedTime)[-1,]
  elapsed <- age[,3] #elapsed time nc==3
  scd <- scd[row(scd) == col(scd) + 1]
  #SCDcrop<-scd[-1,]
  #scd<- diag(SCDcrop)
  roc <- scd/elapsed 
  #roc <- roc*100
  cbind.data.frame(age,scd, roc) #combine extracted columns and remove first row to match with scd
  
}

## apply SCD function function to each core
coreSCD <- lapply(seq_along(coresList), SCD, cores=coresList)
names(coreSCD) <- names(coresList)

#extract dataframes from list
roc <- plyr::ldply(coreSCD, data.frame)

## write this out for use in the GAM modelling script
saveRDS(roc, file = "roc.rds") #this is including all spp

## plot the data
ROCPlot <- ggplot(roc, aes(x = upper_age, y = roc)) +
  geom_point(aes(colour=.id)) +
  facet_wrap(~ .id, ncol = 1, scales = "free") + 
  scale_x_reverse()
ROCPlot


# Function to fit principal curves
fitPcurve <- function(i, cores, axis, method, ...) {
  axis <- axis[i]
  core <- cores[[i]]
  method <- method[i]
  core <- core[ , -which(names(core) %in% c("depth","upper_age", "lower_age", "lake"))] # drop year & depths vars
  core[is.na(core)] <- 0
  core <- core[, colSums(core) > 0] #select only present species
  core <- decostand(core, method="hellinger") #Hellinger transform relative abundance data
  prc <- prcurve(core, vary = TRUE, trace = TRUE, method = method,
                 axis = axis)
  prc
}
## apply fitCurve  function to each core
axis <- c(1, 1, 1, 1, 1, 1, 1, 1) 
method <- c(rep("ca", 8)) 

corePRC <- lapply(seq_along(coresList), fitPcurve,
                  cores = coresList, axis = axis, method = method)
names(corePRC) <- names(coresList)

# extract data from list
scrs <- lapply(corePRC, scores)
scrs <- lapply(scrs, `[`, , 1)
len <- lapply(scrs, length)
nams <- names(scrs)
depths <- lapply(coresList, `[`, , "depth")
age <- lapply(coresList, `[`, ,  "upper_age")
upper_age <- lapply(coresList, `[`, ,  "upper_age")
lower_age <- lapply(coresList, `[`, ,  "lower_age")
df <- data.frame(Lake = rep(nams, times = len),
                 PrC  = unlist(scrs),
                 Depth = unlist(depths),
                 Age = unlist(age),
                 negAge = - unlist(age),
                 upper_age = unlist(upper_age),
                 lower_age = unlist(lower_age))
rownames(df) <- NULL

cores.plot <- ggplot(df, aes(x = Age, y = PrC)) +
  geom_point() + geom_line() +
  #coord_flip() +
  facet_wrap(Lake~., scales = "free")
cores.plot

## write this out for use in the GAM modelling script
saveRDS(df, file = "outputs/PrC-cores-diatoms.rds") #this is including all spp

#############
## POLLEN PrC
#############

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
#llaviucu_pollen <- pollen
llaviucu_pollen <- agropastolarism

##Calculate relative abundance
total <- apply(llaviucu_pollen, 1, sum)
llaviucu_pollen <- llaviucu_pollen / total * 100

#Select spp 
abund <- apply(llaviucu_pollen, 2, max)
n.occur <- apply(llaviucu_pollen>0, 2, sum)
llaviucu_pollen <- llaviucu_pollen[, abund>1 & n.occur >2] #more than 3% of RA and present in >2 samples


# Run Principal Curves
core_hell <- decostand(llaviucu_pollen, method="hellinger")
pollen.prc <- prcurve(core_hell, method = "ca", trace = TRUE, vary = TRUE, penalty = 1.4)

## Extract position on the curve
cores_prc <- scores(pollen.prc, display = "curve")

# Combine dataframe with ages and depths
pollenPrC <- cbind(agedepth, cores_prc)

## Plot Pcurves with depth and ages
pollenPlot <- ggplot(pollenPrC, aes(x = upper_age, y = PrC)) +
  geom_line() + geom_point()
pollenPlot

## save results
write.csv(pollenPrC, "outputs/pollen-PrC.csv", row.names = FALSE)
write.csv(pollenPrC, "outputs/agropastoralism-PrC.csv", row.names = FALSE)

