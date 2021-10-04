######################
#### STRATIPLOTS #####
######################


##loading libraries for functions used
library(analogue) #to join diatom datasets on their common spp
library(rioja) #to merge diatom datasets on their common spp
library(plyr) #allow to join dataframes by common column
library(dplyr) #allow to summarise variables and manipulate multiple dataframes
library(ggplot2) #to make nice plots
library(tidyverse)
library(cluster)

## Prepare Multiproxy dataset for binning and then stratiplotting
# Read relative abundance of diatoms, pollen and agropastoralism
diatRA <- read.csv("data/diatomsRA.csv")[,-1] %>%
  select(-lake, -depth, -lower_age) 
diatRA <- diatRA[order(diatRA$upper_age),] #order time
diatRA_ages <- diatRA[,"upper_age"]
diatRA <- diatRA[,!names(diatRA) %in% "upper_age"]
diatRA <- diatRA[,!names(diatRA) %in% "X"]

pollenRA <- read.csv("data/pollenRA.csv")[,-1] %>% 
  select(-depth, -lower_age)
pollenRA_ages <- pollenRA[,"upper_age"]
pollenRA <- pollenRA[,!(names(pollenRA) %in% "upper_age")]

agropastoralismRA <- read.csv("data/agropastoralismRA.csv")[,-1] %>% 
  gather(key = taxa, value = relative_abundance_percent, -depth, -upper_age, -lower_age)
core_counts_common <- agropastoralismRA
core_counts_common <- pollen_data


  #make it wide
  core_counts_wide <- core_counts_common %>%
    dplyr::select(depth, upper_age, taxa, relative_abundance_percent) %>%
    spread(key = taxa, value = relative_abundance_percent)



# make name vectors to later replace taxon names to assemblage group
nms_list <- data.frame(taxa = c(colnames(diatRA), colnames(pollenRA), colnames(agropastoralismRA)), 
                       group = c(rep("diat",length(diatRA)), rep("pollen",length(pollenRA)), 
                                 rep("agropastoralism",length(agropastoralismRA))))


#make the age categories (60-years bins; median combined age interval) and combine the two datasets
PollenBinned <- binFunc(as.data.frame(pollenRA), as.numeric(pollenRA_ages), 60, 0, 3000) 
AgropastBinned <- binFunc(as.data.frame(agropastoralismRA), as.numeric(pollenRA_ages), 60, 0, 3000)
DiatBinned <- binFunc(as.data.frame(diatRA), as.numeric(diatRA_ages), 60, 0, 3000) 

# cbind dataframes
combinedData <- cbind(PollenBinned, AgropastBinned, DiatBinned)
#varBinInter <- na.approx(combinedData, na.rm = TRUE) #do interpolation between adjacent samples

#gather dataset for plotting
df_long <- combinedData %>% 
  as.data.frame() %>%
  mutate(upper_age=row.names(combinedData)) %>%
  mutate(upper_age=as.numeric(upper_age)) %>%
  gather(key=taxa, value=relative_abundance_percent, -upper_age) %>%
  mutate(assemblage = plyr::mapvalues(taxa, from = nms_list$taxa, to = nms_list$group)) %>%
  mutate(taxa=factor(taxa)) %>%
  #filter(relative_abundance_percent>2) %>%
  arrange(desc(relative_abundance_percent))

levels(df_long$taxa)

## Using tidypaleo R package (https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
theme_set(theme_bw(9))

composite_plot <- ggplot(df_long, aes(x = relative_abundance_percent, y = upper_age)) +
  geom_col_segsh() +
  scale_y_reverse() +
  facet_abundanceh(vars(taxa)) +
  labs(x = "Relative abundance (%)", y = "Cal yr BP")
composite_plot


#read diatom data
mergedCores <- read.csv("data/mergedCores_counts4.csv")[,-1] #with new Fondococha agedepth model
diatoms_save <- mergedCores #save dataframe

#Gather
spp_thin <- diatoms_save %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age,  -lake)#don't gather depths, ages and lake variables


#import dataframe wiht old and new names to group
changes <- read.csv("data/old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
#new1: ecological groups
#new2: harmonized taxonomic names


#spread--> wide format
spp_wide <- spp_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth, lake, upper_age, taxa) %>%
  summarise(count = sum(count)) %>%
  spread(key = taxa, value = count)

#no spread --> long format
spp_long <- spp_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth, lake, upper_age, taxa) %>%
  summarise(count = sum(count))


#filter cores
select <- c("llaviucu")

core_lake <- spp_long %>%
  filter(str_detect(lake, select)) %>% #select lake
  filter(!upper_age==0.0) %>%
  group_by(taxa) %>%
  filter(count > 0) %>% #remove species with 0 counts
  ungroup() %>%
  group_by(depth) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>% #calculate RA
  ungroup()

# filter more abundant taxa; format need to be on long-wide format-->no spreaded 
core_common_taxa <- core_lake %>%
  group_by(taxa) %>%
  summarise(max_rel_abund = max(relative_abundance_percent)) %>%
  filter(max_rel_abund >= 5) %>%
  arrange(max_rel_abund) %>%
  pull(taxa)

# select from initial table
core_counts_common <- core_lake %>%
  filter(taxa %in% core_common_taxa) %>%
  mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
  arrange(taxa)

#make it wide
core_counts_wide <- core_counts_common %>%
  dplyr::select(depth, lake, upper_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent)

length(core_common_taxa)

#do coniss
core_counts_wide[is.na(core_counts_wide)] <- 0

diatHel <- decostand(core_counts_wide[,4:ncol(core_counts_wide)], method="hellinger")
diss <- vegdist(diatHel, method="bray")
clust <- chclust(diss, method="coniss")
bstick(clust)

zones <- cutree(clust, k=4)
locate <- cumsum(rle(zones)$lengths)+1
zones <- core_counts_wide[locate, ][,3]
zones <- zones$upper_age


library(viridis)
seq_palette <- viridis(5)

gg <- ggplot(core_counts_common, aes(x=as.numeric(as.character(upper_age)), y=relative_abundance_percent))
plt <- gg + geom_area(aes(colour=taxa, fill=taxa)) +
  ylab("% total assemblage") +
  xlab("Years cal BP")+
  theme_bw()
plt

############################## END

## plot using analogue Stratiplot (need data input in wide format --> spreaded)
png("Stratiplot.png", width = 11, height = 8, res = 300, units = "in")

Stratiplot(
  core_counts_wide %>% select(-depth, -upper_age),
  core_counts_wide$upper_age,
  ylab = "Cal yr BP", 
  xlab = "Relative abundance (%)",
  # adds padding to the top of the plot
  # to fix cut-off taxa names
  topPad = 10, 
  # make the plot type a "bar" plot
  type = "h", 
  #sort = "wa",
  # add stratigraphic zones from cluster analyses (regime shifts R file)
  #zones = zones,
  # make the bar colour black
  col = "black")

dev.off()  



## merge species and non-species data  
#Use a left-join to add non-species data
llaviucu_xrf <- read.csv("data/llaviucu_xrf.csv")

#read Llaviucu pollen
llaviucu_pollen_ratios <- read.csv("data/llaviucu_pollen.csv") %>% 
  gather(key = taxa, value = relative_abundance_percent, -depth, -age) %>%
  mutate(upper_age=age)

  core_counts_common <- llaviucu_pollen_ratios
  
  #make it long for joining with non-species data
  core_counts_long <- core_counts_common %>%
    select(depth, upper_age, taxa, relative_abundance_percent) 

llaviucu_pollen <- read.csv("data/llaviucu_pollen_raw.csv")

#make it long for joining with non-species data
core_counts_long <- core_counts_common %>%
  select(depth, lake, upper_age, taxa, relative_abundance_percent) 


#llaviucu
core_diat_geochem <- core_counts_long %>%
  left_join(llaviucu_xrf, by = "depth") %>%
  mutate(Mn_Fe = Mn/Fe) %>%
  mutate(Si_Ti = Si/Ti) %>%
  mutate(K_Ti = K/Ti) %>%
  select(taxa, depth, everything())

#table for tidypaleo stratiplot
geochem_data <- core_diat_geochem %>% select(Mn_Fe,Si_Ti,K_Ti,upper_age) %>%
  gather(key=variable,value=value,-upper_age) 

core_diat_geochem_wide <- core_diat_geochem %>%
  select(Mn_Fe, Si_Ti, K_Ti, depth, upper_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent)

#write.csv(core_diat_geochem_wide, "llaviucu_diat_proxy.csv")

png("Stratiplot.png", width = 11, height = 8, res = 300, units = "in")

#varTypes argument needs to specify that the non-species variables should have independently sized axes. 
#This should be a vector with the same number of elements as variables in the plot (I’ve use rep() to repeat “relative” and “absolute” the correct number of times.

# code for diat spp and geochemistry
Stratiplot(
  core_diat_geochem_wide %>% select(-depth, -upper_age),
  core_diat_geochem_wide$upper_age, 
  varTypes = c(rep("absolute", 3), rep("absolute", length(core_common_taxa))), 
  ylab = "Age (cal years BP)", 
  xlab = "Relative abundance (%)",
  topPad = 10, 
  type = c("h"),
  col = "black", 
  #zones = zones
)

dev.off()


## Using tidypaleo R package (https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
library(tidypaleo)
library(patchwork)
theme_set(theme_bw(9))

diat_plot <- ggplot(core_counts_common, aes(x = relative_abundance_percent, y = upper_age)) +
  geom_col_segsh() +
  #geom_lineh() +
  #geom_areah() +
  scale_y_reverse() +
  facet_abundanceh(vars(taxa)) +
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2) +
  labs(x = "Relative abundance (%)", y = "Cal years BP")
#add tephra layers at certain depths
#annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=2082, ymax=2117, alpha=0.7, fill="grey") +
# add CONISS zones
#geom_hline(yintercept = zones, col = "blue", lty = 1, alpha = 0.7) 
diat_plot

ggsave("outputs/pollen_stratplot.png", diat_plot, height = 6, width = 10)


## prepare geochemical data 
# FIRST go to lines 264
geochem_plot <- ggplot(geochem_data, aes(x = value, y = upper_age)) +
  geom_lineh() +
  #geom_point() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(variable)) +
  labs(x = NULL) 

geochem_plot <- geochem_plot +
  geom_lineh_exaggerate(exaggerate_x = 4, col = "grey70", lty = 1)

strat_plt <- wrap_plots(
  diat_plot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  geochem_plot +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(4, 1)
)
strat_plt

ggsave("outputs/llaviucu_stratplot.png", strat_plt, height = 6, width = 10)




