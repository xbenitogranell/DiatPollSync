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
  select(-depth, -lower_age) 
agropastoralismRA <- agropastoralismRA[order(agropastoralismRA$upper_age),]  #order time
agropastoralismRA_ages <- agropastoralismRA[,"upper_age"]
agropastoralismRA <- agropastoralismRA[,!names(agropastoralismRA) %in% "upper_age"]
agropastoralismRA <- agropastoralismRA[,!names(agropastoralismRA) %in% "X"]


# make name vectors to later replace taxon names to assemblage group
nms_list <- data.frame(taxa = c(colnames(diatRA), colnames(pollenRA), colnames(agropastoralismRA)), 
                       group = c(rep("Diatoms",length(diatRA)), rep("Pollen",length(pollenRA)), 
                                 rep("Agropastoralism",length(agropastoralismRA))))


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
  mutate(assemblage=recode(assemblage, 
                    '1'="Agropastoralism",
                    '2'="Diatoms",
                    '3'="Pollen")) %>%
  mutate(taxa=factor(taxa),
         taxa = fct_reorder(taxa, relative_abundance_percent)) %>%
  filter(relative_abundance_percent>4) %>%
  mutate(taxa=reorder(taxa,relative_abundance_percent)) %>%
  group_by(assemblage,taxa) %>%
  arrange(relative_abundance_percent) %>%
  ungroup() %>%
  mutate(taxa = fct_relevel(taxa, "Sporormiella","Cyperaceae","Asteraceae.1","Hedyosmum","Alnus",
                            "Navicula.radiosa","Aulacoseira.alpigena","Cocconeis.placentula.var.placentula","Encyonopsis.subminuta","Hannaea.arcus","Orthoseria.roseana","Ulnaria.cf.ulna.var.amphyrhynchus",
                            "Gomphonema.sp.1.LLAVIUCU.cf.netriviale","Gomphonema.cf.geisslerae", "Gomphonema.sp.2.LLAVIUCU.cf.variostriatum", "Denticula.kuetzingii","Cymbella.cymbiformis.cystula","Encyonopsis.sp.1.LLAVIUCU","Discostella.stelligera",
                            "Fragilaria.cf.capucina","Diatoma.tenuis","Brachysira.microcephala.neoxilis","Staurosira.construens.var.venter","Fragilaria.capuccina.agg","Tabellaria.flocculosa.str.IV","Nupela.sp.1.LLAVIUCU",
                            "Achnanthidium.minutissimum", "Cyathea","Clethra","Symplocos","Polylepis","Weinmannia","Acalypha","Isoetes","Melast.Combret","Podocarpus","Morac.Urtica",
                            "Monolete.psilate","Poaceae")) %>%
  drop_na(relative_abundance_percent)

levels(df_long$taxa)

## Using tidypaleo R package (https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
theme_set(theme_bw(12))
#theme_set(theme_paleo())

assemblage_plot <- ggplot(df_long, aes(x = relative_abundance_percent, y = upper_age, colour=assemblage)) +
  geom_col_segsh(size=1.3) +
  scale_y_reverse() +
  facet_abundanceh(vars(taxa), grouping = vars(assemblage), rotate_facet_labels = 70) +
  #geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2) +
  labs(x = "Relative abundance (%)", y = "cal years BP", colour="Assemblage") +
  #ggtitle("Assemblages") +
  theme (legend.position = "none") 
assemblage_plot

#ggsave("outputs/proxies_stratplot.png", diat_plot, height = 6, width = 10)

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
  mutate(taxa = plyr::mapvalues(taxa, from = changes[,1], to = changes$new_2)) %>%
  group_by(depth, lake, upper_age, taxa) %>%
  summarise(count = sum(count)) %>%
  spread(key = taxa, value = count)

#no spread --> long format
spp_long <- spp_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes[,1], to = changes$new_2)) %>%
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


## merge species and non-species data  
#Use a left-join to add non-species data
llaviucu_xrf <- read.csv("data/llaviucu_xrf.csv")

#read Llaviucu pollen
llaviucu_pollen_ratios <- read.csv("data/llaviucu_pollen.csv") %>% 
  gather(key = taxa, value = relative_abundance_percent, -ï..depth, -age) %>%
  mutate(upper_age=age)

#make it long for joining with non-species data
core_counts_long <- core_counts_common %>%
  select(depth, lake, upper_age, taxa, relative_abundance_percent) 


#llaviucu
core_diat_geochem <- core_counts_long %>%
  left_join(llaviucu_xrf, by = "depth") %>%
  mutate(MnFe = Mn/Fe) %>%
  mutate(SiTi = Si/Ti) %>%
  #mutate(K_Ti = K/Ti) %>%
  select(taxa, depth, everything())

#table for tidypaleo stratiplot
geochem_data <- core_diat_geochem %>% select(MnFe,SiTi,upper_age) %>%
  gather(key=variable,value=value,-upper_age) 

core_diat_geochem_wide <- core_diat_geochem %>%
  select(Mn_Fe, Si_Ti, depth, upper_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent)

# Plot geochemistry using tidypaleo R package
theme_set(theme_bw(12))

geochem_plot <- ggplot(geochem_data, aes(x = value, y = upper_age)) +
  geom_lineh() +
  #geom_line(data = m.pred, aes(x=logMnFe, y=age)) + #first I pre-calculate model fit and make predictions
  #geom_smooth()+
  scale_y_reverse() +
  #geom_point() +
  facet_geochem_gridh(vars(variable)) +
  labs(x = NULL) 
  #ggtitle("XRF")
geochem_plot


## combine plots
library(patchwork)
strat_plt2 <- wrap_plots(
  assemblage_plot + 
    theme(strip.text.y = element_text(size=10),
          strip.text.x = element_text(size=9)),
          #panel.grid = element_blank()),
    geochem_plot +
    theme(strip.text.x = element_text(size = 10),
          axis.ticks.y.left = element_blank(),
          axis.text.y.right = element_text()) +
    labs(y = "cal years BP"),
  nrow = 1,
  widths = c(4, 0.5)
)
strat_plt2

ggsave("outputs/proxies_stratplot_REV.png", strat_plt2, height = 9, width = 11, dpi=150)





