source("scripts/functions_custom.R") #functions to lag datasets and model (Blas) and binning (Seddon)

library(tidypaleo) #remotes::install_github("paleolimbot/tidypaleo")
library(patchwork)
library(zoo)
library(vegan)



## Make it wide for multivariate ordination
df_wide <- df_long %>% spread(key = taxa, value = relative_abundance_percent) 

age_assemblage <- df_wide[,names(df_wide) %in% c("upper_age","assemblage")]
df_wide <- df_wide[,!names(df_wide) %in% c("upper_age","assemblage")]

df_wide[is.na(df_wide)] <- 0

df_hellinger <- decostand(df_wide, method = "hellinger")


## Make multivariate ordination
mod <- decorana(df_hellinger)
plot(mod)

#### Multivariate ordination with multiproxy linearly interpolated dataset
interpolatedData <- read.csv("outputs/principalcurves_ti_interp.csv", row.names=1)

# Prepare data for NMDS
age <- interpolatedData[,names(interpolatedData) %in% c("Age")]
vec <- c(1,2, 101:111) #where NAs are located (Tis)

interpolatedData <- interpolatedData[,!names(interpolatedData) %in% c("Age")]
interpolatedData <- interpolatedData[,!names(interpolatedData) %in% c("elapsedTime_ti")]
interpolatedData <- interpolatedData[,!names(interpolatedData) %in% c("elapsedTime_diat")]
interpolatedData <- interpolatedData[,!names(interpolatedData) %in% c("charcoal")]

interpolatedData2 <- na.omit(interpolatedData)

## Run CONISS
interpl.D <- vegdist(interpolatedData2, "bray")
clust <- chclust(interpl.D, method="coniss")
bstick(clust)
plot(clust)
clustsig <- cutree(clust, k=5)

# Run NMDS
nmds <- metaMDS(interpolatedData2, distance = "jaccard")
plot(nmds)
plot(nmds, display = "sites")
text(nmds, display = "species")

#Extract NMDS site and species scores
scrs_sites <- scores(nmds, display = "sites", choices = 1:2)
scrs_var <- as.data.frame(scores(nmds, display = "species", choices = 1:2))

# Plot
nmds_tbl <- mutate(data.frame(scrs_sites), time=age[-vec])

nmds_plt <- ggplot(nmds_tbl, aes(NMDS1,NMDS2)) + 
  xlab("NMDS1") + ylab("NMDS2") +
  geom_point(aes(colour=time)) +
  geom_text(data=scrs_var, aes(NMDS1, NMDS2, label=rownames(scrs_var)))+
  scale_colour_viridis_c()+
  # geom_vline(aes(xintercept = 0), linetype = "solid", colour="grey") +
  # geom_hline(aes(yintercept = 0), linetype = "solid", colour="grey") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=9))+
  theme_bw()
nmds_plt

## Here create phase space plots
ggplot(data=interpolatedData2, aes(x=agropastPrC, y=diatPrC))+
  geom_point(aes(colour=age[-vec]))+
  scale_colour_viridis_c() +
  theme_bw()

ggplot(data=interpolatedData2, aes(x=agropastPrC, y=pollenPrC))+
  geom_point(aes(colour=age[-vec]))+
  scale_colour_viridis_c() +
  theme_bw()

ggplot(data=interpolatedData2, aes(x=pollenPrC, y=agropastPrC))+
  geom_point(aes(colour=age[-vec]))+
  scale_colour_viridis_c() +
  theme_bw()

ggplot(data=interpolatedData2, aes(x=Ti, y=pollenPrC))+
  geom_point(aes(colour=age[-vec]))+
  scale_colour_viridis_c() +
  theme_bw()

ggplot(data=interpolatedData2, aes(x=Ti, y=agropastPrC))+
  geom_point(aes(colour=age[-vec]))+
  scale_colour_viridis_c() +
  theme_bw()

ggplot(data=interpolatedData2, aes(x=Ti, y=diatPrC))+
  geom_point(aes(colour=age[-vec]))+
  scale_colour_viridis_c() +
  theme_bw()


# DCA
mod_dca<- decorana(interpolatedData2)
mod_dca
summary(mod_dca)

dca.eig <- mod_dca$evals / mod_dca$tot.chi
dca.eig

#plot using ggvegan package
autoplot(mod_dca)

dca_plt <- fortify(mod_dca)

ford <- fortify(dca_plt, axes = 1:2)  # fortify the ordination
take <- c('DCA1', 'DCA2')  # which columns contain the scores we want
species <- subset(ford, Score == 'species')  # take only biplot arrow scores

## multiplier for arrows to scale them to the plot range
mul <- ggvegan:::arrowMul(species[, take],
                          subset(ford, select = take, Score == 'sites'))
species[, take] <- species[, take] * mul  # scale biplot arrows
data_plt <- cbind(ford[1:98,], as.factor(clustsig), age[-vec])
colnames(data_plt)[7] <- c("cluster")
colnames(data_plt)[8] <- c("age")

#plot
dca_plot <- ggplot() +
  geom_point(data = subset(data_plt, Score == 'sites'),
             mapping = aes(x = DCA1, y = DCA2, colour=age)) + 
  geom_segment(data = species,
               mapping = aes(x = 0, y = 0, xend = DCA1, yend = DCA2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  # geom_point(data = subset(data_plt, Score == 'sites'),
  #            aes(x = DCA1, y = DCA2, colour=cluster))+
  # stat_ellipse(data = subset(data_plt, Score == 'sites'), aes(x = DCA1, y = DCA2, colour = cluster))+
  geom_text(data = species, # crudely push labels away arrow heads
            mapping = aes(label = Label, x = DCA1 * 1.1, y = DCA2 * 1.1)) +
  scale_colour_viridis_c()+
  coord_fixed()+
  theme_classic()

ggsave("outputs/DCA_plot.png",
       plot = dca_plot,
       width = 10,
       height=8,
       units="in",
       dpi = 400)


# Extract DCA results
DCA_results <- cbind(data_plt[,c(3:4, 7,8)], interpolatedData2)
cor.r <- corr.test(DCA_results[,1:2], DCA_results[,c(5:ncol(DCA_results))], method = "spearman")$r
                                                 
cor(DCA_results[1:25,1:2])


##
par(mfrow=c(1,2))
plot(mod_dca, display="sites")
plot(mod_dca, display="species", cex=0.6)
par(mfrow=c(1,1))

#Extract CA site and species scores
DCA.scores <- data.frame(DCA1=scores(mod_dca, display = "sites")[,1], 
                         DCA2=scores(mod_dca, display = "sites")[,2], 
                         time=age[-vec],
                         clust=clustsig)

plot(DCA.scores$DCA1, DCA.scores$DCA2, type = "n", xlab = "", ylab = "")
title("Sites")
abline(h=0, col="grey")
abline(v=0, col="grey")

#manual coding for months
segments(DCA.scores$DCA1[-length(DCA.scores$DCA1)],
         DCA.scores$DCA2[-length(DCA.scores$DCA2)],
         DCA.scores$DCA1[-1L],DCA.scores$DCA2[-1L], col=DCA.scores$clust)


lines(DCA.scores[,1:2], col=as.factor(clustsig), lwd = 2)

