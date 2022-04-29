
#Load packages for functions used
library(maps)
library(rworldmap)
library(rnaturalearth)
library(raster)
library(cowplot)
library(magick)
library(ggrepel)
library(cividis)
library(ggsn)
library(tidyverse)
library(viridis)
library(png)

world_df <- map_data("world")
interest <- c("Ecuador")
ecuador <- world_df %>% filter(str_detect(region, interest))

world <- ne_countries(scale = "medium", returnclass = "sf")

southamerica <- ggplot(data = world) +
  geom_sf() +
  geom_polygon(data=ecuador, aes(x=long, y=lat, group=group), fill="blue") +
  coord_sf(ylim=c(-55,15), xlim=c(-85,-35)) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
  axis.text.x=element_blank(),
  axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y=element_blank(),
  axis.ticks.x=element_blank())
southamerica
ggsave("outputs/map_figure1_southamerica.png", southamerica,  height = 8, width = 10, dpi = 300)


#Plot elevation base map (raster library)
#unzip(zipfile = "data/dem.zip", exdir = "data/DEM")
DEM <- raster("data/DEM/dem2.bil")
ext_ecuador<-extent(-82,-76,-5,2.5)
altmod_ecuador <- crop(DEM,ext_ecuador)

#convert the raster to points for plotting
map.p_ecuador <- rasterToPoints(altmod_ecuador)

#Make the points a dataframe for ggplot
df_ecuador <- data.frame(map.p_ecuador)

#Make appropriate column headings
colnames(df_ecuador) <- c("long", "lat", "Elevation")

paramo <- df_ecuador %>% filter(between(Elevation, 3500,3600))
forest <- df_ecuador %>% filter(between(Elevation, 2500,3499))

# paramo_raster <- rasterFromXYZ(paramo[,c(1,2,3)])
# plot(paramo_raster)

# create a df with Lake Llaviucu location, Tomebamba and Paredones
sites <- data.frame(long=c(-79.146823, -79.674633,-79.012812),
                  lat=c(-2.8430, -2.775944,-2.902813),
                  label=c("Lake Llaviucu", "Paredones", "Tomebamba (Cuenca)"))


plt_ecuador <- ggplot(data=df_ecuador, aes(y=lat, x=long))+
  geom_raster(aes(fill=Elevation)) +
  scale_fill_cividis()+
  geom_tile(data=paramo, aes(x=long,y=lat,fill=Elevation), colour="black", alpha=0.2) +
  geom_tile(data=forest, aes(x=long,y=lat,fill=Elevation), colour="forestgreen", alpha=0.3) +
  geom_point(data=sites %>% slice(1), aes(long, lat), shape=19, colour="red",size=4)+
  geom_point(data=sites %>% slice(2,3), aes(long, lat), shape=21, bg="grey", colour="black",size=4)+
  geom_text_repel(data=sites, aes(long, lat, label = label), col="white", box.padding = 0.8,
                  nudge_x = 0.20,
                  nudge_y = 0.6) +
  xlab("Longitude (degrees)")+
  ylab("Latitude (degrees)")+
  theme_bw() +
  theme(# panel.background = element_blank(),
        panel.background = element_rect(color = "grey", fill="white"),
        panel.grid.major = element_line(size=0.1, linetype = "solid", color = "grey"),
        panel.grid.minor = element_line(size=0.1, linetype = "solid", color = "grey"),
        legend.position = "bottom")+
  north(df_ecuador, symbol = 4) + #add north arrow
  scalebar(df_ecuador, dist = 40, dist_unit = "km",
           transform = TRUE, model = "WGS84", st.dist = 0.03, st.size = 3)+
  guides(fill = guide_colourbar(title="Elevation (m)")) +
  annotate(geom="text", x=-80, y=2, label="Pacific Ocean",
           color="black", size=4) +
  coord_equal()
plt_ecuador  
ggsave("outputs/map_figure1_ecuador.png", plt_ecuador,  height = 10, width = 12, dpi = 300)


# composite plot
llaviucu_orto <- readPNG("data/llaviucu_orto.png")
llaviucu_bathymetry <- readPNG("data/llaviucu_bathymetry.png")

composite <- ggdraw() +
  draw_plot(southamerica, x=0.17, y=0.63, width=0.2, height=.4, scale = 0.5) +
  draw_plot(plt_ecuador, x=0.01, y=0.03, scale = 0.8) +
  draw_image(llaviucu_orto, x = 0.55, y = 0.47, width=0.6, height=0.7, scale=0.4) +
  draw_image(llaviucu_bathymetry,  x = 0.55, y = 0.01, width=0.6, height=0.7, scale=0.4) +
  draw_plot_label(c("A", "B", "C", "D"), c(0.20, 0.34, 0.71,0.71), c(0.93, 0.93, 0.93, 0.47), size = 15)
composite
dev.off()

#save study area map (Figure 1)
ggsave("outputs/map_figure1.png", composite,  height = 10, width = 12, dpi = 300)


