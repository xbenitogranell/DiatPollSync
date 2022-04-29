#age-depth models

library(Bchron)
library(tidyverse)

dates <- read.csv('data/cores_dates_Bchron.csv', sep = ",")
head(dates)

core_2009 <- dates[1:4,]
core_2014 <- dates[5:nrow(dates),]

#core 2014 has outliers, exclude them
core_2014 <- core_2014[c(3,5,9,11),]

core_2009$calCurves <- "intcal13"
core_2014$calCurves <- "intcal20"

# run Bchron with pollen core
core_2009_cal13 <- with(core_2009, 
               Bchronology(ages=ages,
                           ageSds=ageSds, 
                           calCurves=calCurves,
                           positions=position, 
                           positionThicknesses=thickness,
                           ids=id, 
                           predictPositions=seq(0,250,by=10)))

plot(core_2009_cal13)

# run Bchron with diatom core
core_2014_cal13 <- with(core_2014, 
                        Bchronology(ages=ages,
                                    ageSds=ageSds, 
                                    calCurves=calCurves,
                                    positions=position, 
                                    positionThicknesses=thickness,
                                    ids=id, 
                                    predictPositions=seq(0,110,by=5)))
plot(core_2014_cal13)

# add new ages with cal20 in the plot of 2009 core
ages_core2009_cal20 <- BchronCalibrate(ages=core_2009$ages, 
                        ageSds=core_2009$ageSds, 
                        positions=core_2009$position, 
                        calCurves='intcal20')

age1 <- BchronCalibrate(ages = c(975),
                        ageSds = c(15),
                        positions = c(12.5),
                        calCurves = 'intcal20')

age2 <- BchronCalibrate(ages = c(1460),
                        ageSds = c(15),
                        positions = c(42.5),
                        calCurves = 'intcal20')


age3 <- BchronCalibrate(ages = c(2960),
                        ageSds = c(30),
                        positions = c(149.0),
                        calCurves = 'intcal20')

age4 <- BchronCalibrate(ages = c(3760),
                        ageSds = c(30),
                        positions = c(221.2),
                        calCurves = 'intcal20')


library(ggridges)
plot(core_2009_cal13) +
  geom_ridgeline(data = as.data.frame(age1$Date1), 
                 aes(x = ageGrid, 
                     y = age1$Date1$positions,
                     height = densities*1000, # Note the 10000 came from trial and error
                     group = 'New date',
                 ),
                 fill = 'grey',
                 colour = 'black') +
  geom_ridgeline(data = as.data.frame(age2$Date1), 
                 aes(x = ageGrid, 
                     y = age2$Date1$positions,
                     height = densities*1000, # Note the 10000 came from trial and error
                     group = 'New date',
                 ),
                 fill = 'grey',
                 colour = 'black') +
  geom_ridgeline(data = as.data.frame(age3$Date1), 
                 aes(x = ageGrid, 
                     y = age3$Date1$positions,
                     height = densities*1000, # Note the 10000 came from trial and error
                     group = 'New date',
                 ),
                 fill = 'grey',
                 colour = 'black') +
  geom_ridgeline(data = as.data.frame(age4$Date1), 
                 aes(x = ageGrid, 
                     y = age4$Date1$positions,
                     height = densities*1000, # Note the 10000 came from trial and error
                     group = 'New date',
                 ),
                 fill = 'grey',
                 colour = 'black') +
  annotate("text", x = 800, y = 17, label = "Age1 IntCal20") +
  annotate("text", x = 1200, y = 50, label = "Age2 IntCal20") +
  annotate("text", x = 3200, y = 160, label = "Age3 IntCal20") +
  annotate("text", x = 4100, y = 230, label = "Age4 IntCal20") +
  labs(title="Core 2009, pollen (IntCal13)", y="Depth (cm)", x="Age (cal years BP)") +
  theme_classic()

# save plot
ggsave("outputs/core2009_intcal13_intcal20.png",
       plot = last_plot(),
       width=8,
       height=6,
       units="in",
       dpi = 400)


# add new ages with cal20 in the plot of 2014 core
ages_core2009_cal20 <- BchronCalibrate(ages=core_2009$ages, 
                                       ageSds=core_2009$ageSds, 
                                       positions=core_2009$position, 
                                       calCurves='intcal20')

age1 <- BchronCalibrate(ages = c(317),
                        ageSds = c(22),
                        positions = c(55.45),
                        calCurves = 'intcal20')

age2 <- BchronCalibrate(ages = c(332),
                        ageSds = c(44),
                        positions = c(59.50),
                        calCurves = 'intcal20')


age3 <- BchronCalibrate(ages = c(2239),
                        ageSds = c(85),
                        positions = c(90),
                        calCurves = 'intcal20')

age4 <- BchronCalibrate(ages = c(2248),
                        ageSds = c(23),
                        positions = c(108.50),
                        calCurves = 'intcal20')


plot(core_2014_cal13) +
  geom_ridgeline(data = as.data.frame(age1$Date1), 
                 aes(x = ageGrid, 
                     y = age1$Date1$positions,
                     height = densities*1000, # Note the 10000 came from trial and error
                     group = 'New date',
                 ),
                 fill = 'grey',
                 colour = 'black') +
  geom_ridgeline(data = as.data.frame(age2$Date1), 
                 aes(x = ageGrid, 
                     y = age2$Date1$positions,
                     height = densities*1000, # Note the 10000 came from trial and error
                     group = 'New date',
                 ),
                 fill = 'grey',
                 colour = 'black') +
  geom_ridgeline(data = as.data.frame(age3$Date1), 
                 aes(x = ageGrid, 
                     y = age3$Date1$positions,
                     height = densities*1000, # Note the 10000 came from trial and error
                     group = 'New date',
                 ),
                 fill = 'grey',
                 colour = 'black') +
  geom_ridgeline(data = as.data.frame(age4$Date1), 
                 aes(x = ageGrid, 
                     y = age4$Date1$positions,
                     height = densities*1000, # Note the 10000 came from trial and error
                     group = 'New date',
                 ),
                 fill = 'grey',
                 colour = 'black') +
  annotate("text", x = 450, y = 65, label = "Age1 IntCal20") +
  annotate("text", x = 600, y = 50, label = "Age2 IntCal20") +
  annotate("text", x = 2300, y = 92, label = "Age3 IntCal20") +
  annotate("text", x = 2300, y = 113, label = "Age4 IntCal20") +
  labs(title="Core 2014, diatoms (IntCal13)", y="Depth (cm)", x="Age (cal years BP)") +
  theme_classic()

# save plot
ggsave("outputs/core2014_intcal13_intcal20.png",
       plot = last_plot(),
       width=8,
       height=6,
       units="in",
       dpi = 400)
