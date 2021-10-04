# Here model PrC of diatoms, pollen, agropastoralism, XRF with GAM, predict values over diatom time-series 
# and use the resulting dataset to perform multivariate GAM and multivariate ordination. This way, predicted data for
# the different variables will have the same time-series


## read in PrC datasets: diatom, pollen and agropastoralism
diatomPrC <- readRDS("outputs/PrC-cores-diatoms.rds") %>%
  filter(Lake=="llaviucu") %>%
  mutate(elapsedTime = abs(upper_age - lower_age))%>%
  mutate(negAge=-Age)
  
## Pollen PrC
pollenPrC <- read.csv("outputs/pollen-PrC.csv") %>%
  mutate(Age=upper_age) %>%
  mutate(negAge=-Age)%>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 
pollenPrC <- pollenPrC[order(pollenPrC$Age),] #order time
mean(diff(pollenPrC$Age))


## Agropastoralism PrC
agropastPrC <- read.csv("outputs/agropastoralism-PrC.csv") %>%
  mutate(Age=upper_age) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) 
agropastPrC <- agropastPrC[order(agropastPrC$Age),] #order time
mean(diff(agropastPrC$Age))

#Arrange diatom PrC in time order
diatomPrC <- diatomPrC[order(diatomPrC$Age),] #order time
diatomPrC <- na.omit(diatomPrC[, c("Age", "PrC", "elapsedTime")])
diatomPrC$Age <- round(diatomPrC$Age, 1)
mean(diff(diatomPrC$Age))

#check the median combined age interval
median(mean(diff(diatomPrC$Age), mean(diff(pollenPrC$Age)))) #20 years


# Diatom PrC GAM
gamPrC_diat <- gam(PrC ~ s(Age, k=30),
            data = diatomPrC, method = "REML", 
            weights = elapsedTime / mean(elapsedTime),
            select = TRUE, family = gaussian(link="identity"),
            na.action = na.omit) #places notes at the deciles of sample ages

gam.check(gamPrC_diat)
appraise(gamPrC_diat)


## Here the synthetic data to predict should be over the range of diatom ages to get the same ages for the synchronous/asynchronous model
diatom_plot_data <- with(diatomPrC, 
                         as_tibble(expand.grid(Age = seq(min(diatomPrC$Age), max(diatomPrC$Age)))))
# Predict over the range of new values
diatom_mod_fit <- predict(gamPrC_diat, 
                          newdata = diatom_plot_data,
                          se.fit = TRUE)
diatom_plot_data$mod_fit <- as.numeric(diatom_mod_fit$fit)

# For one model only
diatom_plot_data <- gather(diatom_plot_data, key=model, value=fit, mod_fit)
diatom_plot_data <- mutate(diatom_plot_data, se= c(as.numeric(diatom_mod_fit$se.fit)),
                           upper = exp(fit + (2 * se)),
                           lower = exp(fit - (2 * se)),
                           fit   = exp(fit))


# Plot
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

diatom_plt <- ggplot(diatom_plot_data) +
  geom_ribbon(aes(x=Age,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= diatomPrC, aes(x = Age, y = PrC), size=0.06) +
  geom_line(aes(x = Age, y = fit, color = model))+
  labs(y = "PrC", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2") +
  scale_colour_brewer(name = "",
                      palette = "Dark2")+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

diatom_plt

# create a name vector for the dataset
diatom_plot_data$variable <- "diatoms"


# Pollen PrC GAM
gamPrC_pollen <- gam(PrC ~ s(Age, k=30),
                   data = pollenPrC, method = "REML", 
                   weights = elapsedTime / mean(elapsedTime), family = gaussian(link="identity"),
                   select = TRUE, na.action = na.omit) #places notes at the deciles of sample ages

gam.check(gamPrC_pollen)
appraise(gamPrC_pollen)

## Here the synthetic data to predict should be over the range of diatom ages to get the same ages for the synchronous/asynchronous model
pollen_plot_data <- with(pollenPrC, 
                         as_tibble(expand.grid(Age = seq(min(diatomPrC$Age), max(diatomPrC$Age)))))
# Predict over the range of new values
pollen_mod_fit <- predict(gamPrC_pollen, 
                          newdata = pollen_plot_data,
                          se.fit = TRUE)
pollen_plot_data$mod_fit <- as.numeric(pollen_mod_fit$fit)

# For one model only
pollen_plot_data <- gather(pollen_plot_data, key=model, value=fit, mod_fit)
pollen_plot_data <- mutate(pollen_plot_data, se= c(as.numeric(pollen_mod_fit$se.fit)),
                           upper = exp(fit + (2 * se)),
                           lower = exp(fit - (2 * se)),
                           fit   = exp(fit))


# Plot
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

pollen_plt <- ggplot(pollen_plot_data) +
  geom_ribbon(aes(x=Age,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= pollenPrC, aes(x = Age, y = PrC), size=0.06) +
  geom_line(aes(x = Age, y = fit, color = model))+
  labs(y = "PrC", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2") +
  scale_colour_brewer(name = "",
                      palette = "Dark2")+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

pollen_plt

# create a name vector for the dataset
pollen_plot_data$variable <- "pollen" 


# Agropastoralism PrC GAM
gamPrC_agropast <- gam(PrC ~ s(Age, k=30),
                     data = agropastPrC, method = "REML", 
                     weights = elapsedTime / mean(elapsedTime),
                     family=gaussian(link="identity"),
                     select = TRUE) #places notes at the deciles of sample ages

gam.check(gamPrC_agropast)
appraise(gamPrC_agropast)

## Here the synthetic data to predict should be over the range of diatom ages to get the same ages for the synchronous/asynchronous model
agropast_plot_data <- with(agropastPrC, 
                         as_tibble(expand.grid(Age = seq(min(diatomPrC$Age), max(diatomPrC$Age)))))
# Predict over the range of new values
agropast_mod_fit <- predict(gamPrC_agropast, 
                          newdata = agropast_plot_data,
                          se.fit = TRUE)
agropast_plot_data$mod_fit <- as.numeric(agropast_mod_fit$fit)

# For one model only
agropast_plot_data <- gather(agropast_plot_data, key=model, value=fit, mod_fit)
agropast_plot_data <- mutate(agropast_plot_data, se= c(as.numeric(agropast_mod_fit$se.fit)),
                           upper = exp(fit + (2 * se)),
                           lower = exp(fit - (2 * se)),
                           fit   = exp(fit))


# Plot
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

pollen_plt <- ggplot(agropast_plot_data) +
  geom_ribbon(aes(x=Age,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= agropastPrC, aes(x = Age, y = PrC), size=0.06) +
  geom_line(aes(x = Age, y = fit, color = model))+
  labs(y = "PrC", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2") +
  scale_colour_brewer(name = "",
                      palette = "Dark2")+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

pollen_plt

# create a name vector for the dataset
agropast_plot_data$variable <- "agropastoralism" 


# Charcoal
#read Pollen data (drop Lycopodium counts)
pollen <- read.csv("data/llaviucu_pollen_counts.csv") %>%
  select(everything(), -contains("Lycopodium"))

#read in Charcoal data
charcoal <- read.csv("data/llaviucu_pollen.csv") %>%
  select(depth,Charcoal.cc.) %>%
  left_join(pollen, by="depth") %>%
  select(depth, Charcoal.cc., upper_age, lower_age) %>%
  na.omit() 

elapsed <- diff(charcoal$upper_age)
charcoal <- charcoal[-nrow(charcoal),] #remove the last observation
charcoal$elapsedTime <- elapsed

charcoal <- charcoal %>% 
  mutate(Age=upper_age) %>%
  mutate(log_charc=log10(Charcoal.cc.+0.25))

# Model GAM
set.seed(10) #set a seed so this is repeatable
charcoal_gam <- gam(log_charc ~ s(Age, k=30),  weights = elapsedTime / mean(elapsedTime),
              data=charcoal,  select = TRUE, family = gaussian(link = "identity"),
              method = "REML")

gam.check(charcoal_gam)
appraise(charcoal_gam)

## Here the synthetic data to predict should be over the range of diatom ages to get the same ages for the synchronous/asynchronous model
charcoal_plot_data <- with(charcoal, 
                           as_tibble(expand.grid(Age = seq(min(diatomPrC$Age), max(diatomPrC$Age)))))
# Predict over the range of new values
charcoal_mod_fit <- predict(charcoal_gam, 
                            newdata = charcoal_plot_data,
                            se.fit = TRUE)
charcoal_plot_data$mod_fit <- as.numeric(charcoal_mod_fit$fit)

# For one model only
charcoal_plot_data <- gather(charcoal_plot_data, key=model, value=fit, mod_fit)
charcoal_plot_data <- mutate(charcoal_plot_data, se= c(as.numeric(charcoal_mod_fit$se.fit)),
                             upper = exp(fit + (2 * se)),
                             lower = exp(fit - (2 * se)),
                             fit   = exp(fit))


# Plot
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

charcoal_plt <- ggplot(charcoal_plot_data) +
  geom_ribbon(aes(x=Age,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= charcoal, aes(x = Age, y = Charcoal.cc.), size=0.06) +
  geom_line(aes(x = Age, y = fit, color = model))+
  #scale_y_continuous(trans = "log10")+
  labs(y = "Charcoal (cc)", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2") +
  scale_colour_brewer(name = "",
                      palette = "Dark2")+
  theme(legend.position = "top",
        strip.text = element_text(size=10))
charcoal_plt

# create a name vector for the dataset
charcoal_plot_data$variable <- "charcoal"


## Model Ti-GAM
# read in XRF Llaviucu data (Majoi=pollen record)
llaviucu_xrf <- read.csv("data/XRFTransformed3000yrs.csv") %>% rename(age=age_calBP) 
xrf_data <- na.omit(llaviucu_xrf[, c("age", "Ti", "Si", "MnFe", "IncCoh")]) %>% 
  mutate(negAge=-age) %>%
  mutate(Age=age) %>%
  mutate(logTi=log10(Ti+0.25)) %>%
  mutate(logMnFe=log10(MnFe+025))

elapsed <- diff(xrf_data$Age)
xrf_data <- xrf_data[-nrow(xrf_data),] #remove the last observation
xrf_data$elapsedTime <- elapsed

#model GAMs 
set.seed(10) #set a seed so this is repeatable
Ti_gam <- gam(Ti ~ s(Age, k=20),  weights = elapsedTime / mean(elapsedTime),
              data=xrf_data,  select = TRUE, family = nb,
              method = "REML")

Si_gam <- gam(Si ~ s(Age, k=20),  weights = elapsedTime / mean(elapsedTime),
              data=xrf_data,  select = TRUE, family = nb,
              method = "REML")

MnFe_gam <- gam(logMnFe ~ s(Age, k=20),  weights = elapsedTime / mean(elapsedTime),
              data=xrf_data,  select = TRUE, family = gaussian(link = "identity"),
              method = "REML")

appraise(MnFe_gam)
gam.check(MnFe_gam)


## Here the synthetic data to predict should be over the range of diatom ages to get the same ages for the synchronous/asynchronous model
xrf_plot_data <- with(xrf_data, 
                     as_tibble(expand.grid(Age = seq(min(diatomPrC$Age), max(diatomPrC$Age)))))

mod <- Ti_gam
mod <- Si_gam
mod <- MnFe_gam

# Predict over the range of new values
xrf_mod_fit <- predict(mod, 
                      newdata = xrf_plot_data,
                      se.fit = TRUE)
xrf_plot_data$mod_fit <- as.numeric(xrf_mod_fit$fit)


# For one model only
xrf_plot_data <- gather(xrf_plot_data, key=model, value=fit, mod_fit)
xrf_plot_data <- mutate(xrf_plot_data, se= c(as.numeric(xrf_mod_fit$se.fit)),
                       upper = exp(fit + (2 * se)),
                       lower = exp(fit - (2 * se)),
                       fit   = exp(fit))

# Plot
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

xrf_plt <- ggplot(xrf_plot_data) +
  geom_ribbon(aes(x=Age,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= xrf_data, aes(x = Age, y = MnFe), size=0.06) +
  #scale_y_continuous(trans = "log10")+
  geom_line(aes(x = Age, y = fit, color = model))+
  labs(y = "Ti (counts)", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2") +
  scale_colour_brewer(name = "",
                      palette = "Dark2")+
  theme(legend.position = "top",
        strip.text = element_text(size=10))
xrf_plt

# create a name vector for the dataset
#xrf_plot_data$variable <- "Ti"
xrf_plot_data$variable <- "Si"
xrf_plot_data$variable <- "MnFe"


## Combine all GAM-inferred variables
multiproxy <- rbind(diatom_plot_data,pollen_plot_data,agropast_plot_data,
                    charcoal_plot_data,ti_plot_data,xrf_plot_data) %>%
  mutate(variable=factor(variable))

## trim to multiproxy dataset to diatom dataset time visits (response variable with more than time series)
multiproxy_trimmed <- multiproxy[multiproxy$Age %in% diatomPrC$Age,]

#plot multiproxy time series
ggplot(multiproxy_trimmed) + 
  geom_ribbon(aes(x=Age,
                  ymin = lower,
                  ymax = upper,
                  fill = variable),
              alpha=0.2)+
  geom_line(aes(x = Age, y = fit, color = variable))+
  labs(y = "PrC GAM", x = "Age (cal yr BP)") +
  facet_wrap(~variable, scales = "free")+
  scale_fill_brewer(name = "", palette = "Dark2") +
  scale_colour_brewer(name = "",
                      palette = "Dark2")+
  theme(legend.position = "bottom",
        strip.text = element_text(size=10)) +
  theme_bw()

# Save plot
ggsave("outputs/PrCGAMinferredmultiproxy.png",
       plot = last_plot(),
       width = 10,
       height=8,
       units="in",
       dpi = 400)

## Make it wide for multivariate ordination
multiproxy_wide <- multiproxy_trimmed %>% 
  select(-c("model", "se", "upper", "lower")) %>%
  spread(key = variable, value = fit) 

# Make multivariate ordination
rownames(multiproxy_wide) <- multiproxy_wide$Age
elapsed <- diff(multiproxy_wide$Age)
multiproxy_wide <- multiproxy_wide[-nrow(multiproxy_wide),] #remove the last observation
multiproxy_wide$elapsedTime <- elapsed

#save dataset
write.csv(multiproxy_wide, "outputs/multiproxy_gam_inferred.csv")


# Run NMDS
age_nmds <- multiproxy_wide[,names(multiproxy_wide) %in% c("Age")]
multiproxy_wide <- multiproxy_wide[,!names(multiproxy_wide) %in% c("Age")]
multiproxy_wide <- multiproxy_wide[,!names(multiproxy_wide) %in% c("elapsedTime")]

nmds <- metaMDS(multiproxy_wide)
plot(nmds)
plot(nmds, display = "species")
text(nmds, display = "species")

#Extract NMDS site and species scores
scrs_sites <- scores(nmds, display = "sites", choices = 1:2)
scrs_var <- as.data.frame(scores(nmds, display = "species", choices = 1:2))

# Plot
nmds_tbl <- mutate(data.frame(scrs_sites), time=age_nmds$Age)

nmds_plt <- ggplot(nmds_tbl, aes(NMDS1,NMDS2)) + 
  xlab("NMDS1") + ylab("NMDS2") +
  geom_point(aes(colour=time)) +
  geom_text(data=scrs_var, aes(NMDS1, NMDS2, label=rownames(scrs_var)))+
  scale_colour_viridis_c()+
  geom_vline(aes(xintercept = 0), linetype = "solid", colour="grey") +
  geom_hline(aes(yintercept = 0), linetype = "solid", colour="grey") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=9))+
  theme_bw()+
  ggtitle("NMDS with interpolated dataset")
nmds_plt


## Here create phase space plots
ggplot(data=multiproxy_wide, aes(x=agropastoralism, y=diatoms))+
  geom_point(aes(colour=Age))+
  scale_colour_viridis_c()

ggplot(data=multiproxy_wide, aes(x=charcoal, y=diatoms))+
  geom_point(aes(colour=Age))+
  scale_colour_viridis_c()
  
ggplot(data=multiproxy_wide, aes(x=Ti, y=pollen))+
  geom_point(aes(colour=Age))+
  scale_colour_viridis_c()

ggplot(data=multiproxy_wide, aes(x=Ti, y=agropastoralism))+
  geom_point(aes(colour=Age))+
  scale_colour_viridis_c()

ggplot(data=multiproxy_wide, aes(x=MnFe, y=diatoms))+
  geom_point(aes(colour=Age))+
  scale_colour_viridis_c()

ggplot(data=multiproxy_wide, aes(x=MnFe, y=agropastoralism))+
  geom_point(aes(colour=Age))+
  scale_colour_viridis_c()

ggplot(data=multiproxy_wide, aes(x=pollen, y=diatoms))+
  geom_point(aes(colour=Age))+
  scale_colour_viridis_c()


#### Multivariate GAM-inferred values multiproxy
multiproxy_wide <- read.csv("outputs/multiproxy_gam_inferred.csv")

mod1 <- gam(diatoms ~ s(Age, k=20, bs="ad") + s(Ti, k=10) 
            + s(pollen, k=10) + s(agropastoralism, k=10),
            data = multiproxy_wide, method = "REML", 
            weights = elapsedTime / mean(elapsedTime),
            select = TRUE, family = gaussian(link="identity"),
            na.action = na.omit,
            knots = list(negAge=quantile(multiproxy_wide$Age, seq(0,1, length=10)))) #places notes at the deciles of sample ages

plot(mod1, page=1, scale = 0)
gam.check(mod1)
summary(mod1)


## Predict
predGam <- cbind(multiproxy_wide, 
                 data.frame(predict.gam(mod1, multiproxy_wide, 
                                        type = "terms" , se.fit = TRUE)))

#plot
var <- predGam$fit.s.agropastoralism.
se.var <- predGam$se.fit.s.agropastoralism.

predGamPlt <- ggplot(predGam, aes(x = Age, y = var)) +
  geom_line() +
  geom_ribbon(aes(ymin = var + (2 * se.var), ymax = var - (2 * se.var), alpha = 0.1)) +
  #scale_x_reverse() +
  labs(y = "Diatom GAM PrC", x = "Cal yr BP", title = "")+
  theme(legend.position = "none")+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()+ theme(legend.position = "none")+
  ggtitle(expression("Agropastoralism" %->% "Diatoms"))

predGamPlt

ggsave("outputs/GAMinferredPrC_covariates.png",
       plot = predGamPlt ,
       width = 10,
       height=8,
       units="in",
       dpi = 400)

