################################################
#######   HGAM-based derivative models  ########
################################################

#######  diatom assemblages

#read dataframe with diatom absolute counts
mergedCores <- read.csv("data/mergedCores_counts4.csv")[-1] #read dataframe with diatom absolute counts including Fondococha

#import dataframe wiht old and new names to group
changes <- read.csv("data/old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
#new1: ecological groups
#new2: harmonized taxonomic names

#filter cores
select <- c("llaviucu")

# this is to select the most common species for a single lake   
core <- mergedCores %>% 
  filter(str_detect(lake, select))

agedepth <- core[, names(core) %in% c("depth", "upper_age", "lower_age", "lake")]
diat <- core[, !names(core) %in% c("depth", "upper_age", "lower_age", "lake")]
diat[is.na(diat)] <- 0

#Select most common species 
criteria <- 0.5 #% of the total samples

n.occur <- apply(diat>0, 2, sum)
diat_red <- diat[, n.occur > (dim(diat)[1])*criteria] #
diat <- cbind(agedepth, diat_red)

##This creates the dataset containing the most common species for a single lake core
#Gather
spp_thin <- diat %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age,  -lake) #don't gather depths, ages and lake variables

diat_data <- spp_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth) %>%
  mutate(total_sample = sum(count)) %>% 
  filter(!total_sample == "0") %>% #this is to remove empty samples
  filter(!upper_age == 0) %>% #this is to remove ages == 0 (triumfo and fondodocha record)
  mutate(log_total_counts = log10(total_sample+1)) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  mutate(negAge = -upper_age) %>%
  mutate(AgeCE = upper_age*(-1)+1950) %>%
  #filter(AgeCE >= "0") %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) %>%
  ungroup() %>%
  group_by(lake) %>%
  mutate(spp = factor(taxa)) 

levels(diat_data$spp)

#model S HGAM : similar smootheness between groups (spp) without global smooth 
#Eric's CHANGE: using negative binomial instead of Poisson, as there seems to be evidence of overdispersion in counts, and it results in smoother trends for each species
set.seed(10) #set a seed so this is repeatable

diatom_gam_S <- gam(count ~ s(negAge, spp, k=20, bs="fs") + offset(log_total_counts),
                    weights = elapsedTime / mean(elapsedTime),
                    data=diat_data, family = nb, 
                    knots = list(negAge=quantile(diat_data$negAge, seq(0,1, length=10))), #places notes at the deciles of sample ages
                    method = "REML")

gam.check(diatom_gam_S)
draw(diatom_gam_S)

#model I HGAM: different smootheness for each taxa without global smooth
diatom_gam_I<- gam(count ~ s(negAge, by=spp, k=20, bs="fs") +
                     s(spp, bs="re") + offset(log_total_counts),
                   weights = elapsedTime / mean(elapsedTime),
                   data=diat_data, family = nb,
                   knots = list(negAge=quantile(diat_data$negAge, seq(0,1, length=10))), #places notes at the deciles of sample ages
                   method = "REML")

gam.check(diatom_gam_I)
draw(diatom_gam_I)

#Compare different model fits using AIC
AIC_table <- AIC(diatom_gam_S, diatom_gam_I)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

#Create synthetic data to predict over a range of ages
diat_plot_data <- with(diat_data, as_tibble(expand.grid(negAge = seq(min(diat_data$negAge), max(diat_data$negAge)),
                                                        spp = factor(levels(diat_data$spp)),
                                                        log_total_counts = mean(log_total_counts))))

diat_modS_fit <- predict(diatom_gam_S, 
                         newdata = diat_plot_data,
                         se.fit = TRUE)

diat_modI_fit <- predict(diatom_gam_I,
                         newdata = diat_plot_data,
                         se.fit = TRUE)

#non-shared trends
diat_plot_data$modS_fit <- as.numeric(diat_modS_fit$fit)
diat_plot_data$modI_fit <- as.numeric(diat_modI_fit$fit)

# comparing non-shared trends
diat_plot_data <- gather(diat_plot_data, key=model, value=fit, modS_fit, modI_fit)

diat_plot_data <- mutate(diat_plot_data, se= c(as.numeric(diat_modS_fit$se.fit),
                                               as.numeric(diat_modI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))

# For only one model
diat_plot_data <- gather(diat_plot_data, key=model, value=fit, modI_fit)
diat_plot_data <- mutate(diat_plot_data, se= c(as.numeric(diat_modI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))

#Plot the model output for non-shared trends, with means plus standard deviations for each model.
diat_plot_model_labels <- paste("Model", c("S", "I"))
diat_plot_model_labels <- factor(diat_plot_model_labels, levels = diat_plot_model_labels)

#non-shared trends
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

diat_plot <- ggplot(diat_plot_data) +
  facet_wrap(~spp, nrow = 4,scales = "free_y")+
  geom_ribbon(aes(x=negAge,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= diat_data, aes(x = negAge, y = count), size=0.06) +
  geom_line(aes(x = negAge, y = fit, color = model))+
  labs(y = "Absolute counts", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2",
                    labels = diat_plot_model_labels) +
  scale_colour_brewer(name = "",
                      palette = "Dark2", labels = diat_plot_model_labels)+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

diat_plot

##Derivatives and posterior distribution simulation
# Eric's code
set.seed(10) #set a seed so this is repeatable
n_sims = 250

n_length = 140

years <- seq(min(diat_plot_data$negAge),
             max(diat_plot_data$negAge),
             length.out = n_length)

#model <- diatom_gam_S
model <- diatom_gam_I

#pred <- diat_modS_fit
pred <- diat_modI_fit

# Generate multivariate normal simulations of model coefficients
random_coefs <- t(rmvn(n_sims, mu = coef(model),V = vcov(model)))

confint_sims <- crossing(spp=unique(diat_plot_data$spp),
                         negAge = seq(min(diat_plot_data$negAge),
                                      max(diat_plot_data$negAge),
                                      length.out = n_length),
                         log_total_counts=0)

map_pred_sims <- predict(model,
                         confint_sims,
                         type = "lpmatrix") %*% random_coefs %>%
  as_data_frame() %>%
  bind_cols(confint_sims)%>%
  gather(key = simulation, value = pred, -negAge, -log_total_counts,-spp)


#specifying the step size for numerical derivative calculations
delta = 0.01

#calculating the predicted value for the current year plus delta
step_ahead_fits = confint_sims %>%
  mutate(negAge = negAge+delta)%>%
  predict(model, 
          ., type = "lpmatrix") %*% random_coefs 


#calculating the predicted value for the current year minus delta
step_behind_fits = confint_sims %>%
  mutate(negAge = negAge-delta)%>%
  predict(model,
          ., type = "lpmatrix") %*% random_coefs 


#using the predicted values for year plus and minus delta to calculate
#derivatives for each species for each simulation
derivs <- calc_1st_deriv(step_behind_fits,step_ahead_fits,delta = delta)%>%
  as_data_frame()%>%
  bind_cols(confint_sims)%>%
  gather(key = simulation,value = deriv, -spp,-negAge, -log_total_counts)

#Creating summaries of derivatives for each simulation for each year
deriv_summaries <- derivs %>%
  group_by(negAge,simulation)%>%
  summarize(deriv_mean = mean(deriv),
            deriv_sd = sd(deriv))%>%
  group_by(negAge)%>% #turning derivative summaries into 95% confidence intervals
  select(-simulation)%>%
  summarize_all(.funs = list(lower = ~quantile(.,probs = 0.025),
                             upper = ~quantile(.,probs = 0.975),
                             med   = ~quantile(.,probs = 0.5)))

## save derivative summaries for later use
write.csv(deriv_summaries, "outputs/diatom-derivatives_24timeSteps.csv", row.names = FALSE)

########################################
####### HGAM pollen assemblages ########
########################################

# read data (drop Lycopodium counts)
pollen <- read.csv("data/llaviucu_pollen_counts.csv") %>%
  select(everything(), -contains("Lycopodium"))

# Read pollen ratios (Majoi's data) and merge charcoal data
pollen <- read.csv("data/llaviucu_pollen.csv") %>%
  select(depth,Charcoal.cc.) %>%
  left_join(pollen, by="depth")
  
#select human disturbance pollen taxa
agropastolarism_indicators <- c("Zea", "Hedyosmum", "Phaseolus", "Ipomoea", "Rumex", 
                                "Alnus", "Cyperaceae", "Cecropia", "Asteracea", "sporormiella","Charcoal.cc.")
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


#Select most common species (only for pollen dataset)
criteria <- 0.6 #% of the total samples

#this should be counts
n.occur <- apply(llaviucu_pollen>0, 2, sum)
pollen_red <- llaviucu_pollen[, n.occur > (dim(llaviucu_pollen)[1])*criteria] #

#assign to subset
pollen <- cbind(agedepth, pollen_red)
#pollen <- cbind(agedepth, agropastolarism)

#Gather
spp_thin <- pollen %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age) #don't gather depths, ages and lake variables

pollen_data <- spp_thin %>%
  group_by(depth) %>%
  mutate(total_sample = sum(count)) %>% 
  filter(!total_sample == "0") %>% #this is to remove empty samples
  filter(!upper_age == 0) %>% #this is to remove ages == 0 (triumfo and fondodocha record)
  mutate(log_total_counts = log10(total_sample+1)) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  mutate(negAge = -upper_age) %>%
  mutate(AgeCE = upper_age*(-1)+1950) %>%
  mutate(elapsedTime = abs(upper_age - lower_age)) %>%
  ungroup() %>%
  mutate(spp = factor(taxa)) 

levels(pollen_data$spp)

#model S HGAM : similar smootheness between groups (spp) without global smooth 
set.seed(10) #set a seed so this is repeatable

pollen_gam_S <- gam(count ~ s(negAge, spp, k=30, bs="fs") + offset(log_total_counts),
                    weights = elapsedTime / mean(elapsedTime),
                    data=pollen_data, family = nb, 
                    method = "REML")

gam.check(pollen_gam_S)
draw(pollen_gam_S)

#model I HGAM: different smootheness for each taxa without global smooth
pollen_gam_I<- gam(count ~ s(negAge, by=spp, k=30, bs="fs") +
                     s(spp, bs="re") + offset(log_total_counts),
                   weights = elapsedTime / mean(elapsedTime),
                   data=pollen_data, family = nb,
                   method = "REML")

gam.check(pollen_gam_I)
draw(pollen_gam_I)

#Compare different model fits using AIC
AIC_table <- AIC(pollen_gam_S, pollen_gam_I)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))
AIC_table

#Create synthetic data to predict over a range of ages

## Here the synthetic data to predict should be over the range of diatom ages to get the same ages for the synchronous/asynchronous model
## see negAge argument
pollen_plot_data <- with(pollen_data, as_tibble(expand.grid(negAge = seq(min(diat_data$negAge), max(diat_data$negAge)),
                                                            spp = factor(levels(pollen_data$spp)),
                                                            log_total_counts = mean(log_total_counts))))

pollen_modS_fit <- predict(pollen_gam_S, 
                           newdata = pollen_plot_data,
                           se.fit = TRUE)

pollen_modI_fit <- predict(pollen_gam_I,
                           newdata = pollen_plot_data,
                           se.fit = TRUE)

#non-shared trends
pollen_plot_data$modS_fit <- as.numeric(pollen_modS_fit$fit)
pollen_plot_data$modI_fit <- as.numeric(pollen_modI_fit$fit)

# comparing non-shared trends
pollen_plot_data <- gather(pollen_plot_data, key=model, value=fit, modS_fit, modI_fit)
pollen_plot_data <- mutate(pollen_plot_data, se= c(as.numeric(pollen_modS_fit$se.fit),
                                                   as.numeric(pollen_modI_fit$se.fit)),
                           upper = exp(fit + (2 * se)),
                           lower = exp(fit - (2 * se)),
                           fit   = exp(fit))

# For one model only
pollen_plot_data <- gather(pollen_plot_data, key=model, value=fit, modI_fit)
pollen_plot_data <- mutate(pollen_plot_data, se= c(as.numeric(pollen_modI_fit$se.fit)),
                           upper = exp(fit + (2 * se)),
                           lower = exp(fit - (2 * se)),
                           fit   = exp(fit))

#Plot the model output for non-shared trends, with means plus standard deviations for each model.
pollen_plot_model_labels <- paste("Model", c("S", "I"))
pollen_plot_model_labels <- factor(pollen_plot_model_labels, levels = pollen_plot_model_labels)

#non-shared trends
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

pollen_plt <- ggplot(pollen_plot_data) +
  facet_wrap(~spp, nrow = 4,scales = "free_y")+
  geom_ribbon(aes(x=negAge,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= pollen_data, aes(x = negAge, y = count), size=0.06) +
  geom_line(aes(x = negAge, y = fit, color = model))+
  labs(y = "Absolute counts", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2",
                    labels = pollen_plot_model_labels) +
  scale_colour_brewer(name = "",
                      palette = "Dark2", labels = pollen_plot_model_labels)+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

pollen_plt

ggsave("outputs/llaviucu_HGAM_agropastoralism.png",
       plot = pollen_plt ,
       width = 10,
       height=8,
       units="in",
       dpi = 400)


## save model results for later use
write.csv(pollen_plot_data, "outputs/pollen-natural-HGAMs-fitted-values.csv", row.names = FALSE)
write.csv(pollen_plot_data, "outputs/pollen-agrospastolarism-HGAMs-fitted-values.csv", row.names = FALSE)


##Derivatives and posterior distribution simulation
## Eric's code
#Function for calculating first derivatives of time series given a before,
#after, and delta step size
calc_1st_deriv = function(fit_before, fit_after,delta) {
  (fit_after-fit_before)/(2*delta)
}

set.seed(10) #set a seed so this is repeatable
n_sims = 250

n_length = 140

years <- seq(min(pollen_plot_data$negAge),
             max(pollen_plot_data$negAge),
             length.out = n_length)

#model <- pollen_gam_S
model <- pollen_gam_I

#pred <- pollen_modS_fit
pred <- pollen_modI_fit

# Generate multivariate normal simulations of model coefficients
random_coefs <- t(rmvn(n_sims, mu = coef(model),V = vcov(model)))

confint_sims <- crossing(spp=unique(pollen_plot_data$spp),
                         negAge = seq(min(pollen_plot_data$negAge),
                                      max(pollen_plot_data$negAge),
                                      length.out = n_length),
                         log_total_counts=0)

map_pred_sims <- predict(model,
                         confint_sims,
                         type = "lpmatrix") %*% random_coefs %>%
  as_data_frame() %>%
  bind_cols(confint_sims)%>%
  gather(key = simulation, value = pred, -negAge, -log_total_counts,-spp)


#specifying the step size for numerical derivative calculations
delta = 0.01

#calculating the predicted value for the current year plus delta
step_ahead_fits = confint_sims %>%
  mutate(negAge = negAge+delta)%>%
  predict(model, 
          ., type = "lpmatrix") %*% random_coefs 


#calculating the predicted value for the current year minus delta
step_behind_fits = confint_sims %>%
  mutate(negAge = negAge-delta)%>%
  predict(model,
          ., type = "lpmatrix") %*% random_coefs 


#using the predicted values for year plus and minus delta to calculate
#derivatives for each species for each simulation
derivs <- calc_1st_deriv(step_behind_fits,step_ahead_fits,delta = delta)%>%
  as_data_frame()%>%
  bind_cols(confint_sims)%>%
  gather(key = simulation,value = deriv, -spp,-negAge, -log_total_counts)

#Creating summaries of derivatives for each simulation for each year
deriv_summaries <- derivs %>%
  group_by(negAge,simulation)%>%
  summarize(deriv_mean = mean(deriv),
            deriv_sd = sd(deriv))%>%
  group_by(negAge)%>% #turning derivative summaries into 95% confidence intervals
  select(-simulation)%>%
  summarize_all(.funs = list(lower = ~quantile(.,probs = 0.025),
                             upper = ~quantile(.,probs = 0.975),
                             med   = ~quantile(.,probs = 0.5)))

#Plotting mean rate of change plus the 95% CI
mean_plot <- deriv_summaries %>%
  ggplot(aes(x = negAge, 
             y = deriv_mean_med, 
             ymin = deriv_mean_lower,
             ymax = deriv_mean_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("Cal years BP") +
  theme_bw()
mean_plot 

#Plotting standard deviation of rate of change plus the 95% CI
sd_plot <- deriv_summaries %>%
  ggplot(aes(x = negAge, 
             y = deriv_sd_med, 
             ymin=deriv_sd_lower,
             ymax=deriv_sd_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("Cal years BP") +
  theme_bw()
sd_plot


## save derivative summaries for later use
write.csv(deriv_summaries, "outputs/pollen-natural-derivatives.csv", row.names = FALSE)
write.csv(deriv_summaries, "outputs/pollen-agropastoralism-derivatives.csv", row.names = FALSE)

# save derivative summaries at 24 time steps
write.csv(deriv_summaries, "outputs/pollen-agropastoralism-derivatives_24timeSteps.csv", row.names = FALSE)
write.csv(deriv_summaries, "outputs/pollen-natural-derivatives_24timeSteps.csv", row.names = FALSE)

#######################
## Ti GAM analysis ####
#######################

# read in XRF Llaviucu data (Majoi=pollen record)
llaviucu_xrf <- read.csv("data/XRFTransformed3000yrs.csv") %>% rename(age=age_calBP) 
Ti_xrf <- na.omit(llaviucu_xrf[, c("age", "Ti")]) %>% 
            mutate(negAge=-age) %>%
            mutate(logTi=log10(Ti+0.25))

elapsed <- diff(Ti_xrf$age)
Ti_xrf <- Ti_xrf[-nrow(Ti_xrf),] #remove the last observation
Ti_xrf$elapsedTime <- elapsed


#model GAM: 
set.seed(10) #set a seed so this is repeatable

Ti_gam <- gam(Ti ~ s(negAge, k=20),  weights = elapsedTime / mean(elapsedTime),
                    data=Ti_xrf,  select = TRUE, family = nb, 
                    method = "REML")

pacf(residuals(Ti_gam)) # indicates AR1
plot(Ti_gam, scale = 0)
gam.check(Ti_gam)
summary(Ti_gam)


## Here the synthetic data to predict should be over the range of diatom ages to get the same ages for the synchronous/asynchronous model
## see negAge argument
ti_plot_data <- with(Ti_xrf, 
                     as_tibble(expand.grid(negAge = seq(min(diat_data$negAge), max(diat_data$negAge)))))
                                                            

# Predict over the range of new values
Ti_mod_fit <- predict(Ti_gam, 
                           newdata = ti_plot_data,
                           se.fit = TRUE)


#
ti_plot_data$mod_fit <- as.numeric(Ti_mod_fit$fit)


# For one model only
ti_plot_data <- gather(ti_plot_data, key=model, value=fit, mod_fit)
ti_plot_data <- mutate(ti_plot_data, se= c(as.numeric(Ti_mod_fit$se.fit)),
                           upper = exp(fit + (2 * se)),
                           lower = exp(fit - (2 * se)),
                           fit   = exp(fit))

# Plot
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

Ti_plt <- ggplot(ti_plot_data) +
  geom_ribbon(aes(x=negAge,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= Ti_xrf, aes(x = negAge, y = Ti), size=0.06) +
  geom_line(aes(x = negAge, y = fit, color = model))+
  labs(y = "Ti (counts)", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2") +
  scale_colour_brewer(name = "",
                      palette = "Dark2")+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

Ti_plt

# Save plot
ggsave("outputs/Ti_gam_model.png",
       plot = Ti_plt ,
       width = 10,
       height=8,
       units="in",
       dpi = 400)



##Derivatives and posterior distribution simulation
#Function for calculating first derivatives of time series given a before,
#after, and delta step size
calc_1st_deriv = function(fit_before, fit_after,delta) {
  (fit_after-fit_before)/(2*delta)
}

set.seed(10) #set a seed so this is repeatable
n_sims = 250

n_length = 140

years <- seq(min(ti_plot_data$negAge),
             max(ti_plot_data$negAge),
             length.out = n_length)

model <- Ti_gam
pred <- Ti_mod_fit

# Generate multivariate normal simulations of model coefficients
random_coefs <- t(rmvn(n_sims, mu = coef(model),V = vcov(model)))

confint_sims <- crossing(negAge = seq(min(ti_plot_data$negAge),
                                      max(ti_plot_data$negAge),
                                      length.out = n_length))

map_pred_sims <- predict(model,
                         confint_sims,
                         type = "lpmatrix") %*% random_coefs %>%
  as_data_frame() %>%
  bind_cols(confint_sims)%>%
  gather(key = simulation, value = pred, -negAge)


#specifying the step size for numerical derivative calculations
delta = 0.01

#calculating the predicted value for the current year plus delta
step_ahead_fits = confint_sims %>%
  mutate(negAge = negAge+delta)%>%
  predict(model, 
          ., type = "lpmatrix") %*% random_coefs 


#calculating the predicted value for the current year minus delta
step_behind_fits = confint_sims %>%
  mutate(negAge = negAge-delta)%>%
  predict(model,
          ., type = "lpmatrix") %*% random_coefs 


#using the predicted values for year plus and minus delta to calculate
#derivatives for each species for each simulation
derivs <- calc_1st_deriv(step_behind_fits,step_ahead_fits,delta = delta)%>%
  as_data_frame()%>%
  bind_cols(confint_sims)%>%
  gather(key = simulation,value = deriv, -negAge)

#Creating summaries of derivatives for each simulation for each year
deriv_summaries <- derivs %>%
  group_by(negAge)%>%
  summarize(deriv_mean = mean(deriv),
            deriv_sd = sd(deriv))%>%
  group_by(negAge)%>% #turning derivative summaries into 95% confidence intervals
  #select(-simulation)%>%
  summarize_all(.funs = list(lower = ~quantile(.,probs = 0.025),
                             upper = ~quantile(.,probs = 0.975),
                             med   = ~quantile(.,probs = 0.5)))

#Plotting mean rate of change plus the 95% CI
mean_plot <- deriv_summaries %>%
  ggplot(aes(x = negAge, 
             y = deriv_mean_med, 
             ymin = deriv_mean_lower,
             ymax = deriv_mean_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("Cal years BP") +
  theme_bw()
mean_plot 

#Plotting standard deviation of rate of change plus the 95% CI
sd_plot <- deriv_summaries %>%
  ggplot(aes(x = negAge, 
             y = deriv_sd_med, 
             ymin=deriv_sd_lower,
             ymax=deriv_sd_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("Cal years BP") +
  theme_bw()
sd_plot

# Save derivative results
write.csv(deriv_summaries, "outputs/Ti-gam-derivatives_24timeSteps.csv", row.names = FALSE)


## Read HGAM-derivatives results for plotting
deriv_summaries_pollen <- read.csv("outputs/pollen-natural-derivatives.csv") %>%
  mutate(proxy="pollen") 

deriv_summaries_All_pollen <- read.csv("outputs/pollen-agropastoralism-derivatives.csv") %>%
  mutate(proxy="agropastoralism") %>%
  bind_rows(deriv_summaries_pollen)

#Plotting mean rate of change plus the 95% CI
mean_plot <- deriv_summaries_All_pollen %>%
  ggplot(aes(x = negAge, 
             y = deriv_mean_med, 
             ymin = deriv_mean_lower,
             ymax = deriv_mean_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  facet_wrap(~proxy, scales = "free")+
  xlab("Cal years BP") +
  theme_bw()
mean_plot 


#Plotting standard deviation of rate of change plus the 95% CI
sd_plot <- deriv_summaries_All_pollen %>%
  ggplot(aes(x = negAge, 
             y = deriv_sd_med, 
             ymin=deriv_sd_lower,
             ymax=deriv_sd_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  facet_wrap(~proxy, scales = "free")+
  xlab("Cal years BP") +
  theme_bw()
sd_plot

mean_sd_pl_pollen <- plot_grid(mean_plot, sd_plot, align = "v", nrow = 2, labels="auto")
mean_sd_pl_pollen

ggsave("outputs/llaviucu_HGAM_pollen.png",
       plot = mean_sd_pl_pollen ,
       width = 10,
       height=8,
       units="in",
       dpi = 400)



## Read in Diatom HGAM derivatives and merge with pollen derivatives
deriv_summaries_all <- read.csv("outputs/diatom-derivatives.csv") %>%
  mutate(proxy="diatoms") %>%
  bind_rows(deriv_summaries_All_pollen)

write.csv(deriv_summaries_all, "outputs/diatom_pollen_derivatives.csv")

#Plotting mean rate of change plus the 95% CI
mean_plot <- deriv_summaries_all %>%
  ggplot(aes(x = negAge, 
             y = deriv_mean_med, 
             ymin = deriv_mean_lower,
             ymax = deriv_mean_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("") +
  facet_wrap(~proxy, scales="free")+
  ggtitle("")+
  theme_bw()
mean_plot 


#Plotting standard deviation of rate of change plus the 95% CI
sd_plot <- deriv_summaries_all %>%
  ggplot(aes(x = negAge, 
             y = deriv_sd_med, 
             ymin=deriv_sd_lower,
             ymax=deriv_sd_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("Cal years BP") +
  facet_wrap(~proxy, scales="free")+
  theme_bw()
sd_plot


# natural, agropastoralism and diatoms
library(cowplot)
mean_sd_pl_all <- plot_grid(mean_plot, sd_plot, align = "hv", 
                            nrow=2,ncol = 1, axis="tblr", labels="auto")
mean_sd_pl_all

# save plot
ggsave("outputs/llaviucu_HGAMderivative_GI_pollen_diatom_all_sync.png",
       plot = mean_sd_pl_all,
       width = 10,
       height=8,
       units="in",
       dpi = 400)




