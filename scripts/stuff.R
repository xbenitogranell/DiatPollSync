
#Here I'm trying to create a new dataset to predict pollen over the same range of diatom ages
# see negAge vector

#Create synthetic data to predict over a range of ages
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


##Derivatives and posterior distribution simulation
## Eric's code
#Function for calculating first derivatives of time series given a before,
#after, and delta step size
calc_1st_deriv = function(fit_before, fit_after,delta) {
  (fit_after-fit_before)/(2*delta)
}

set.seed(10) #set a seed so this is repeatable
n_sims = 250

years <- seq(min(pollen_plot_data$negAge),
             max(pollen_plot_data$negAge),
             length.out = 50)

#model <- pollen_gam_S
model <- pollen_gam_I

#pred <- pollen_modS_fit
pred <- pollen_modI_fit

# Generate multivariate normal simulations of model coefficients
random_coefs <- t(rmvn(n_sims, mu = coef(model),V = vcov(model)))

confint_sims <- crossing(spp=unique(pollen_plot_data$spp),
                         negAge = seq(min(pollen_plot_data$negAge),
                                      max(pollen_plot_data$negAge),
                                      length.out = 50),
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




