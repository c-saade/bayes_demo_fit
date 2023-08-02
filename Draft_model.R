## Draft model ##
# A standalone R-script to fit and diagnose a demographic model
# (competitive Lotka-Volterra) to experimental bi-specific data.

# Cleaning environment
rm(list=ls())

# loading packages ############################################################

# data wrangling
library(dplyr)

# figures
library(ggplot2)
theme_set(theme_light())
library(ggpubr)

# models
library(rstan)
library(deSolve)
library(coda)

dir.create("out", showWarnings = FALSE) # creating a directory to store the outputs

# Loading and preparing data ######################################################

density_data = read.table("./Data/density_data.csv", header = T, sep = ',')

# filtering the data
filtered_data <- 
  density_data %>%
  filter(community == 'Ble_Tet') %>% # keeping only the bi-specific time series (Blepharisma sp. + Tetrahymena sp.)
  select(c(measure_point, hours, Ble, Tet)) %>% # keeping only time variables and species densities
  group_by(measure_point) %>% # averaging each measurement over replicates
  summarise(t = mean(hours), n1 = mean(Ble), n2 = mean(Tet)) %>%
  arrange(t)

# taking a look at the time series
panel_a = ggplot(data = filtered_data) +
  geom_point(mapping = aes(x = t, y = n1), color = 'blue') +
  xlab('Time (h)') +
  ylab('Density of Blepharisma sp. (indiv./mL)')

panel_b = ggplot(data = filtered_data) +
  geom_point(mapping = aes(x = t, y = n2), color = 'red') +
  xlab('Time (h)') +
  ylab('Density of Tetrahymena sp. (indiv./mL)')

ggsave(('out/time_series.pdf'), ggarrange(panel_a, panel_b, nrow = 2), width = 5, height = 6)

# turning data into a list to give to the Rstan sampler.
data = list(n  = nrow(filtered_data), # number of observations
            t = filtered_data$t, # vector of times
            n1 = filtered_data$n1, # vector of density of species 1 (Ble)
            n2 = filtered_data$n2) # vector of density of species 2 (Tet)

## Declaring and compiling Stan model ###################################################

model_str = '
functions{
  real[] odemodel(real t, real[] N, real[] p, real[] x_r, int[] x_i){
    // p[1]=r1, p[2] = r2, p[3] = K1, p[4] = K2, p[5] = alpha1, p[6] = alpha2 
    real dNdt[2]; 
    dNdt[1] = p[1]*N[1]*(1 - (N[1] + p[5]*N[2])/p[3]);
    dNdt[2] = p[2]*N[2]*(1 - (p[6]*N[1] + N[2])/p[4]);
    return dNdt;
  }
}

data{
  int n; // number of observations
  real t[n]; // time
  real n1[n]; // observations n1
  real n2[n]; // observations n2
}

parameters{
  real<lower=0> r1; // growth rate
  real<lower=0> r2; // growth rate
  real<lower=0> K1; // carrying capacity
  real<lower=0> K2; // carrying capacity
  real<lower=0> alpha1; // comp term
  real<lower=0> alpha2; // comp term
  real<lower=0> n10sim; // initial density n1
  real<lower=0> n20sim; // initial density n2
  real<lower=0> sdev1;
  real<lower=0> sdev2;
}

model{
  real p[6]; // vector of parameters for the ODE
  real simval[n-1,2]; // simulated values, matrix. dim1 = time without t0, dim2 = dim_ODE = 2 (S = 1, I = 2)
  
  // priors 
  r1 ~ lognormal(-3,1);
  r2 ~ lognormal(-2,1);
  K1 ~ lognormal(5, 0.8);
  K2 ~ lognormal(8, 0.8);
  alpha1 ~ lognormal(-3,1);
  alpha2 ~ lognormal(2.5,1);
  n10sim ~ normal(n1[1],5);
  n20sim ~ normal(n2[1],10);
  sdev1 ~ gamma(1, 1);
  sdev2 ~ gamma(1, 1);
  
  // parameters for integrator
  p[1] = r1;
  p[2] = r2;
  p[3] = K1;
  p[4] = K2;
  p[5] = alpha1;
  p[6] = alpha2;

  // integrate ODE
  simval = integrate_ode_rk45(odemodel, {n10sim, n20sim}, t[1], t[2:n], p, rep_array(0.0,0), rep_array(0,0));
  // likelihood
  n1[1] ~ normal(n10sim, sdev1);
  n2[1] ~ normal(n20sim, sdev2);
  for (i in 2:n){
    n1[i] ~ normal(simval[i-1, 1], sdev1);
    n2[i] ~ normal(simval[i-1, 2], sdev2);
  }
}

generated quantities{
  real r_ratio = r1/r2;
  real alpha_ratio = alpha1/alpha2;
}
'

# compiling the model
model = stan_model(model_code=model_str)


# Declaring model options and initial conditions ##################################

# stan options
chains = 4 # number of parallel chains
options(mc.cores = chains) # number of core used (check that you have at least 3)

# number of total iterations and warm-up steps
iter   =  4000
warmup =  2000

# initial values for sampling 
init=rep(list(list(r1=0.01,
                   r2=0.1,
                   K1 = 150,
                   K2 = 4000,
                   alpha1 = 1,
                   alpha2 = 1,
                   n10sim=5,
                   n20sim=500,
                   sdev1 = 1,
                   sdev2 = 1
))
,chains)


# Fitting model ####################################################################

# Note that the fit can take some time, depending on your hardware.
# You can skip the fitting by commenting out the next two commands
# ('fit = sampling...' and 'save(fit...)') and uncommenting
# load("./out/fit_posterior.RData")

fit = sampling(model,
                  data=data,
                  iter=iter,
                  warmup=warmup,
                  chains=chains,
                  init=init,
                  control = list(adapt_delta = 0.98, max_treedepth=12),
                  refresh=100,
                  seed = 123)

save(fit, file="./out/fit_posterior.RData")

#load("./out/fit_posterior.RData")

# Model diagnostics #############################################################
params = c("r1","r2", "K1", "K2", "alpha1", "alpha2", "r_ratio", "alpha_ratio")

# Checking posterior properties (rhat should be as close to 1 as possible)
print(fit)

# saving a plot of the chains
samples=As.mcmc.list(fit)
pdf('out/chains.pdf')
plot(samples[, params])
dev.off()

# saving pair plot of parameters
pdf('out/pair.pdf')
pairs(fit, pars=params)
dev.off()

## Posterior predictions  ####################################################

# declaring the lotka-volterra model
ode.model = function(t,N,p){
  r1 = p$r1
  r2 = p$r2
  K1 = p$K1
  K2 = p$K2
  alpha1 = p$alpha1
  alpha2 = p$alpha2
  dn1 = r1*N[1]*(1-(N[1]+alpha1*N[2])/K1)
  dn2 = r2*N[2]*(1-(alpha2*N[1]+N[2])/K2)
  return(list(c(dn1, dn2)))
}


# sampling 1000 values from the posterior distribution and
# computing the predicted dynamics for each
n_post = 1000
times = seq(min(data$t), max(data$t), length.out = 200)
posteriors = as.matrix(fit)
for (k in 1:n_post){
  par = posteriors[sample(1:nrow(posteriors), 1),]
  sim = ode(c(par['n10sim'], par['n20sim']),
            times, ode.model, list(r1 = par['r1'],
                                   r2= par['r2'],
                                   K1 = par['K1'],
                                   K2 = par['K2'],
                                   alpha1 = par['alpha1'],
                                   alpha2 = par['alpha2']),
            method = 'ode45')
  
  temp  = data.frame(time = sim[,1], n1 = sim[,2], n2 = sim[,3], id = k)
  
  if (k == 1) {
    predictions = temp
  } else {
    predictions = rbind(predictions, temp)
  }
  
}

# plotting the 1000 predictions over the data
panel_a = ggplot(filtered_data) +
  geom_point(mapping = aes(x = t, y = n1), color = 'blue') +
  geom_line(data = predictions, mapping = aes(x = time, y = n1, group = id), color = 'blue', alpha = 0.01) +
  xlab('') +
  ylab('Density of Blepharisma sp. (indiv./mL)')

panel_b = ggplot(filtered_data) +
  geom_point(mapping = aes(x = t, y = n2), color = 'red') +
  geom_line(data = predictions, mapping = aes(x = time, y = n2, group = id), color = 'red', alpha = 0.01) +
  xlab('Time (h)') +
  ylab('Density of Tetrahymena sp. (indiv./mL)')

ggsave('out/posterior_predictions.pdf', ggarrange(panel_a, panel_b, nrow = 2), width = 5, height = 6)

# plotting the median and 90% quantiles
panel_a = ggplot(filtered_data) +
  geom_point(mapping = aes(x = t, y = n1), color = 'blue') +
  stat_summary(data = predictions, mapping = aes(x = time, y = n1),
               fun.min = function(x) quantile(x, 0.05),
               fun.max = function(x) quantile(x, 0.95),
               geom = 'ribbon', fill = 'blue', alpha = 0.2) +
  stat_summary(data = predictions, mapping = aes(x = time, y = n1),
               fun = median,
               geom = 'line', color = 'blue') +
  xlab('') +
  ylab('Density of Blepharisma sp. (indiv./mL)')

panel_b = ggplot(filtered_data) +
  geom_point(mapping = aes(x = t, y = n2), color = 'red') +
  stat_summary(data = predictions, mapping = aes(x = time, y = n2),
               fun.min = function(x) quantile(x, 0.05),
               fun.max = function(x) quantile(x, 0.95),
               geom = 'ribbon', fill = 'red', alpha = 0.2) +
  stat_summary(data = predictions, mapping = aes(x = time, y = n2),
               fun = median,
               geom = 'line', color = 'red') +
  xlab('Time (h)') +
  ylab('Density of Tetrahymena sp. (indiv./mL)')

ggsave('out/posterior_predictions_quantiles.pdf', ggarrange(panel_a, panel_b, nrow = 2), width = 5, height = 6)