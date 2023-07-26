## A stan model using mono, bi and trispecific time series to fit demographic parameters on Tet, Col and Ble.

rm(list=ls())

# load packages ##############################################################
library(rstan)
library(deSolve)
library(coda)
library(dplyr)

dir.create("out", showWarnings = FALSE) ## create a directory to store outputs

## Declare Stan model ########################################################

model_str = '
functions{
  real[] odemodel(real t, real[] N, real[] p, real[] x_r, int[] x_i){
    // p[1]=r1, p[2] = r2, p[3] = alpha11, p[4] = alpha21, p[5] = alpha12, p[6] = alpha22 
    real dNdt[2]; 
    dNdt[1] = p[1]*N[1] - p[3]*N[1]^2 - p[4]*N[1]*N[2];
    dNdt[2] = p[2]*N[2] - p[5]*N[1]*N[2] - p[6]*N[2]^2;
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
  real<lower=0> alpha11; // density dependence species 1
  real<lower=0> alpha21; // density dependence species 2 on 1
  real<lower=0> alpha12; // density dependence species 1 on 2
  real<lower=0> alpha22; // density dependence species 2
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
  alpha11 ~ lognormal(-9, 1);
  alpha21 ~ lognormal(-9, 1);
  alpha12 ~ lognormal(-8, 1);
  alpha22 ~ lognormal(-8, 1);
  n10sim ~ normal(n1[1],5);
  n20sim ~ normal(n2[1],10);
  sdev1 ~ gamma(1, 1);
  sdev2 ~ gamma(1, 1);
  
  // parameters for integrator
  p[1] = r1;
  p[2] = r2;
  p[3] = alpha11;
  p[4] = alpha21;
  p[5] = alpha12;
  p[6] = alpha22;

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
}
'

s_model = stan_model(model_code=model_str) ## compiling the model


density_data = read.table("./Data/density_data.csv", header = T, sep = ',') ## reading data

filtered_data <- 
density_data %>%
  filter(community == 'Ble_Tet') %>%
  select(c(measure_point, hours, Ble, Tet)) %>%
  group_by(measure_point) %>%
  summarise(t = mean(hours), n1 = mean(Ble), n2 = mean(Tet)) %>%
  arrange(t)

ggplot(data = filtered_data) +
  geom_point(mapping = aes(x = t, y = n1))

ggplot(data = filtered_data) +
  geom_point(mapping = aes(x = t, y = n2))


data = list(n  = nrow(filtered_data),
            
            t = filtered_data$t,
            
            n1 = filtered_data$n1,
            n2 = filtered_data$n2)

# stan options
chains = 3
rstan_options(auto_write = TRUE)
options(mc.cores = chains)

iter   =  1000
warmup =  200
thin   =     1

# initial values for sampling 
init=rep(list(list(r1=0.01,
                   r2=0.1,
                   alpha11 = 10^-4,
                   alpha21 = 10^-4,
                   alpha12 = 10^-3,
                   alpha21 = 10^-3,
                   n10sim=5,
                   n20sim=100,
                   sdev1 = 1,
                   sdev2 = 1
))
,chains)


# run model and print result
fit_obs = sampling(s_model,
                   data=data,
                   iter=iter,
                   warmup=warmup,
                  thin=thin,
                  chains=chains,
                  init=init,
                  control = list(adapt_delta = 0.9, max_treedepth=12),
                  refresh=10
)
 
save(fit_obs, file="./out/fit_posterior.RData")

#load("./out/fit_posterior.RData")

print(fit_obs)

samples=As.mcmc.list(fit)
params = c("r1","r2", "alpha11", "alpha21", "alpha12", "alpha21")
plot(samples[, params])

pairs(fit_obs, pars=params)

ode.model = function(t,N,p){
  r1 = p$r1
  r2 = p$r2
  alpha11 = p$alpha11
  alpha21 = p$alpha21
  alpha12 = p$alpha12
  alpha22 = p$alpha21
  dn1 = r1*N[1] - alpha11*N[1]**2 - alpha21*N[1]*N[2]
  dn2 = r2*N[2] - alpha12*N[1]*N[2] - alpha22*N[2]**2
  return(list(c(dn1, dn2)))
}


posteriors = as.matrix(fit_obs)

n_post = 1000
times = seq(min(data$t), max(data$t), length.out = 200)
for (k in 1:n_post){
  par = posteriors[sample(1:nrow(posteriors), 1),]
  sim = ode(c(par['n10sim'], par['n20sim']),
            times, ode.model, list(r1 = par['r1'],
                                   r2= par['r2'],
                                   alpha11 = par['alpha11'],
                                   alpha21 = par['alpha21'],
                                   alpha12 = par['alpha12'],
                                   alpha22 = par['alpha22']),
            method = 'ode45')
  
  temp  = data.frame(time = sim[,1], n1 = sim[,2], n2 = sim[,3], id = k)
  
  if (k == 1) {
    predictions = temp
  } else {
    predictions = rbind(predictions, temp)
  }
  
}


ggplot(filtered_data) +
  geom_point(mapping = aes(x = t, y = n1), color = 'blue') +
  geom_point(mapping = aes(x = t, y = n2), color = 'red') +
  geom_line(data = predictions, mapping = aes(x = time, y = n1, group = id), color = 'blue', alpha = 0.01) +
  geom_line(data = predictions, mapping = aes(x = time, y = n2, group = id), color = 'red', alpha = 0.01) +
  theme_classic()

ggplot(filtered_data) +
  geom_point(mapping = aes(x = t, y = n1), color = 'blue') +
  geom_point(mapping = aes(x = t, y = n2), color = 'red') +
  stat_summary(data = predictions, mapping = aes(x = time, y = n1),
               fun.min = function(x) quantile(x, 0.05),
               fun.max = function(x) quantile(x, 0.95),
               geom = 'ribbon', fill = 'blue', alpha = 0.2) +
  stat_summary(data = predictions, mapping = aes(x = time, y = n1),
               fun = median,
               geom = 'line', color = 'blue') +
  stat_summary(data = predictions, mapping = aes(x = time, y = n2),
               fun.min = function(x) quantile(x, 0.05),
               fun.max = function(x) quantile(x, 0.95),
               geom = 'ribbon', fill = 'red', alpha = 0.2) +
  stat_summary(data = predictions, mapping = aes(x = time, y = n2),
               fun = median,
               geom = 'line', color = 'red') +
  theme_classic()
