## prior predictive distribution
rm(list=ls())

# load packages and data #####################################################
library(rstan)
library(deSolve)
library(coda)
library(dplyr)

density_data = read.table("./Data/density_data.csv", header = T, sep = ',') ## reading data

filtered_data <- 
  density_data %>%
  filter(community == 'Ble_Tet') %>%
  select(c(measure_point, hours, Ble, Tet)) %>%
  group_by(measure_point) %>%
  summarise(t = mean(hours), n1 = mean(Ble), n2 = mean(Tet)) %>%
  arrange(t)

# functions ##################################################################
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

draw_random_parameters = function(){
  params['r1'] = rlnorm(1, -3, 1)
  params['r2'] = rlnorm(1, -2, 1)
  params['K1'] = rlnorm(1, 5, 0.8)
  params['K2'] = rlnorm(1, 8, 0.5)
  params['alpha1'] = rlnorm(1, -3, 1)
  params['alpha2'] = rlnorm(1, 2.5, 1)
  params['n10sim'] = abs(rnorm(1, 1.64,5))
  params['n20sim'] = rnorm(1, 510,10)
  return(params)
}

n_prior_pred = 200
times = seq(min(filtered_data$t), max(filtered_data$t), length.out = 200)
for (k in 1:n_prior_pred){
  print(k)
  par = draw_random_parameters()
  sim = ode(unlist(c(par['n10sim'], par['n20sim'])),
            times, ode.model, par,
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
  geom_point(mapping = aes(x = t, y = n2/10), color = 'red') +
  geom_line(data = predictions, mapping = aes(x = time, y = n1, group = id), color = 'blue', alpha = 0.01) +
  geom_line(data = predictions, mapping = aes(x = time, y = n2/10, group = id), color = 'red', alpha = 0.01) +
  theme_classic() +
  ylim(0, 1000)

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

