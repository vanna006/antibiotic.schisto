### basic schisto model for Hannah's experiment

library(deSolve)

#stuff from Chris Hoover
#Agrochem model parameters

area = 200

parameters <- c (
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = area*1.5,      # Human population at site of interest (based on 300 people at 200m^2 water contact site from Sokolow PNAS)
  #Om = 1 / sqrt(area),# degree of overlap between water contamination, snail, and human habitats
  
  # Snail reproductive parameters
  f_N = 0.6,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
  #   from Woolhouse & Chandiwana et al. 1990 
  K_N = 50*area,   # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  chi = 1.332, #fecundity compensation
  rho = 0.007, #castration
    
  # Snail mortality parameters
  mu_N = 1/60,        # Natural mortality rate of large snails (deaths/snail/day; mean lifespan ~ 60 days (~2 months))
  mu_I = 1/10-1/60,   # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  # miracidia parameters
  m = 432,           # Miracidial shedding rate per adult female worm assuming 0.36 eggs/mL urine and 1200 mL urine per person per day
  v = 0.084,         # Egg viability of S. haematobium (i.e. miracidia/egg) from Halstead et al
  pi_M = 6.22,       # mean Miracidia-hrs per day in agrochemical-free water
  
  # cercariae parameters
  theta = 15.6, # 1320,      # mean daily cercarial shedding rate of patently infected snails @25C (cercariae/I-snail/day); Pfluger 1984
  pi_C = 14.21,        # mean cercariae-hrs in agrochemical-free water
  
  # transmission parameters
  beta = 6.666667e-08,       # Human-to-snail infection probability in reference area (exposed snails/susceptible snail/miracidi-hr/day); 
  #     divided by 24 to account for hourly scale of miracidial survival; from fitting procedure in Halstead et al
  sigma = 0.025,         # Latent period for exposed snails (infectious snails/exposed snail/day))
  lambda = 7.5e-8, # Snail-to-human infection probability in reference area (adult worms/cercariae-hr); 
  #     divided by 24 to account for hourly scale of miracidial survival; from fitting procedure in Halstead et al
  
  # Schisto mortality parameters
  mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  mu_H = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
  )
state<-c(S = 8000, E = 1000, I = 200, W = 2)
  #base model code that contains no agrochemical influences or predators
  base_mod =function(t, state, parameters) { 
    with(as.list(c(state, parameters)),{
      
      #Dynamic variables
      N = S + E + I                                       # Total number of snails
      
      #Schistosome larval concentration equations
      M = 0.5*W*H*m*v*pi_M                                # Number of female worms
      C = theta*I*pi_C                                    # Cercariae concentration
      
      #differential equations
      
      dSdt = f_N*(1-(N/K_N))*(S + (chi*E) + (rho*I)) - mu_N*S - beta*M*S    # Susceptible snails
      
      dEdt = beta*M*S - mu_N*E - sigma*E                  # Exposed snails
      
      dIdt = sigma*E - (mu_N + mu_I)*I                    # Infected snails
      
      dWdt = lambda*C - (mu_W+mu_H)*W                     # Mean worm burden in human population
      
      
      return(list(c(dSdt, dEdt, dIdt, dWdt)))
      
    })
  }  
  
  times <- seq(0, 6000, by = 0.01)
  out <- ode(y = state, times = times, func = base_mod, parms = parameters)
  plot(out)
  
  
### model with antibiotic parameters
  area = 200
  
  parameters <- c (
    # Location parameters
    A = area,          # Area of site of interest, m^2
    H = area*1.5,      # Human population at site of interest (based on 300 people at 200m^2 water contact site from Sokolow PNAS)
    #Om = 1 / sqrt(area),# degree of overlap between water contamination, snail, and human habitats
    
    # Snail reproductive parameters
    f_N = 0.81,        # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
    #   from Woolhouse & Chandiwana et al. 1990 
    K_N = 50*area,   # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
    chi = 1.174, #fecundity compensation
    rho = 0.115, #castration
    
    # Snail mortality parameters
    mu_N = 0.015,        # Natural mortality rate of large snails (deaths/snail/day; mean lifespan ~ 60 days (~2 months))
    mu_I = 0.107,   # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
    
    # miracidia parameters
    m = 432,           # Miracidial shedding rate per adult female worm assuming 0.36 eggs/mL urine and 1200 mL urine per person per day
    v = 0.084,         # Egg viability of S. haematobium (i.e. miracidia/egg) from Halstead et al
    pi_M = 6.22,       # mean Miracidia-hrs per day in agrochemical-free water
    
    # cercariae parameters
    theta = 22.2,  # 1880,      # mean daily cercarial shedding rate of patently infected snails @25C (cercariae/I-snail/day); Pfluger 1984
    pi_C = 14.21,        # mean cercariae-hrs in agrochemical-free water
    
    # transmission parameters
    beta = 8.66e-7,       # Human-to-snail infection probability in reference area (exposed snails/susceptible snail/miracidi-hr/day); 
    #     divided by 24 to account for hourly scale of miracidial survival; from fitting procedure in Halstead et al
    sigma = 0.025,         # Latent period for exposed snails (infectious snails/exposed snail/day))
    lambda = 7.5e-8, # Snail-to-human infection probability in reference area (adult worms/cercariae-hr); 
    #     divided by 24 to account for hourly scale of miracidial survival; from fitting procedure in Halstead et al
    
    # Schisto mortality parameters
    mu_W = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
    mu_H = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
  )
  state<-c(S = 8000, E = 1000, I = 200, W = 2)
  #base model code that contains no agrochemical influences or predators
  base_mod =function(t, state, parameters) { 
    with(as.list(c(state, parameters)),{
      
      #Dynamic variables
      N = S + E + I                                       # Total number of snails
      
      #Schistosome larval concentration equations
      M = 0.5*W*H*m*v*pi_M                                # Number of female worms
      C = theta*I*pi_C                                    # Cercariae concentration
      
      #differential equations
      
      dSdt = f_N*(1-(N/K_N))*(S + (chi*E) + (rho*I)) - mu_N*S - beta*M*S    # Susceptible snails
      
      dEdt = beta*M*S - mu_N*E - sigma*E                  # Exposed snails
      
      dIdt = sigma*E - (mu_N + mu_I)*I                    # Infected snails
      
      dWdt = lambda*C - (mu_W+mu_H)*W                     # Mean worm burden in human population
      
      
      return(list(c(dSdt, dEdt, dIdt, dWdt)))
      
    })
  }  
  
  times <- seq(0, 6000, by = 0.01)
  out2 <- ode(y = state, times = times, func = base_mod, parms = parameters)
  plot(out2)
  
  
### create some plots with the output.
Sout <- out2[,2]/out[,2]
Eout <- out2[,3]/out[,3]
Iout <- out2[,4]/out[,4]
Wout <- out2[,5]/out[,5]
Nout <- (out2[,2] + out2[,3] + out2[,4]) / (out[,2] + out[,3] + out[,4])


S <- cbind(out[,1],Sout)
E <- cbind(out[,1],Eout)
I <- cbind(out[,1],Iout)
W <- cbind(out[,1],Wout)
N <- cbind(out[,1],Nout)

par(mfrow=c(2,2))
plot(S, type = 'l', lwd = 2, col = 'black', xlab = 'Days', ylab = 'Relative number of susceptible snails')
plot(E, type = 'l', lwd = 2, col = 'blue', xlab = 'Days', ylab = 'Relative number of exposed snails')
plot(I, type = 'l', lwd = 2, col = 'orange', xlab = 'Days', ylab = 'Relative number of infected snails')
plot(W, type = 'l', lwd = 2, col = 'red', xlab = 'Days', ylab = 'Relative worm burden')

par(mfrow=c(1,1))
plot(N, type = 'l', lwd = 2, col = 'black', xlab = 'Days', ylab = 'Relative snail population')


#code for smoothed lines
library(ggplot2)
ds <- as.data.frame(S)
de <- as.data.frame(E)
di <- as.data.frame(I)
dw <- as.data.frame(W)
dn <- as.data.frame(N)

ggplot(ds, aes(ds[,1], ds[,2])) + geom_point() +
  geom_smooth(method = "loess", span = 0.3, se = FALSE) 







