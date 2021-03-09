#script originator: Daniel Suh
#date: 3/4/21


#latin hypercube sampling for eqs from eqs_ode.Rmd

library(lhs)
library(sensitivity)

de<-function(t,x,params){
  S <- x[1]
  E <- x[2]
  Ik <- x[3]
  Iu <- x[4]
  Rk <- x[5]
  Ru <- x[6]
  V <- x[7]
  with(as.list(params),{
    dS <- -beta1*S*Ik - beta2*S*Iu - phi*gamma - (1-phi)*gamma*(S/(S+E+Iu+Ru)) + lambda - mu*S
    dE <- beta1*S*Ik + beta2*S*Iu - psi*theta*E - (1-psi)*theta*E - (1-phi)*gamma*(E/(S+E+Iu+Ru)) - mu*E
    dIk <- psi*theta*E - alpha*Ik - mu1*Ik
    dIu <- (1-psi)*theta*E - alpha*Iu - (1-phi)*gamma*(Iu/(S+E+Iu+Ru)) - mu1*Iu
    dRk <- alpha*Ik - mu*Rk
    dRu <- alpha*Iu - (1-phi)*gamma*(Ru/(S+E+Iu+Ru)) - mu*Ru
    dV <- gamma - mu*V
    res<-c(dS,dE,dIk,dIu,dRk,dRu,dV)
    list(res)
  })}
maxTime <- 20.0 # time
times<-seq(0,maxTime,by=0.5) # how often this calculates
# notes on params
# beta =      <- contact transmission rate
# lambda =    <- constant growth rate
# mu =        <- density-dependent mortality rate
# phi =       <- range: 0-1 proportion vaccines targeted to known susceptible individuals
# psi =       <- range: 0-1 probability infected individuals are identified
# alpha =     <- recovery rate
# gamma =     <- number of daily vaccines

params<-c(beta1=0.000003,
          beta2=0.00003,
          lambda=1/2,
          mu=1/10,
          mu1=1/5,
          phi=0.5,
          psi=0.3,
          theta=1/5,
          alpha=1/7,
          gamma=10000)  # model parameters

xstart<-c(S=5*10^6,
          E=1.6*10^5,
          Ik=1.3*10^5,
          Iu=4.9*10^5,
          Rk=6.6*10^5,
          Ru=2.7*10^6,
          V=0)  # initial conditions


h <- 1000
set.seed(8878896)
lhs <- maximinLHS(h,10)


beta1.min=0.00000003
beta1.max=0.0000003
beta2.min=0.00000003
beta2.max=0.0000003
lambda.min=5
lambda.max=20
mu.min=1/1000
mu.max=1/100
mu1.min=1/100
mu1.max=1/10
phi.min=0
phi.max=1
psi.min=0
psi.max=1
theta.min=1/7
theta.max=1/4
alpha.min=1/7
alpha.max=1/14
gamma.min=10000
gamma.max=50000


params.set <- cbind(
  beta1 = lhs[,1]*(beta1.max-beta1.min)+beta1.min,
  beta2 = lhs[,2]*(beta2.max-beta2.min)+beta2.min,
  lambda = lhs[,3]*(lambda.max-lambda.min)+lambda.min,
  mu = lhs[,4]*(mu.max-mu.min)+mu.min,
  mu1 = lhs[,5]*(mu1.max-mu1.min)+mu1.min,
  psi = lhs[,6]*(psi.max-psi.min)+psi.min,
  theta = lhs[,7]*(theta.max-theta.min)+theta.min,
  alpha = lhs[,8]*(alpha.max-alpha.min)+alpha.min,
  gamma = lhs[,9]*(gamma.max-gamma.min)+gamma.min)


levels <- 11

h2 <- 1000

R0 <- 2
threshold <- (1-1/R0)*as.numeric(xstart[1])
antibody.prop <- seq(0,1,0.1)
j <- 1  
data <- data.frame(matrix(rep(NA,levels*h2*11),nrow=levels*h2))
for(i in 1:h2){
  for (k in 1:length(antibody.prop)){
    
    data[j,1:10] <- params <- as.list(c(params.set[i,], phi=antibody.prop[k]))
    output <- as.data.frame(ode(xstart, times, de, params))
    output$herd <- output$S<=threshold
    tmp <- output %>% filter(herd==1)
    data[j,11] <- tmp$time[1]
    j <- j+1
    
  }
}

names(data) <- c(names(params),'herd.threshold')

data %>% mutate(phi = as.factor(phi)) %>% ggplot(., aes(x=phi, y=herd.threshold)) + geom_boxplot(notch = T)


bonferroni.alpha = 0.05/10
prcc <- pcc(data[,1:9], data[,11], nboot=1000, rank=T, conf=1-bonferroni.alpha)
print(prcc)
