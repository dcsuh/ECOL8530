## Daniel Suh
## Feb. 4 2021
##ODE's for COVID-19 vaccination protocol


library(tidyverse)
library(magrittr)
library(deSolve)


de<-function(t,x,params){
  S <- x[1]  
  I <- x[2]
  R <- x[3]
  A <- x[4]
  V <- x[5]
  with(as.list(params),{
    dS <- -beta*S*I - epsilon*S
    dI <- beta*S*I - alpha*I
    dR <- alpha*I + phi*A
    dA <- epsilon*S - phi*A - (1-phi)*A
    dV <- (1-phi)*A
    res<-c(dS,dI,dR,dA,dV)
    list(res)
  })}
maxTime <- 20.0 # time is in years - run model for this time
times<-seq(0,maxTime,by=0.01) # how long we run the model for
# notes on params
# beta =       <- contact transmission rate
# phi =      <- range:0-1 probability infected
# alpha =     <- recovery rate
# epsilon =       <- antibody testing rate
params<-c(beta=0.03,
          phi=0.1,
          alpha=6,
          epsilon=1)  # model parameters

xstart<-c(S=100,
          I=1,
          R=0,
          A=0,
          V=0)  # initial conditions

output<-as.data.frame(lsoda(xstart,times,de,params)) # tells computer to solve (integrate) equations

output %>% ggplot(.,aes(x=time))+
  geom_line(aes(y=S,col="Susceptible"))+
  geom_line(aes(y=I,col="Infected"))+
  geom_line(aes(y=R,col="Recovered"))+
  geom_line(aes(y=A,col="Antibody"))+
  geom_line(aes(y=V,col="Vaccinated"))+
  scale_colour_manual(values = c("red","forestgreen","blue","yellow","orange"))+
  theme(legend.position = c(0.75,0.5))+
  labs(y="N",x="Time",col="Population")+
  ggtitle("Transmission Model for COVID19 Vaccinations")+
  theme(plot.title = element_text(hjust=0, size=10))
