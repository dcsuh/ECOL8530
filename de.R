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
  V <- x[4]
  with(as.list(params),{
    dS <- -beta*S*I - epsilon*(phi*S) - epsilon*(1-phi)*S
    dI <- beta*S*I - alpha*I
    dR <- alpha*I + epsilon*phi*S
    dV <- epsilon*(1-phi)*S
    res<-c(dS,dI,dR,dV)
    list(res)
  })}
maxTime <- 20.0 # time
times<-seq(0,maxTime,by=0.1) # how often this calculates
# notes on params
# beta =      <- contact transmission rate
# phi =       <- range:0-1 probability seropositive
# alpha =     <- recovery rate
# epsilon =   <- testing rate
params<-c(beta=0.0003,
          phi=0.5,
          alpha=1/7,
          epsilon=0.1)  # model parameters

xstart<-c(S=10000,
          I=1,
          R=0,
          V=0)  # initial conditions

output<-as.data.frame(lsoda(xstart,times,de,params)) # tells computer to solve (integrate) equations
output$X <- output$R + output$V
output %>% ggplot(.,aes(x=time))+
  geom_line(aes(y=S,col="Susceptible"))+
  geom_line(aes(y=I,col="Infected"))+
  geom_line(aes(y=R,col="Recovered"))+
  geom_line(aes(y=V,col="Vaccinated"))+
  geom_line(aes(y=X,col="X"))+
  scale_colour_manual(values = c("red","forestgreen","blue","orange","black"))+
  theme(legend.position = c(0.75,0.5))+
  labs(y="N",x="Time",col="Population")+
  ggtitle("Transmission Model for Antibody-Assisted Vaccinations at an antibody testing rate of 0.1 and probabiliy of previous infection at 0.5")+
  theme(plot.title = element_text(hjust=0, size=10))+
  geom_hline(yintercept=0.5*10000)


de_1<-function(t,x,params){
  S <- x[1]  
  I <- x[2]
  R <- x[3]
  V <- x[4]
  with(as.list(params),{
    dS <- -beta*S*I - epsilon*S
    dI <- beta*S*I - alpha*I
    dR <- alpha*I
    dV <- epsilon*S
    res<-c(dS,dI,dR,dV)
    list(res)
  })}
maxTime <- 20.0 # time
times<-seq(0,maxTime,by=0.1) # how often this calculates
# notes on params
# beta =      <- contact transmission rate
# phi =       <- range:0-1 probability seropositive
# alpha =     <- recovery rate
# epsilon =   <- vaccination rate
params<-c(beta=0.0003,
          phi=0.5,
          alpha=1/7,
          epsilon=0.05)  # model parameters

xstart<-c(S=10000,
          I=1,
          R=0,
          V=0)  # initial conditions

output_1<-as.data.frame(lsoda(xstart,times,de_1,params)) # tells computer to solve (integrate) equations

output_1$X <- output_1$R + output_1$V
output_1 %>% ggplot(.,aes(x=time))+
  geom_line(aes(y=S,col="Susceptible"))+
  geom_line(aes(y=I,col="Infected"))+
  geom_line(aes(y=R,col="Recovered"))+
  geom_line(aes(y=V,col="Vaccinated"))+
  geom_line(aes(y=X,col="X"))+
  scale_colour_manual(values = c("red","forestgreen","blue","orange","black"))+
  theme(legend.position = c(0.75,0.5))+
  labs(y="N",x="Time",col="Population")+
  ggtitle("Transmission Model for blind Vaccinations at a vaccination rate of 0.05")+
  theme(plot.title = element_text(hjust=0, size=10))+
  geom_hline(yintercept=0.5*10000)

