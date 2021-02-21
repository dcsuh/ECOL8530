## Daniel Suh
## Feb. 4 2021
##ODE's for COVID-19 vaccination protocol


library(tidyverse)
library(magrittr)
library(deSolve)


de<-function(t,x,params){
  S <- x[1]
  E <- x[2]
  Ik <- x[3]
  Iu <- x[4]
  Rk <- x[5]
  Ru <- x[6]
  V <- x[7]
  with(as.list(params),{
    dS <- -beta*S*Ik - beta*S*Iu - epsilon*(1-phi)*S
    dE <- beta*S*Ik + beta*S*Iu - upsilon*psi*E - upsilon*(1-psi)*E
    dIk <- upsilon*psi*E - alpha*Ik
    dIu <- upsilon*(1-psi)*E - alpha*Iu
    dRk <- alpha*Ik + epsilon*phi*Ru
    dRu <- alpha*Iu - epsilon*phi*Ru
    dV <- epsilon*(1-phi)*S
    res<-c(dS,dE,dIk,dIu,dRk,dRu,dV)
    list(res)
  })}
maxTime <- 100.0 # time
times<-seq(0,maxTime,by=0.1) # how often this calculates
# notes on params
# beta =      <- contact transmission rate
# phi =       <- range:0-1 probability seropositive
# alpha =     <- recovery rate
# epsilon =   <- antibody testing rate
# upsilon =   <- infected testing rate
# psi =       <- range: 0-1 probability infected
params<-c(beta=0.0003,
          phi=0.5,
          alpha=1/7,
          epsilon=0.1,
          upsilon=0.1,
          psi=0.5)  # model parameters

R0 <-(params[1]/params[3]) #calculate R0 with current parameters

xstart<-c(S=10000,
          E=0,
          Ik=1,
          Iu=0,
          Rk=0,
          Ru=0,
          V=0)  # initial conditions

output<-as.data.frame(lsoda(xstart,times,de,params)) # tells computer to solve (integrate) equations
output %<>% rowwise() %>% mutate(total = sum(c(S,E,Ik,Iu,Rk,Ru,V)))
output$X <- output$Rk + output$V
output$HI_threshold <- (1-1/R0)*xstart[1] #add in herd immunity threshold
output %>% ggplot(.,aes(x=time))+
  geom_line(aes(y=S,col="Susceptible"))+
  geom_line(aes(y=E,col="Exposed"))+
  geom_line(aes(y=Ik,col="Infected_k"))+
  geom_line(aes(y=Iu,col="Infected_u"))+
  geom_line(aes(y=Rk,col="Recovered_k"))+
  geom_line(aes(y=Ru,col="Recovered_u"))+
  geom_line(aes(y=V,col="Vaccinated"))+
  geom_line(aes(y=X,col="X"))+
  scale_colour_manual(values = c("red","orange","yellow","green","blue","purple","gray","black"))+
  theme(legend.position = c(0.75,0.4))+
  labs(y="N",x="Time",col="Population")+
  ggtitle("Transmission Model for Antibody-Assisted Vaccinations")+
  theme(plot.title = element_text(hjust=, size=10))+
  geom_hline(yintercept=output$HI_threshold)


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

