############ GIBBS SAMPLER WITH JEFFREY'S PRIOR #################

library(HDInterval)
setwd("K:/Departments/STA/SCC/Budhathoki")
library(ars)
x1<-read.csv("sim1.csv",header = T)
attach(x1)
xdata=Xa
k=149
n1=length(Xa)

gamma1 = 2

parm = list(xdata=xdata, n1=n1, k=k, gamma1=gamma1)

# we need a function that computes log(f(u)),for given u. f(u) is proportional to the density we want to sample from.
###### log of full conditional for theta > 0 ######

f1 = function(x,xdata,n1,k,gamma1){
  (2*k*gamma1-n1*gamma1-1)*log(x) - sum((x/xdata[1:k])^gamma1) - sum((xdata[(k+1):n1]/x)^gamma1)
}

#### fprima = d/du log(f(u)) ####

fprima1 = function(x,xdata,n1,k,gamma1){
  ((2*k*gamma1-n1*gamma1-1)/x) - (gamma1*x^(gamma1-1))*sum((1/xdata[1:k])^gamma1) + gamma1*(x^(-gamma1-1))*sum((xdata[(k+1):n1])^gamma1)                         
}





##### log of full conditional for gamma > 0
theta1 = 8


f2 = function(x,xdata,n1,k,theta1){
  (2*k*x-n1*x-1)*log(theta1) + n1*log(x) + x*sum(log(xdata[(k+1):n1])) - x*sum(log(xdata[1:k]))-sum((theta1/xdata[1:k])^x)-
    sum((xdata[(k+1):n1]/theta1)^x)
}


fprima2 = function(x,xdata,n1,k,theta1){
  (2*k-n1)*log(theta1) + (n1/x) + sum(log(xdata[(k+1):n1])) - sum(log(xdata[1:k]))-sum(((theta1/xdata[1:k])^x)*log(theta1/xdata[1:k]))-
    sum(((xdata[(k+1):n1]/theta1)^x)*log((xdata[(k+1):n1]/theta1)))
}



# install.packages('tidyverse')
library(tidyverse)
x1 <- read.csv("sim1.csv")   #choose sim.csv
attach(x1)
n <- length(Xa)
gammahat <-rep(NA,n)
thetahat <- rep(NA,n)
yi <- log(Xa)
for (m in 1:(n-1)){
  
  x1<-x1 %>%
    mutate(X1i=-log(log(((300+0.4)*2.718)/(2*(a-0.3)))))
  x1<-x1 %>%
    mutate(X2i=log(log(((300+0.4)*2.718)/(2*(300-a+0.7)))))
  x1<-x1 %>%
    mutate(yi)
  
  c1=subset(x1,a<=m,select=c(X1i,yi))
  
  sum(c1$X1i)
  
  c2=subset(x1,a>m,select=c(X2i,yi))
  
  sum(c2$X2i)
  xbar=(sum(c1$X1i)+sum(c2$X2i))/300
  xbar
  ybar=sum(yi)/300
  ybar
  c1<-c1 %>%
    mutate(up1=(X1i-xbar)*(yi-ybar))
  c2<-c2 %>%
    mutate(up2=(X2i-xbar)*(yi-ybar))
  up=sum(c1$up1)+sum(c2$up2)
  c1<-c1 %>%
    mutate(up11=(X1i-xbar)^2)
  c2<-c2 %>%
    mutate(up22=(X2i-xbar)^2)
  lw=sum(c1$up11)+sum(c2$up22)
  amhat=up/lw
  gammahat[m]=1/amhat
  
  
  bmhat=ybar-amhat*xbar
  thetahat[m]=exp(bmhat)
  
}
cbind(x1[-n],thetahat)
K=rep(NA,n)
Z<-as.data.frame(x1)
Xa
for (i in 1:(n-1)){
  if (Xa[i] < thetahat[i] & thetahat[i] < Xa[i+1]){
    K[i]<-1
  }else{
    K[i]<-0
  }
}
K
m_position=which(K %in% 1)
m_position
thetahat[m_position]
gammahat[m_position]


##Inference on both gamma and theta simultaneously
set.seed(1994)
xdata <- Xa
theta1 <- 8#thetahat[m_position]    #initial values, these are logarithm of LS estimates
gamma1 <- 2 #gammahat[m_position]
N=10000
k = 149#m_position
n1=300



parm<-list(xdata=xdata,theta1=theta1,gamma1=gamma1,k=k,n1=n1)
theta1mc <-gamma1mc <-alpha1 <-alpha2 <-rep(NA, N) # Vectors to store results
#last= 7.0   #for theta1 
#last2=1.5  #for gamma1
set.seed(9199)
for (j in 1:N){
  ##For gamma1
  cand2_100 = ars(5, f2, fprima2, x= 2, m=1, lb=TRUE, xlb=0, ub = TRUE, xub=100, xdata=parm$xdata, n1=parm$n1, k=parm$k, theta1=parm$theta1)
  gamma1mc[j]<-gamma1<-cand2_100[5]
  parm<- list(xdata=xdata,theta1=theta1, gamma1=gamma1,k=k, n1=n1)
  
  ##For theta1
  
  cand_100 = ars(5, f1, fprima1, x=7, m=1, lb=TRUE, xlb=0, ub=TRUE, xub=100, xdata=parm$xdata, n1=parm$n1, k=parm$k, gamma1=parm$gamma1)
  theta1mc[j]<-theta<-cand_100[5]
  parm<- list(xdata=xdata, theta1=theta1, gamma1=gamma1,k=k, n1=n1)
}

jx<-2001:N
# When do we accept and when do we reject ???????????????????????????
#sum(theta1mc[jx])/length(jx)        #acceptance rate for eta.2  
#sum(gamma1mc[jx])/length(jx)        #acceptance rate for eta.1  
#gammamc <- exp(eta.2mc)
#thetamc <- exp(eta.1mc)
plot(gamma1mc[jx],type = "l")
plot(theta1mc[jx],type = "l")
mean(gamma1mc[jx])
mean(theta1mc[jx])
quantile(gamma1mc[jx],c(0.025,0.975))
quantile(theta1mc[jx],c(0.025,0.975))
hdi(gamma1mc[jx])
hdi(theta1mc[jx])
