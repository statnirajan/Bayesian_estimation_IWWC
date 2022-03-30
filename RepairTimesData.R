######## ############
######An algorithm to evaluate least squares estimators for Repair Times data
#install.packages("statsr")
#library(readxl)
library(dplyr)
library(statsr)
setwd("/Users/nirajanbudhathoki/Documents/RA")

x1<-read.csv("RepairTimesData.csv",header = T)
attach(x1)
n <- length(Xa)
gammahat <-rep(NA,n)
thetahat <- rep(NA,n)
yi <- log(Xa)
for (m in 1:(n-1)){
  x1<-x1 %>%
    mutate(X1i=-log(log(((119+0.4)*2.718)/(2*(a-0.3)))))
  x1<-x1 %>%
    mutate(X2i=log(log(((119+0.4)*2.718)/(2*(119-a+0.7)))))
  x1<-x1 %>%
    mutate(yi)
  
  c1=subset(x1,a<=m,select=c(X1i,yi))
  
  sum(c1$X1i)
  
  c2=subset(x1,a>m,select=c(X2i,yi))
  
  sum(c2$X2i)
  xbar=(sum(c1$X1i)+sum(c2$X2i))/119
  xbar
  ybar=sum(yi)/119
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
m_position                    #59
thetahat[m_position]          #8.07917  #LS Estimate of Theta
gammahat[m_position]          #0.78345  #LS Estimate of Gamma

#### Bayesian Estimation with Jeffrey's prior ##################################
##Inference on both gamma and theta simultaneously
set.seed(10004)
xdata=Xa
gamma=0.78345   #initial values from LS estimates
theta=8.07917
N=10000
m= m_position
n=length(Xa)
#log of full conditional/posterior for theta>0
lnf.theta<-function(theta,xdata,n,m,gamma){
  ((2*m*gamma-n*gamma-1)*log(theta))-sum((theta/xdata[1:m])^gamma)-sum((xdata[(m+1):n]/theta)^gamma)
}
#log of the full conditional/posterior for gamma>0
lnf.gamma<-function(gamma,xdata,n,m,theta){
  ((2*m*gamma-n*gamma-1)*log(theta))+(n*log(gamma))+(gamma*log(prod(xdata[(m+1):n])))-(gamma*log(prod(xdata[1:m])))-
    sum((theta/xdata[1:m])^gamma)-sum((xdata[(m+1):n]/theta)^gamma)
}

parm<-list(xdata=xdata,theta=theta,gamma=gamma,n=n,m=m)
rmc1 <- rmc2<- gammamc <-thetamc <-alpha1 <-alpha2 <-rep(NA, N) # Vectors to store results
last= 0.5   #for gamma
last2=7.0  #for theta
set.seed(9099)
for (j in 1:N){
  ##For gamma
  repeat{ cand<-rnorm(1, mean=last, sd=sqrt(0.01))  # Proposing normal distribution
  if (cand>0) break }
  alpha1[j] <-exp(lnf.gamma(cand, xdata=parm$xdata, n=parm$n, m=parm$m, theta=parm$theta)-
                    lnf.gamma(last,xdata=parm$xdata, n=parm$n, m=parm$m, theta=parm$theta))
  
  rmc1[j]<-r1<-rbinom( n = 1,size = 1, prob = min(alpha1[j],1))
  if(r1==1) last <-cand
  gammamc[j]<-gamma<-last
  parm<- list(xdata=xdata, gamma=gamma, theta=theta, n = n, m=m)
  ##For theta
  
  repeat{ cand2<-rnorm(1, mean=last2, sd=sqrt(2))  # Proposing normal distribution
  if (cand2>0) break }
  alpha2[j] <-exp(lnf.theta(cand2, xdata=parm$xdata, n=parm$n, m=parm$m, gamma=parm$gamma)-
                    lnf.theta(last2,xdata=parm$xdata, n=parm$n, m=parm$m, gamma=parm$gamma))
  
  rmc2[j]<-r2<-rbinom( n = 1,size = 1, prob = min(alpha2[j],1))
  if(r2==1) last2 <-cand2
  thetamc[j]<-theta<-last2
  parm<- list(xdata=xdata, gamma=gamma, theta=theta, n = n, m=m)
}
jx<-2001:N
sum(rmc1[jx])/length(jx)        #acceptance rate for gamma  45.4%
sum(rmc2[jx])/length(jx)        #acceptance rate for theta  46.9%
plot(gammamc[jx],type="l")
plot(thetamc[jx],type="l")
mean(gammamc[jx])               #posterior estimate of gamma is 0.79597
mean(thetamc[jx])               #posterior estimate of theta is 7.94524
quantile(gammamc[jx],c(0.025,0.975))
quantile(thetamc[jx],c(0.025,0.975))

###ASSESSING GOODNESS OF FITS ###########
###ANDERSON_DARLING TEST STATISTICS #######
n <- length(Xa)
theta=  thetahat[m_position]      # mean(thetamc[jx])          thetahat[m_position]                     
gamma=  gammahat[m_position]    # mean(gammamc[jx])          #gammahat[m_position]       
cdf=rep(NA,n)
for (i in 1:n){
  if (Xa[i]<=theta){
    cdf[i]=(exp(1)/2)*exp(-(theta/Xa[i])^gamma)
  } else{
    cdf[i]=1-((exp(1)/2)*exp(-(Xa[i]/theta)^gamma))
  }
}
ln.cdf1=log(cdf)
decr_cdf=sort(cdf,decreasing = TRUE)
ln.cdf2=log(1-(decr_cdf))
sum1=ln.cdf1+ln.cdf2
S=sum((((2*a)-1)/n)*sum1)
AD=-n-S
print(AD)   #0.5622 LS   #0.6312 Bayes

# Using goftest library in R since it gives p-value as well.
library(goftest)
mycdf <- function(Xa,theta,gamma){
  cdf=rep(NA,n)
  for (i in 1:n){
    if (Xa[i]<=theta){
      cdf[i]=(exp(1)/2)*exp(-(theta/Xa[i])^gamma)
    } else{
      cdf[i]=1-((exp(1)/2)*exp(-(Xa[i]/theta)^gamma))
    }
  }
  return(cdf)
}
mycdf(Xa,theta,gamma)
ad.test(Xa,mycdf,theta,gamma)   # Test statistic: 0.56224, p-value = 0.684
ad.test(Xa,mycdf,theta=mean(thetamc[jx]),gamma=mean(gammamc[jx]))

#### KOLMOGOROV-SMIRNOV TEST STATISTICS ######
ks.test(Xa,mycdf,theta,gamma)
ks.test(Xa,mycdf,theta=mean(thetamc[jx]),gamma=mean(gammamc[jx]))

# This is given as Max|Fobs(Xi)- Fexp(Xi)|####
n=length(Xa)
theta= mean(thetamc[jx])           #thetahat[m_position]
gamma= mean(gammamc[jx])           # gammahat[m_position]  
table(Xa)
t.1=as.data.frame(table(Xa))
t.1<-t.1%>%
  mutate(Cum.freq=cumsum(t.1$Freq))
obs.dist.fxn=t.1$Cum.freq/n
exp.dist.fxn=rep(NA,length(t.1$Xa))
attach(t.1)
Xn=as.numeric(Xa)
k=length(Xn)
for (i in 1:k){
  if (Xn[i]<=theta){
    exp.dist.fxn[i]=(exp(1)/2)*exp(-(theta/Xn[i])^gamma)
  } else{
    exp.dist.fxn[i]=1-((exp(1)/2)*exp(-(Xn[i]/theta)^gamma))
  }
}
abs.diff=abs(obs.dist.fxn-exp.dist.fxn)
print(max(abs.diff))  #0.0658 LS   #0.0666 Bayes



##Akaike Information Criteria#####
# AIC = -2*log.like + 2*p  where p is #of parameters estimated, 2 in our case.
n=119
m=58
gamma= mean(gammamc[jx])            #gammahat[m_position]        
theta= mean(thetamc[jx])            #thetahat[m_position]      
xdata=as.numeric(x1$Xa)
log.like=n*log(exp(1)/2)-n*log(gamma)-sum(log(xdata[1:n]))+gamma*(sum(log(xdata[(m+1):n])))-gamma*(sum(xdata[1:m]))+
  ((2*m*gamma-n*gamma)*log(theta))-sum((theta/xdata[1:m])^gamma)-sum((xdata[(m+1):n]/theta)^gamma)
AIC=-2*log.like + 2*2
print(AIC)         #978.176 LS      #989.996 Bayes


## Standard Errors
amhat = 1/gammahat[m_position]
bmhat = log(thetahat[m_position])
pred.y1 = amhat*x1$X1i[1:59] + bmhat
pred.y2 = amhat*x1$X2i[60:119] + bmhat
pred = c(pred.y1,pred.y2)
pred
res = x1$yi - pred
sig = sqrt((sum(res[1:59]^2) + sum (res[60:119]^2))/117)

l.x1=x1$X1i[1:59]-mean(c(x1$X1i[1:59],x1$X2i[60:119]))
s1 = sum(l.x1^2)
l.x2=x1$X2i[60:119]-mean(c(x1$X1i[1:59],x1$X2i[60:119]))
s2 = sum(l.x2^2)
sqrt(sum(s1,s2))
se_amhat = sig/sqrt(sum(s1,s2))  # This is SE(amhat)
(se_gammahat = se_amhat/(amhat^2)) #This is SE(gammamhat) # 0.00796

# SE(thetahat)

se_bmhat = sqrt((sum((x1$X1i[1:59])^2) + sum((x1$X2i[60:119])^2)) / (n * sum(s1,s2)))*sig
(se_thetahat = exp(bmhat) * se_bmhat)   #This is SE(thetahat) # 0.07548









