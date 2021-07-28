install.packages("ggplot2")
library(ggplot2)
install.packages("nlstools")
library(nlstools)
install.packages("lmvar")
library(lmvar)
install.packages("DAAG")
library(DAAG)
install.packages("dplyr")
library(dplyr)

setwd('/Users/mahsa/Desktop/MinimaxEstimationOfLargestSparseMean/MinimaxEstimationOfLargestSparseMean')
#Data = read.csv('All_Logit_n100.txt',header=TRUE, sep = ',')
Data = read.csv('All_Logit_n200.txt',header=TRUE, sep = ',')
Data_MSE <- readRDS(file="MinimaxMSECurves2000.vs.K.RDS")

### Modeling MSE and Epsilon=Kappa 
MSE = Data_MSE[,3]
EPS = Data_MSE[,2]
L = Data_MSE[,1]
Df <- data.frame(Date=as.Date(character()),
                       File=character(), 
                       User=character(), 
                       stringsAsFactors=FALSE) 
for (i in unique(L))
{
  loc1 = which(L==i)
  mse_filtered = MSE[loc1] ###MSE for each L frames separately
  eps_filtered = EPS[loc1] ###EPS for each L frames separately
  model = splinefun(x=mse_filtered, y=eps_filtered , method = "fmm", ties = mean)
  
  ### Finding the values in the experiments for the same L value
  loc2 = which(Data$L==i)
  exp_delta = Data$Delta[loc2]
  exp_kappa = Data$Kappa[loc2]
  exp_kappa_f_loc = (exp_kappa<=1 & exp_kappa>=0)
  exp_delta = exp_delta[exp_kappa_f_loc]
  exp_kappa = exp_kappa[exp_kappa_f_loc]
  
  ### Predicting the theoretical kappa values for the experimental deltas
  newMSE = exp_delta*sqrt(i)
  kappa_pred = model(x=newMSE)
  print(kappa_pred)
  res1 = kappa_pred - exp_kappa ### theory - observation
  s = length(res1)

  newDf <- data.frame(i*(c(1,rep(1,s-1))),exp_delta,exp_kappa,kappa_pred,res1)
  Df <- rbind(newDf,Df)
  
}

#### Plotting stuff
NColor<-factor(as.numeric(Df$i....c.1..rep.1..s...1...))
pl1<-ggplot(data=Df)+geom_point(aes(x=Df$exp_delta*sqrt(Df$i....c.1..rep.1..s...1...),y=Df$res1,colour=NColor))+theme_bw()+ggtitle('Residuals versus Delta*sqrt(L) for n=200 and different L')
print(pl1)

loc = which(Df$res1<=1 & Df$res1>=0)
newDf = Df[loc,]
NColor2<-factor(as.numeric(newDf$i....c.1..rep.1..s...1...))
pl2<-ggplot(data=newDf)+geom_point(aes(x=newDf$exp_delta,y=newDf$res1,colour=NColor2))+theme_bw()+ggtitle('Residuals versus Delta*sqrt(L) for n=200 and different L')
print(pl2)

