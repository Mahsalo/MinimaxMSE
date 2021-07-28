install.packages("aod")
install.packages("ggplot2")
library(aod)
library(ggplot2)
install.packages("nlstools")
library(nlstools)
install.packages("car")
library(car)
require(minpack.lm)
library(ggpubr)
theme_set(theme_pubr())
install.packages("tidyverse")
library(tidyverse)
install.packages("lmvar")
library(lmvar)
library(MASS)
install.packages("DAAG")
library(DAAG)
install.packages("dplyr")
library(dplyr)

### Loading the logit data for each (L,l) pairs and Dense Sensing Matrix
### The dimensions is changing in [100,200,400,800]

setwd("/Users/mahsa/Desktop/")

############################# USING MMSE FORMULA #####################################
LogitData <- read.csv("Dense_Logit_L1l1.txt", header=TRUE, sep=",")

head(LogitData)
attach(LogitData)
names(LogitData)

setwd("/Users/mahsa/Desktop/MinimaxEstimationOfLargestSparseMean/MinimaxEstimationOfLargestSparseMean")
eps_mse_T = read.csv("EPS_MSE_T.csv",header= TRUE)
MSE = eps_mse_T[,"MSE"]
EPS = eps_mse_T[,"EPS"]
model = splinefun(x=MSE, y=EPS , method = "fmm", ties = mean)
kappa_pred = model(x=Delta)
res = Kappa - kappa_pred
RES = abs(res)<1
R = res[RES]
D = Delta[RES]
dim = n[RES]

#####

name = c("dim","Delta","Residuals")
newData = data.frame(dim,D,R)
colnames(newData) = name

RESmat = matrix(0,nrow=length(dim),ncol=length(unique(n)))
delt = matrix(0,nrow=length(dim),ncol=length(unique(n)))
###Finding  residuals dimension-wise
i=1
for( u in unique(dim))
{
  loc1=which (newData[,"dim"]==u)
  rs=newData[loc1,"Residuals"]
  print(summary(rs))
  RESmat[1:length(rs),i]=rs
  delt[1:length(rs),i]=newData[loc1,"Delta"]
  i=i+1
}
###modeling based on resdiuals and dimensions are not separated
y = newData[,"Residuals"]
x1 = newData[,"Delta"]
x2 = newData[,"dim"]
complete_model <- nlsLM(data=newData,formula=y~(a+b*x1+c*x1^2+d/(x2)^e),start=list(a=0.1,b=0.1,c=0.1,d=0.1,e=1),control = nls.control(maxiter = 1000)) 
x3 = 1/x2
model1<-lm(data=newData,formula=y~x1+x3+x1*x3)
model2<-lm(data=newData,formula=y~x1)
model3<-lm(data=newData,formula=y~x1+x3)
c = anova(model1,model2,model3)


sum_complete <- summary(complete_model)
Resid_complete=sum(residuals(complete_model)^2)
NColor<-factor(as.numeric(newData$dim))

pl1<-ggplot(data=newData)+geom_point(aes(x=fitted(complete_model),y=resid(complete_model),colour=NColor))+theme_bw()
print(pl1)

### Fit polynomials separately
install.packages("gridExtra")
library(gridExtra)
### n=100
n1= 100
res1 = RESmat[,1]
l1 = res1>0
res1 = res1[l1]
d1 = delt[,1]
data1 = data.frame(res1,d1[l1])
poly1 <- lm(res1~poly(d1[l1],2))

p1<-ggplot(data=data1)+geom_point(aes(x=d1[l1],y=res1))+theme_bw()
p1<-p1+geom_point(aes(x=d1[l1],y=fitted(poly1)))+ggtitle('n=100')
print(p1)

### n=200
n2= 200
res2 = RESmat[,2]
l2 = res2>0
res2 = res2[l2]
d2 = delt[,2]
data2 = data.frame(res2,d2[l2])
poly2 <- lm(res2~poly(d2[l2],2))

p2<-ggplot(data=data2)+geom_point(aes(x=d2[l2],y=res2))+theme_bw()
p2<-p2+geom_point(aes(x=d2[l2],y=fitted(poly2)))+ggtitle('n=200')
print(p2)
### n=400
n3= 400
res3 = RESmat[,3]
l3 = res3>0
res3 = res3[l3]
d3 = delt[,3]
data3 = data.frame(res3,d1[l3])
poly3 <- lm(res3~poly(d3[l3],2))

p3<-ggplot(data=data3)+geom_point(aes(x=d3[l3],y=res3))+theme_bw()
p3<-p3+geom_point(aes(x=d3[l3],y=fitted(poly3)))+ggtitle('n=400')
print(p3)
### n=800
n4= 800
res4 = RESmat[,4]
l4 = res4>0
res4 = res4[l4]
d4 = delt[,4]
data4 = data.frame(res4,d4[l4])
poly4 <- lm(res4~poly(d4[l4],2))

p4<-ggplot(data=data4)+geom_point(aes(x=d4[l4],y=res4))+theme_bw()
p4<-p4+geom_point(aes(x=d4[l4],y=fitted(poly4)))+ggtitle('n=800')
print(p4)
######
grid.arrange(p1, p2,p3,p4, 
             top = "Actual residuals versus delta for different dimensions",
             nrow = 1)

#### Checking the statistical significance of two models

x3 = 1/x2
lin_model1<-lm(y~x3,data=newData)
lin_model2<-lm(y~x3+I(x3^(1.2)),data=newData)
anova(lin_model1,lin_model2)

