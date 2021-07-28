
Y <- read.csv("Data_all_logits.txt",header=TRUE)
X <- readRDS(file="MinimaxMSECurves.vs.K.RDS")

Y.k <- unique(Y$L)
X.k <-unique(X[,"k"])
K <- intersect(Y.k,X.k)

par(mfrow=c(2,2))
for(k in K){
	l.min = min(Y$l[Y$L==k])
	Y.delta = Y$Delta[Y$L==k & Y$l==l.min]
	Y.kappa = Y$Kappa[Y$L==k & Y$l==l.min]
	X.k = X[,"k"]
	X.eps   = X[X.k==k,"eps"]
	X.mse   = X[X.k==k,"mse"]
	plot (Y.delta,Y.kappa,col="red")
	lines(X.mse  ,X.eps  ,col="green")
	print(c(k,l.min))
}
#
par(mfrow=c(1,2))
first=TRUE
for(k in K){
	l.min = min(Y$l[Y$L==k])
	Y.delta = Y$Delta[Y$L==k & Y$l==l.min]
	Y.kappa = Y$Kappa[Y$L==k & Y$l==l.min]
	X.k = X[,"k"]
	X.eps   = X[X.k==k,"eps"]
	X.mse   = X[X.k==k,"mse"]
	if(first) { plot (Y.delta,Y.kappa,col="red");first=FALSE} else{ lines (Y.delta,Y.kappa,col="red")} 
	print(c(k,l.min))
}
first=TRUE
for(k in K){
	l.min = min(Y$l[Y$L==k])
	Y.delta = Y$Delta[Y$L==k & Y$l==l.min]
	Y.kappa = Y$Kappa[Y$L==k & Y$l==l.min]
	X.k = X[,"k"]
	X.eps   = X[X.k==k,"eps"]
	X.mse   = X[X.k==k,"mse"]
	if(first) { plot (X.mse  ,X.eps  ,col="green"); first=FALSE} else{ lines (X.mse  ,X.eps  ,col="green")} 
	print(c(k,l.min))
}
par(mfrow=c(1,1))
first=TRUE
prox = function(v,d){ which(abs(d-v) == min(abs(d-v))) }
colrs <- c("red","blue","green","cyan","magenta","gold","black","aqua")
for(i.k in 1:length(K)){
	k <- K[i.k];
	clr <- colrs[i.k];
	l.min = min(Y$l[Y$L==k])
	Y.delta = Y$Delta[Y$L==k & Y$l==l.min]
	Y.kappa = Y$Kappa[Y$L==k & Y$l==l.min]
	X.k = X[,"k"]
	X.eps   = X[X.k==k,"eps"]
	X.mse   = X[X.k==k,"mse"]
	Y.mse   = Y.delta;
	for(i in 1:length(Y.kappa)){
		j = prox(Y.kappa[i],X.eps);
		Y.mse[i] = X.mse[j]
	}
	subst = Y.delta < 0.4;
	if(first) { plot (Y.delta[subst],Y.mse[subst]/sqrt(k) ,col=clr,main="agreement of transformed MSE with PT",ylab="MSE/sqrt(L)",xlab="Empirical PT"); abline(0,1,col="black"); first=FALSE} else{ lines (Y.delta[subst] ,Y.mse[subst]/sqrt(k)  ,col=clr)} 
	print(summary(lm(Y.mse/sqrt(k) ~ Y.delta,subset=subst)))
}
dev.print(postscript,file="Agreement.ps")



