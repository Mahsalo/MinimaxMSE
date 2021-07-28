# 

k.list <- c(1,2,5,7,10);
tab.list <- list();
for(k in k.list){
	nMonte <- 500;
	Xk <- rep(NaN,nMonte);
	for( iMonte in 1:nMonte){
		z <- rnorm(k);
		xk <- max(z);
		Xk[iMonte] <- xk;
	}
	mu.k <- mean(Xk);  # expected largest k-sample of a standard normal
	var.k <- var(Xk);  # variance largest k-sample of a standard normal

    if(k==1){
    	imse.0 <- function(T0,K0=1){
		    mse.0 <-function(x,T0,K0) {(x-T0)^2*(x>T0)*dnorm(x)}
		    integrate(function(x){mse.0(x,T0=T0,K0=1)},T0,max(T0+1,4))$value 
	    }
    } else{
		imse.0 <- function(T0,K0){
			mse.0 <-function(x,T0,K0) {(x-T0)^2*(x>T0)*K0*pnorm(x)^(K0-1)*dnorm(x)}
			integrate(function(x){mse.0(x,T0=T0,K0=K0)},T0,max(T0+1,4))$value 
	    }
    }

tiny = 1.3-4;
nGrid=100;
epsGrid <- seq(from=1/(nGrid+1),to=nGrid/(2*nGrid+1),length=nGrid);
tab <- matrix(0,nGrid,3)
colnames(tab) <- c("eps","T.eps","mse.eps")
TGrid <- seq(from=0,to=3*(mu.k+sqrt(var.k)),length=3*nGrid);
for( iGrid in 1:nGrid){
	eps <- epsGrid[iGrid];
	mse.epsf <- function(eps,T0,K=k) {(1-eps)*imse.0(T0,K0=K) + eps*(T0^2 +1)};
	mse0 <- mu.k^2 +var.k; T0 = 0;
	for(j in 1:length(TGrid)){ T1 <- TGrid[j]; mse1<- mse.epsf(eps,T1); if(mse1 <  mse0){mse0 = mse1; T0 <- T1}}
	tab[iGrid,] <- c(eps,T0,mse0)
}
tab.list[[sprintf("k%i",k)]] = tab;
}
par(mfrow=c(1,2))
par(oma=c(2,2,2,2))
colrs <- c("red","blue","green","cyan","magenta","gold","black","aqua")
for(i.k in 1:length(k.list)){
	k <- k.list[i.k];
tab = tab.list[[sprintf("k%i",k)]];
clr <- colrs[i.k];
if(k == 1){
	plot(tab[,"eps"],tab[,"mse.eps"],ylab="mse",xlab="eps",col=clr,ylim=c(0,1))
    mtext(text=sprintf("minimax MSE and threshold for largest mean of k"),side=3,line=0,outer=TRUE)
} else {
	lines(tab[,"eps"],tab[,"mse.eps"],col=colrs[i.k])	
}
}
legend(.25,.25,names(tab.list),col=colrs[1:length(tab.list)],pch=c("o",".",".",".","."),lwd=2)
for(i.k in 1:length(k.list)){
	k <- k.list[i.k];
tab = tab.list[[sprintf("k%i",k)]];
clr <- colrs[i.k];
if(k == 1){
	plot(tab[,"eps"],tab[,"T.eps"],ylab="thresh",xlab="eps",col=clr)
} else {
	lines(tab[,"eps"],tab[,"T.eps"],col=colrs[i.k])	
}
}
dev.print(postscript,file=sprintf("minimax-largest-sparse-mean-of-k-varying-k.eps"))
