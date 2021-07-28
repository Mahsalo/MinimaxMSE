# 

k <- 3;  # how many frames in maximum
nMonte <- 500;
Xk <- rep(NaN,nMonte);
for( iMonte in 1:nMonte){
	z <- rnorm(k);
	xk <- max(z);
	Xk[iMonte] <- xk;
}
mu.k <- mean(Xk);  # expected largest k-sample of a standard normal
var.k <- var(Xk);  # variance largest k-sample of a standard normal

tiny = 1.3-4;
nGrid=100;
epsGrid <- seq(from=1/(nGrid+1),to=nGrid/(2*nGrid+1),length=nGrid);
tab <- matrix(0,nGrid,6)
colnames(tab) <- c("eps","T.eps","mse.eps","T2.eps","mse2.eps","mse3.eps")
TGrid <- seq(from=0,to=mu.k,length=nGrid);
for( iGrid in 1:nGrid){
	eps <- epsGrid[iGrid];
	a = (1-2*eps);
	b = -2*(1-eps)*mu.k;
	c =  (1-eps)*(mu.k^2+var.k) - eps;
	mse.epsf <- function(eps,T,mu=mu.k,v=var.k) {(1-eps)*((mu - T)^2 + v) + eps  * (T^2 +1)};
	#T.eps <- ifelse(abs(a) > tiny,(-b + sqrt(b^2-4*a*c))/(2*a), (mu.k^2 + var.k) /(4*mu.k));
	#mse.eps <- (1-eps)*((mu.k - T.eps)^2 + var.k) + eps  * (T.eps^2 +1);
	mse0 <- mu.k^2 +var.k; T0 = 0;
	for(j in 1:nGrid){ T1 <- TGrid[j]; mse1<- mse.epsf(eps,T1); if(mse1 <  mse0){mse0 = mse1; T0 <- T1}}
	T2 = (1-eps)*mu.k;
	mse2 <- (1-eps)*(eps^2 + 1-eps)*mu.k^2 + (1-eps)*var.k + eps;
	mse3 <- mse.epsf(eps,T2)
	tab[iGrid,] <- c(eps,T0,mse0,T2,mse2,mse3)
}
par(mfrow=c(1,2))
par(oma=c(2,2,2,2))
plot(tab[,"eps"],tab[,"mse.eps"],ylab="mse",xlab="eps")
mtext(text=sprintf("minimax MSE and threshold for largest mean of k, k=%i",k),side=3,line=0,outer=TRUE)
lines(tab[,"eps"],tab[,"mse2.eps"],col='green')
lines(tab[,"eps"],tab[,"mse3.eps"],col='red')
plot(tab[,"eps"],tab[,"T.eps"],ylab="thresh",xlab="eps")
lines(tab[,"eps"],tab[,"T2.eps"],col='green')
dev.print(postscript,file=sprintf("minimax-largest-sparse-mean-of-k%i.eps",k))
