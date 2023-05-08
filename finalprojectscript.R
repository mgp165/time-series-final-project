#calculations
#Total quarterly production of cement in Australia (in millions of tonnes), first quarter 1958 - fourth quarter 2016 (n=232).
library(forecast)
data = read.table("tutoria15681581.txt",header=FALSE,sep=' ') 
data=data[,2]  

z = ts(data, start=c(1958,1), frequency=4)

#plot the data

plot(data, type='l', axes=FALSE, xlab="", ylab='cement(in millions of tons)', main='Total quarterly production of cement in Australia (in millions of tonnes), first quarter 1958 - fourth quarter 2015')
box()
axis(1, at=seq(0,232,by=24), labels=seq(1958,2015,by=6))
axis(2, at=seq(0,3,by=0.5),labels=seq(0,3,by=0.5))  

#comentar que vemos trend y seasonal y varianza distinta

ggseasonalplot(z, main='cemento ')


#box cox transformation

# x: observed time series
# n.tilde: number of observations per group

n.tilde = 4



Plot.var=function(x,n.tilde)  
{
  A.temp=matrix(x,ncol=n.tilde,byrow=T)     # Each row correspond to a group of n.tilde consecutive observations 
  xbar.vec=apply(A.temp,1,mean)       # vector of means 
  s.vec=sqrt(apply(A.temp,1,var))     # vector of standard deviations 
  plot(log(xbar.vec),log(s.vec),xlab="log(mean)",ylab="log(standard deviation)",main=paste("Number of observations per group:",n.tilde),type="p")
}


#ver plot y comprobar que la relacion es lineal para hacer la transformacion

Plot.var(z, n.tilde)


# x: observed time series
# n.tilde: number of observations per group


BoxCox=function(x,n.tilde)
{
  A.temp=matrix(x,ncol=n.tilde,byrow=T)    # Each row correspond to a group of n.tilde consecutive observations 
  xbar.vec=apply(A.temp,1,mean)      # vector of means 
  s.vec=sqrt(apply(A.temp,1,var))    # vector of standard deviations 
  lambda=1-lm(log(s.vec)~log(xbar.vec))$coefficients[2]
  yt=(x^lambda-1)/lambda
  lambda1=round(lambda,2)
  ts.plot(yt,xlab="",ylab="",main=paste("Transformed series obtained using the Box-Cox transformation with lambda=",lambda1),type="l")
}
# x: observed time series
# n.tilde: number of observations per group

BoxCox.out=function(x,n.tilde)
{
  A.temp=matrix(x,ncol=n.tilde,byrow=T)   
  xbar.vec=apply(A.temp,1,mean)       
  s.vec=sqrt(apply(A.temp,1,var))
  lambda=1-lm(log(s.vec)~log(xbar.vec))$coefficients[2]
  yt=(x^lambda-1)/lambda
  #return(list(yt,lambda))
  return(yt)
}

BoxCox(data, n.tilde)
X.tilde = BoxCox.out(data, n.tilde)

#now the variance is constant as we can see

Acf(X.tilde,main="ACF of the logarithm of the registration series",xlab="Lag",ylab="ACF")

#decreases linearly slowly 
# we apply differencing operator

Wt.1 = diff(X.tilde, lag=1, differences=1)
ts.plot(Wt.1)
#the mean of Wt.1 is constant and zero so the process is stationary

#as Wt.1 is stationary, X.tilde is integrated of order 1


#plot ggseasonal para comprobar que es estacionaria


Acf(Wt.1, main="ACF of the first difference of the logarithm of the registration series",xlab="Lag",ylab="ACF")

#as we see acf non zero on seasonal lags, we apply seasonal differencing operator (s=4)

Wt.2 = diff(Wt.1, lag=4, differences = 1)

ts.plot(Wt.2)

acf.2=Acf(Wt.2, 100)

plot(acf.2)
points(seq(4,56,by=4), acf.2[seq(4,56,by=4)],type="h",col="blue",lwd=1.5)
points(seq(1:3),type="h",col="red",lwd=1.5)

acf.2[seq(4,60,by=4)]

prueba=c(0.364,0.219,0.088,0.068,0.145,0.033,0.164,0.132,0.053,0.016,0.05,0.075,0.04,0.002,0.009)
plot(prueba,type='h')
