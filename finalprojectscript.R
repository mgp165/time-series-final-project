#calculations
#Total quarterly production of cement in Australia (in millions of tonnes), first quarter 1958 - fourth quarter 2015 (n=232).
library(forecast)
data = read.table("tutoria15681581.txt",header=FALSE,sep=' ') 
data=data[,2]  

z = ts(data, start=c(1958,1), frequency=4)

#plot the data

plot(data, type='l', axes=FALSE, xlab="", ylab='cement(in millions of tons)', main='Total quarterly production of cement in Australia, first quarter 1958 - fourth quarter 2015')
box()
axis(1, at=seq(0,232,by=24), labels=seq(1958,2015,by=6))
axis(2, at=seq(0,3,by=0.5),labels=seq(0,3,by=0.5))  

#comentar que vemos trend y seasonal y varianza distinta

ggseasonplot(z, main='cemento ')

#la periodicidad es efectivamente 4

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
#summary(lm(log(apply(matrix(z,ncol=n.tilde,byrow=T),1,mean)))~log(sqrt(apply(matrix(z,ncol=n.tilde,byrow=T),1,var))))

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

plot(Wt.2,type='l')
ggseasonplot(ts(Wt.2,frequency=4))
#con ggseasonal ver estacionariedad


acf.2=Acf(Wt.2, 100)
plot(acf.2$acf[-1],type='h')
points(seq(4,100,by=4), acf.2$acf[seq(5,101,by=4)],type="h",col="blue",lwd=1.5)
#mostraremos con abs para ver mejor las caracteristicas
plot(abs(acf.2$acf[-1]),type='h')
points(seq(4,100,by=4), abs(acf.2$acf[seq(5,101,by=4)]),type="h",col="blue",lwd=1.5)

#si queremos ver solamente la ACF de la seasonal part sin valor absoluto...
abline(h=0)
plot(seq(4,100,by=4),acf.2$acf[seq(5,101,by=4)],type="h",col="blue",lwd=1.5)

#vemos sinusoide en seasonal, en regular vemos decreasing geomertically
#por tanto, la seasonal puede ser arma o ar, la regular puede ser arma o ar

#pacf igual

pacf.2=Pacf(Wt.2, 100)


plot(pacf.2$acf,type='h')
points(seq(4,100,by=4), pacf.2$acf[seq(4,101,by=4)],type="h",col="blue",lwd=1.5)

plot(abs(pacf.2$acf),type='h')
points(seq(4,100,by=4), abs(pacf.2$acf[seq(4,101,by=4)]),type="h",col="blue",lwd=1.5)

#parte seasonal: 4 coefs no cero,resto si ==> AR4

#seleccion modelos en la hoja de peñas, foto por whatsapp



mod1 = Arima(data,order=c(1,1,1),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
mod2 = Arima(data,order=c(2,1,1),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
mod3 = Arima(data,order=c(1,1,0),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
mod4 = Arima(data,order=c(2,1,0),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)


AIC.BIC=rbind(c(mod1$"aic",mod1$"aicc",mod1$"bic"),
              c(mod2$"aic",mod2$"aicc",mod2$"bic"),
              c(mod3$"aic",mod3$"aicc",mod3$"bic"),
              c(mod4$"aic",mod4$"aicc",mod4$"bic"))
dimnames(AIC.BIC)=list(c("Model 1","Model 2", "Model 3", 'Model 4'),c("AIC","AICc","BIC"))
AIC.BIC


n.min=100  # minimal sample size  required to estimate a model
n.data=length(data)
MSE.vec1=numeric(n.data-n.min) # Mean Squared Error model 1
MAE.vec1=numeric(n.data-n.min) # Mean Absolute Error model 1
MSE.vec2=numeric(n.data-n.min) # Mean Squared Error model 2
MAE.vec2=numeric(n.data-n.min) # Mean Absolute Error model 2
MSE.vec3=numeric(n.data-n.min) # Mean Squared Error model 3
MAE.vec3=numeric(n.data-n.min) # Mean Absolute Error model 3
MSE.vec4=numeric(n.data-n.min) # Mean Squared Error model 4
MAE.vec4=numeric(n.data-n.min) # Mean Absolute Error model 4

for (k in 1:(n.data-n.min))
{
  fit.mod1=Arima(data[1:(n.min+k-1)],order=c(1,1,1),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)   
  forcast.mod1 <- forecast(fit.mod1, h=1)[['mean']]
  fit.mod2=Arima(data[1:(n.min+k-1)],order=c(2,1,1),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
  forcast.mod2 <- forecast(fit.mod2, h=1)[['mean']]
  fit.mod3=Arima(data[1:(n.min+k-1)],order=c(1,1,0),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
  forcast.mod3 <- forecast(fit.mod3, h=1)[['mean']]
  fit.mod4=Arima(data[1:(n.min+k-1)],order=c(2,1,0),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
  forcast.mod4 <- forecast(fit.mod4, h=1)[['mean']]
  
  MSE.vec1[k]=(data[(n.min+k)]-forcast.mod1)^2
  MAE.vec1[k]=abs(data[(n.min+k)]-forcast.mod1)
  MSE.vec2[k]=(data[(n.min+k)]-forcast.mod2)^2
  MAE.vec2[k]=abs(data[(n.min+k)]-forcast.mod2)
  MSE.vec3[k]=(data[(n.min+k)]-forcast.mod3)^2
  MAE.vec3[k]=abs(data[(n.min+k)]-forcast.mod3)
  MSE.vec4[k]=(data[(n.min+k)]-forcast.mod4)^2
  MAE.vec4[k]=abs(data[(n.min+k)]-forcast.mod4)
}

# Summary Results 
Results.mat=rbind(
  apply(cbind(MSE.vec1,MSE.vec2,MSE.vec3,MSE.vec4)[-280,],2,mean),
  apply(cbind(MAE.vec1,MAE.vec2,MAE.vec3,MAE.vec4)[-280,],2,mean)
)
dimnames(Results.mat)=list(c("MSE","MAE"),c("Model 1","Model 2", "Model 3",'Model 4'))

Results.mat

