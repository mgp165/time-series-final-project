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


acf.2=Acf(Wt.2, 58)
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


#los mas creibles a priori
mod1 = Arima(data,order=c(1,1,1),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
mod2 = Arima(data,order=c(2,1,1),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
mod3 = Arima(data,order=c(1,1,0),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
mod4 = Arima(data,order=c(2,1,0),seasonal=list(order=c(4,1,0),period=4),lambda=-0.13)
#añadimos mas 
mod5 = Arima(data,order=c(2,1,0),seasonal=list(order=c(2,1,2),period=4),lambda=-0.13)
mod6 = Arima(data,order=c(2,1,1),seasonal=list(order=c(2,1,2),period=4),lambda=-0.13)
mod7 = Arima(data,order=c(2,1,0),seasonal=list(order=c(2,1,0),period=4),lambda=-0.13)
mod8 = Arima(data,order=c(2,1,1),seasonal=list(order=c(2,1,0),period=4),lambda=-0.13)



AIC.BIC=rbind(c(mod1$"aic",mod1$"aicc",mod1$"bic"),
              c(mod2$"aic",mod2$"aicc",mod2$"bic"),
              c(mod3$"aic",mod3$"aicc",mod3$"bic"),
              c(mod4$"aic",mod4$"aicc",mod4$"bic"),
              c(mod5$"aic",mod5$"aicc",mod5$"bic"),
              c(mod6$"aic",mod6$"aicc",mod6$"bic"),
              c(mod7$"aic",mod7$"aicc",mod7$"bic"),
              c(mod8$"aic",mod8$"aicc",mod8$"bic"))

dimnames(AIC.BIC)=list(c("Model 1","Model 2", "Model 3", 'Model 4', "Model 5","Model 6", "Model 7", 'Model 8'),c("AIC","AICc","BIC"))
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
MSE.vec5=numeric(n.data-n.min) # Mean Squared Error model 5
MAE.vec5=numeric(n.data-n.min) # Mean Absolute Error model 5
MSE.vec6=numeric(n.data-n.min) # Mean Squared Error model 6
MAE.vec6=numeric(n.data-n.min) # Mean Absolute Error model 6
MSE.vec7=numeric(n.data-n.min) # Mean Squared Error model 7
MAE.vec7=numeric(n.data-n.min) # Mean Absolute Error model 7
MSE.vec8=numeric(n.data-n.min) # Mean Squared Error model 8
MAE.vec8=numeric(n.data-n.min) # Mean Absolute Error model 8

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
  fit.mod5=Arima(data[1:(n.min+k-1)],order=c(2,1,0),seasonal=list(order=c(2,1,2),period=4),lambda=-0.13)   
  forcast.mod5 <- forecast(fit.mod5, h=1)[['mean']]
  fit.mod6=Arima(data[1:(n.min+k-1)],order=c(2,1,1),seasonal=list(order=c(2,1,2),period=4),lambda=-0.13)
  forcast.mod6 <- forecast(fit.mod6, h=1)[['mean']]
  fit.mod7=Arima(data[1:(n.min+k-1)],order=c(2,1,0),seasonal=list(order=c(2,1,0),period=4),lambda=-0.13)
  forcast.mod7 <- forecast(fit.mod7, h=1)[['mean']]
  fit.mod8=Arima(data[1:(n.min+k-1)],order=c(2,1,1),seasonal=list(order=c(2,1,0),period=4),lambda=-0.13)
  forcast.mod8 <- forecast(fit.mod8, h=1)[['mean']]
  
  MSE.vec1[k]=(data[(n.min+k)]-forcast.mod1)^2
  MAE.vec1[k]=abs(data[(n.min+k)]-forcast.mod1)
  MSE.vec2[k]=(data[(n.min+k)]-forcast.mod2)^2
  MAE.vec2[k]=abs(data[(n.min+k)]-forcast.mod2)
  MSE.vec3[k]=(data[(n.min+k)]-forcast.mod3)^2
  MAE.vec3[k]=abs(data[(n.min+k)]-forcast.mod3)
  MSE.vec4[k]=(data[(n.min+k)]-forcast.mod4)^2
  MAE.vec4[k]=abs(data[(n.min+k)]-forcast.mod4)
  MSE.vec5[k]=(data[(n.min+k)]-forcast.mod5)^2
  MAE.vec5[k]=abs(data[(n.min+k)]-forcast.mod5)
  MSE.vec6[k]=(data[(n.min+k)]-forcast.mod6)^2
  MAE.vec6[k]=abs(data[(n.min+k)]-forcast.mod6)
  MSE.vec7[k]=(data[(n.min+k)]-forcast.mod7)^2
  MAE.vec7[k]=abs(data[(n.min+k)]-forcast.mod7)
  MSE.vec8[k]=(data[(n.min+k)]-forcast.mod8)^2
  MAE.vec8[k]=abs(data[(n.min+k)]-forcast.mod8)
}

# Summary Results 

Results.mat=rbind(
  apply(cbind(MSE.vec1,MSE.vec2,MSE.vec3,MSE.vec4,MSE.vec5,MSE.vec6,MSE.vec7,MSE.vec8)[-132,],2,mean),
  apply(cbind(MAE.vec1,MAE.vec2,MAE.vec3,MAE.vec4,MAE.vec5,MAE.vec6,MAE.vec7,MAE.vec8)[-132,],2,mean)
)
dimnames(Results.mat)=list(c("MSE","MAE"),c("Model 1","Model 2", "Model 3", 'Model 4', "Model 5","Model 6", "Model 7", 'Model 8'))

Results.mat


#modelos del 1 al 6 buenos 
#nos quedamos con el 4 por AIC, BIC, AICc
#nos quedamos con 5 MSEy 6 por MAE


#esto nos lo saltamos de momento, vamos a diagnosis de residuos de estos 3
#si hay uno claramente mejor nos lo quedamos, si no, vamos a lo de lab2 de hacer realizaciones y comparar, si de nuevo
#todos nos diesen lo mismo, predecimos con los 3 y ya, puede ocurrir que tengamos 3 buenos modelos, incluso podemos 
#comentar que para predicciones reales podriamos hacer una democracia de modelos

#plot acf y pacf teoricas 
library(polynom)
acfreal = Acf(Wt.2,40)$acf
pacfreal = Pacf(Wt.2, 40)$acf


coef1 = fit.mod1$coef
ar1reg = polynomial(c(1, -coef1[1]))
ar1seas = polynomial(c(1, 0,0,0,-coef1[3],0,0,0,-coef1[4], 0,0,0,-coef1[5], 0,0,0,-coef1[6]))
ar.vec1 = ar1reg*ar1seas
acfmod1 = ARMAacf(ar=-coefficients(ar.vec1)[-1],ma=c(coef1[2]),40)
plot(acfmod1[-1],type='l',xlim=c(0,40),ylim=c(-0.5,0.5))
abline(h=0)
par(new=TRUE)
plot(acfreal[-1],type='l',col='blue',xlim=c(0,40),ylim=c(-0.5,0.5))
pacfmod1 = ARMAacf(ar=-coefficients(ar.vec1)[-1],ma=c(coef1[2]),pacf=TRUE, 40)
plot(pacfmod1,type='h',xlim=c(0,40),ylim=c(-1,1))
abline(h=0)
par(new=TRUE)
plot(pacfreal,type='h',col='blue',xlim=c(0,40),ylim=c(-1,1))

plot(acfreal[-1]-acfmod1[-1],type='l')
plot(pacfreal-pacfmod1,type='l')


coef2 = fit.mod2$coef
ar2reg = polynomial(c(1,-coef2[1],-coef2[2]))
ar2seas =  polynomial(c(1, 0,0,0,-coef2[4],0,0,0,-coef2[5], 0,0,0,-coef2[6], 0,0,0,-coef2[7]))
ar.vec2 = ar2reg*ar2seas
acfmod2 = ARMAacf(ar=-coefficients(ar.vec2)[-1],ma=c(coef2[3]),40)
plot(acfmod2[-1],type='l',xlim=c(0,40),ylim=c(-0.5,0.5))
abline(h=0)
par(new=TRUE)
plot(acfreal[-1],type='l',col='blue',xlim=c(0,40),ylim=c(-0.5,0.5))
pacfmod2 = ARMAacf(ar=-coefficients(ar.vec2)[-1],ma=c(coef2[3]),pacf=TRUE, 40)
plot(pacfmod2[-1],type='l',xlim=c(0,40),ylim=c(-1,1))
abline(h=0)
par(new=TRUE)
plot(pacfreal,type='l',col='blue',xlim=c(0,40),ylim=c(-1,1))

plot(acfreal-acfmod2[-1],type='l')
plot(pacfreal[-1]-pacfmod2[-1],type='l')


coef3 = fit.mod3$coef
ar3reg = polynomial(c(1,-coef3[1]))
ar3seas =  polynomial(c(1, 0,0,0,-coef3[2],0,0,0,-coef3[3], 0,0,0,-coef3[4], 0,0,0,-coef3[5]))
ar.vec3 = ar3reg*ar3seas
acfmod3 = ARMAacf(ar=-coefficients(ar.vec3)[-1],ma=0,40)
pacfmod3 = ARMAacf(ar=-coefficients(ar.vec3)[-1],ma=0,pacf=TRUE,40)

plot(acfreal[-1]-acfmod3[-1],type='l')
plot(pacfreal[-1]-pacfmod3[-1],type='l')


coef4 = fit.mod4$coef
ar4reg = polynomial(c(1,-coef4[1],-coef4[2]))
ar4seas =  polynomial(c(1, 0,0,0,-coef4[3],0,0,0,-coef4[4], 0,0,0,-coef4[5], 0,0,0,-coef4[6]))
ar.vec4 = ar4reg*ar4seas
acfmod4 = ARMAacf(ar=-coefficients(ar.vec4)[-1],ma=0,40)
pacfmod4 = ARMAacf(ar=-coefficients(ar.vec4)[-1],ma=0,pacf=TRUE,40)

plot(acfreal[-1]-acfmod4[-1],type='l')
plot(pacfreal[-1]-pacfmod4[-1],type='l')


coef5 = fit.mod5$coef
ar5reg = polynomial(c(1,-coef5[1],-coef5[2]))
ar5seas =  polynomial(c(1, 0,0,0,-coef5[3],0,0,0,-coef5[4]))
ar.vec5 = ar5reg*ar5seas
acfmod5 = ARMAacf(ar=-coefficients(ar.vec5)[-1],ma=c(0,0,0,coef5[5],0,0,0,coef5[6]),40)
pacfmod5 = ARMAacf(ar=-coefficients(ar.vec5)[-1],ma=c(0,0,0,coef5[5],0,0,0,coef5[6]),pacf=TRUE,40)

plot(acfreal[-1]-acfmod5[-1],type='l')
plot(pacfreal[-1]-pacfmod5[-1],type='l')


coef6 = fit.mod6$coef
ar6reg = polynomial(c(1,-coef6[1],-coef6[2]))
ar6seas =  polynomial(c(1, 0,0,0,-coef6[4],0,0,0,-coef6[5]))
ar.vec6 = ar6reg*ar6seas
ma6reg = polynomial(c(1,coef6[3]))
ma6seas = polynomial(c(1, 0,0,0,coef6[6],0,0,0,coef6[7]))
ma.vec6 = ma6reg*ma6seas
acfmod6 = ARMAacf(ar=-coefficients(ar.vec6)[-1],ma=coefficients(ma.vec6)[-1],40)
pacfmod6 = ARMAacf(ar=-coefficients(ar.vec6)[-1],ma=coefficients(ma.vec6)[-1],pacf=TRUE,40)

plot(acfreal[-1]-acfmod6[-1],type='l')
plot(pacfreal-pacfmod6[-1],type='l')



diferencias_acfs_concatenadas = c(acfreal[-1]-acfmod1[-1],acfreal[-1]-acfmod2[-1],acfreal[-1]-acfmod3[-1],acfreal[-1]-acfmod4[-1],acfreal[-1]-acfmod5[-1],acfreal[-1]-acfmod6[-1])
ggseasonplot(ts(diferencias_acfs_concatenadas,frequency=40))
diferencias_pacfs_concatenadas = c(pacfreal-pacfmod1,pacfreal-pacfmod2,pacfreal-pacfmod3,pacfreal-pacfmod4,pacfreal-pacfmod5,pacfreal-pacfmod6)
ggseasonplot(ts(diferencias_pacfs_concatenadas,frequency=40))

#####

#diagnosis of residuals

require(astsa)
source('Diagnostic.R')
require(portes)
require(nortest)
#MODELO 4

S.ACF(mod4$residuals)
NonParametric.Tests(mod4$residuals)
Check.normality(mod4$residuals)
My.Ljung.Box(mod4$residuals, 7)


#MODELO 5

S.ACF(mod5$residuals)
NonParametric.Tests(mod5$residuals)
Check.normality(mod5$residuals)
My.Ljung.Box(mod5$residuals,7)

#MODELO 6

S.ACF(mod6$residuals)
NonParametric.Tests(mod6$residuals)
Check.normality(mod6$residuals)
My.Ljung.Box(mod6$residuals,8)

#con el S.ACF obtenemos lo mismo para los 3 
#con el check normality lo mismo
#con el nonparametric tenemos un poco de evidencia para descartar mod6
#con el ljung box nos quedamos con el 4 (k=8 por lo del paper linkeado en el lab6)


##fitting and predictions


plot(seq(1,232), data,ylab="registrations",xlab="",type="l",main="titutlo",lwd=2)
points(seq(1,232),forecast(mod4)$fitted,col="red",type="l",lwd=2)

n.data=length(data)
fit.mod4=Arima(data[1:224],order=c(2,1,0),seasonal=list(order=c(4,1,0),period=4 ),lambda=-0.13)   
forcast.mod4 <- forecast(fit.mod4, h=8) #predicciones

plot(seq(180,232), data[180:232],ylab="Registrations",xlab="",type="l",main="titulo ",ylim=c(1.4,3.5))
points(seq(225,232),forcast.mod4$mean,type="l",col="red") # predictions
points(seq(225,232),forcast.mod4$lower[,2],type="l",col="blue") # lower CI
points(seq(225,232),forcast.mod4$upper[,2],type="l",col="blue") # up CI 
#en vez de h=8 meterle 12 

Predic.mod4=forecast(mod4,32)
Predic.mod4

plot(Predic.mod4,200)
