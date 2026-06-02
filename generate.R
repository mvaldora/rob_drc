## --------------------------------------------------
## Generates data following a dose response model
##

generateDRM_C0<-function(doses,theta,sigmas,ns){
  x<-y<-c()
  ndosis<-length(doses)
  n<-ns[1]
  for(i in 1:n){
    y<-c(y,freg3(doses,theta)+rnorm(rep(1,ndosis),rep(0,ndosis),sigmas))
  }
  ans<-list()
  ans$x<-rep(doses, n)
  ans$y<-y
  ans
}


## para cada n, 10 % de las observaciones se cambian por un valor y = 1.5 IQR+ 3ercuartil
contaminate<-function(Dat, cont, epsilon){
  doses <- unique(Dat$x)
  ns <- length(Dat$x)/length(doses)
  nsouts <- round(ns*epsilon)
  if(nsouts>0){
  for(i in 1:length(doses)){
    yaux <- Dat$y[Dat$x==doses[i]]
    desvioestiqr <- IQR(yaux)
    contnueva <- quantile(yaux,0.75)+ cont * desvioestiqr
    yaux[1:nsouts] <- contnueva + rnorm(nsouts,0,desvioestiqr)
    Dat$y[Dat$x==doses[i]] <- yaux
  }
  }
  ans<-list()
  ans$x<-Dat$x
  ans$y<-Dat$y
  ans
}

if(FALSE){
  theta<-c(0.32,   4.24, 112.78,   0.01 )
  doses<-c(0,  48,  96, 120, 240, 360, 480, 960)
  sigmas<-c(0.091028621, 0.088337152, 0.075790617, 0.083957533, 0.037103477,
                     0.005550276, 0.005237775, 0.005616060)
  ns<-c(144, 144, 144, 144, 144, 144, 144, 144) 
  Dat <- generateDRM_C0(doses,theta,sigmas,ns)
  Dat2 <- contaminate(Dat, 5, 0.1)
  plot(Dat2$x,Dat2$y,col=2)  
  points(Dat$x,Dat$y)
}