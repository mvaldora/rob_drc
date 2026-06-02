rm(list = ls())
library(robustbase)
library(RobStatTM)
library(tidyverse)
#library(statmod)
library(drc)
source("robustDRC_functions.R")
source("derivatives.R")
source("generate.R")

#Data<-"Hep2"
Data<-"a549"

if(Data=="Hep2"){
DataMix <- read.csv2("Data/GF_CF_Hep2.csv")
DataSub1 <- read.csv2("Data/GF_Hep2.csv")
DataSub2 <- read.csv2("Data/CF_Hep2.csv")
}
if(Data=="a549"){
DataMix <- read.csv2("Data/GF_CF_A549.csv")
DataSub1 <- read.csv2("Data/GF_a549.csv")
DataSub2 <- read.csv2("Data/CF_a549.csv")
}

names(DataMix)
names(DataSub1)
names(DataSub2)

## change column names  of 
## Dat1, Dat2, Datmix1 and Datmix2, if necessary.


Dat1<-list()
if(Data=="Hep2"){
Dat1$x <-  DataSub1$CONCMARCH
Dat1$y <- DataSub1$MTT
}else if(Data=="a549"){
  Dat1$x <- DataSub1$CONC.MARCH
  Dat1$y <- DataSub1$MTT
}else if(Data == "Hep2Blanco"){
  Dat1$x <-  DataSub1$CONCMIPA
  Dat1$y <- DataSub1$MTT
}
est1all <- estimateDRM(Dat1)
concentra <- as.numeric(names(table(Dat1$x)))

summary(drm(Dat1$y ~ Dat1$x, fct= LL.4()))
  Dat2<-list()
Dat2$x<-DataSub2$CONC
Dat2$y<-DataSub2$MTT
est2all<-estimateDRM(Dat2)
concentra<-as.numeric(names(table(Dat2$x)))


Datmix1<-list()
Datmix1$x<-DataMix$CONC.MARCH
Datmix1$y<-DataMix$MTT
estmix1all<-estimateDRM(Datmix1)

Datmix2<-list()
Datmix2$x<-DataMix$CONC.CIPER
Datmix2$y<-DataMix$MTT
estmix2all<-estimateDRM(Datmix2)


est1<-est1all$coefficients
est2<-est2all$coefficients
estmix1<-estmix1all$coefficients
estmix2<-estmix2all$coefficients

est_ritz_weighted<-function(Dat){
concentra <- unique(Dat$x)
ll <- length(concentra)
s <- m<- 0
for (i in 1:ll) {
  s[i] <- sd((Dat$y[Dat$x == concentra[i]]))
  m[i] <- mad((Dat$y[Dat$x == concentra[i]]))
}
sds <- madi <- 0
for (i in 1:ll) {
  sds <- sds + s[i] * (Dat$x == concentra[i])
  madi <- madi + m[i] * (Dat$x == concentra[i])
}
drm(Dat$y~ Dat$x,fct=LL.4(),weights = 1 / sds)
}

alpha<-0.05/2
CIs_ritz_con_pesos<-function(Dat){
ajuste_ritz<-est_ritz_weighted(Dat)  
std_errors<-summary(ajuste_ritz)[[3]][,2]
CIs_ritz <-cbind(ajuste_ritz$coefficients,
                        ajuste_ritz$coefficients-qnorm((1-alpha/2))*std_errors,
                        ajuste_ritz$coefficients+qnorm((1-alpha/2))*std_errors)
colnames(CIs_ritz) <- c("estimate", "lower", "upper")
CIs_ritz[c(3,1,4,2),]
}


est1clas<-est_ritz_weighted(Dat1)$coefficients[c(3,1,4,2)]
est2clas<-est_ritz_weighted(Dat2)$coefficients[c(3,1,4,2)]
estmix1clas<-est_ritz_weighted(Datmix1)$coefficients[c(3,1,4,2)]
estmix2clas<-est_ritz_weighted(Datmix2)$coefficients[c(3,1,4,2)]

est1classp<-drm(Dat1$y~ Dat1$x,fct=LL.4())$coefficients[c(3,1,4,2)]
est2classp<-drm(Dat2$y~ Dat2$x,fct=LL.4())$coefficients[c(3,1,4,2)]
estmix1classp<-drm(Datmix1$y~ Datmix2$x,fct=LL.4())$coefficients[c(3,1,4,2)]
estmix2classp<-drm(Datmix1$y~ Datmix2$x,fct=LL.4())$coefficients[c(3,1,4,2)]



plot(Dat1$x,Dat1$y,xlab="Dose",ylab="Response")
grilla<-seq(min(Dat1$x),max(Dat1$x))
lines(grilla,freg3(grilla,est1clas),col=4)
lines(grilla,freg3(grilla,est1),col=3)
lines(grilla,freg3(grilla,est1classp),col=5)


CIs_ritz_mix1<-CIs_ritz_con_pesos(Datmix1)
CIs_ritz_mix2<-CIs_ritz_con_pesos(Datmix2)



## change zero doses to very small ones

doses <- unique(Dat1$x)
zerodoses<-Dat1$x == 0
Dat1$x[zerodoses] <-doses[2]/10 
intsdat1<-confintDRC(Dat1,alpha=alpha)

zerodoses<-Dat2$x == 0
doses <- unique(Dat2$x)
Dat2$x[zerodoses] <-doses[2]/10 
intsdat2<-confintDRC(Dat2,alpha=alpha)
zerodoses<-Datmix1$x == 0
doses <- unique(Datmix1$x)
Datmix1$x[zerodoses] <-doses[2]/10 
intsdatmix1<-confintDRC(Datmix1,alpha=alpha)
zerodoses<-Datmix2$x == 0
doses <- unique(Datmix2$x)
Datmix2$x[zerodoses] <-doses[2]/10 
Datmix2$x<-round(Datmix2$x,10)
intsdatmix2<-confintDRC(Datmix2,alpha=alpha)

par(mfrow=c(1,2))
i<-5
ff<-freg3
edsrob1 <- (EDS(Dat1, est1clas, i / 10, ff))
edsrob2 <- (EDS(Dat2, est2clas, i / 10, ff))
edsrob3 <- (EDS(Datmix1, estmix1clas, i / 10, ff))
edsrob4 <- (EDS(Datmix2, estmix2clas, i / 10, ff))
plot(
  c(0, edsrob1),
  c(edsrob2, 0),
  type = "l",
  xlab = "march",
  ylab = "ciper",
  xlim = c(0, max(edsrob1, edsrob3) ),
  ylim = c(0, max(edsrob2, edsrob4) ),
  main = paste("ED", 50, "LS"),
  xaxt="n",yaxt="n")
points(edsrob3, edsrob4,pch=19)


polygon(c(CIs_ritz_mix1[3,2],CIs_ritz_mix1[3,3],CIs_ritz_mix1[3,3],CIs_ritz_mix1[3,2]),
        c(CIs_ritz_mix2[3,2],CIs_ritz_mix2[3,2],CIs_ritz_mix2[3,3],CIs_ritz_mix2[3,3]))



edsrob1 <- (EDS(Dat1, est1, i / 10, ff))
edsrob2 <- (EDS(Dat2, est2, i / 10, ff))
edsrob3 <- (EDS(Datmix1, estmix1, i / 10, ff))
edsrob4 <- (EDS(Datmix2, estmix2, i / 10, ff))
    plot(
      c(0, edsrob1),
      c(edsrob2, 0),
      type = "l",
      xlab = "march",
      ylab = "ciper",
      xlim = c(0, max(edsrob1, edsrob3) ),
      ylim = c(0, max(edsrob2, edsrob4) ),
      main = paste("ED", 50, "ROB"),xaxt="n",yaxt="n")
points(edsrob3, edsrob4,pch=19)

polygon(c(intsdatmix1[3,3],intsdatmix1[3,4],intsdatmix1[3,4],intsdatmix1[3,3]),
        c(intsdatmix2[3,3],intsdatmix2[3,3],intsdatmix2[3,4],intsdatmix2[3,4]))


