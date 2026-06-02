funw1 <- function(u)
{
  if (!missing(u))
    min(1 / (u ^ 2), 3 - 3 * u ^ 2 + u ^ 4)
  else
    NA
}


funw0 <- function(u)
{
  funw1(u / 1.564) / (1.564 ^ 2)
}

#' computes the 4-parameter logistic regression function
#'
#' @param z A real number indicating the dose
#' @param tita A 4-dimensional vector of parameters
#' @return A real number
#' @export
freg3 <- function(z, tita)
{
  y <- tita[4] + (tita[1] - tita[4]) / (1 + (z / tita[3]) ^ tita[2])
  return(y)
}
#' computes the 4-parameter logistic regression function
#'
#' @param z A real number indicating the dose
#' @param tita1 Parameter indicating the limit of \code{funcparas} when x tends to 0
#' @param tita2 Parameter indicating the slope of the logistic function
#' @param tita3 Parameter corresponding to EC50
#' @param tita4 Parameter indicating the limit of \code{funcparas} when x tends to infinity
#' @return A real number
#' @export
funparas <- function(z, tita1, tita2, tita3, tita4)
{
  y = tita4 + (tita1 - tita4) / (1 + (z / tita3) ^ tita2)
  return(y)
}
######################################
#
# GM-estimador
#
######################################
#' computes the loss function of the robust estimator for heteroscedastic 4-parameter logistic model
#'
#' @param theta A 4-dimensional vector of paramters
#' @param x1 A vector of doses of dimension n
#' @param y A vector of responses of dimension n
#' @param ww The standard deviations of the responses for each dose
#' @param sigma.est A preliminary estimator of sigma, the standard deviation of the residuals
#' @param ce A calibrating constant 
#' @return A function
#' @export

gmest <- function(theta, x1, y, ww, sigma.est, ce) {
  vv = rep(0, length(y))
  if (theta[3] <= 0.00000001) {
    pp = theta[4] * rep(1, length(y))
  } else{
    pp = theta[4] + (theta[1] - theta[4]) /
      (1 + (x1 / theta[3]) ^ theta[2])
  }
  aa <- (y - pp) / (ww * sigma.est*ce )
  vv<-rho(aa,1)
#  for (i in 1:length(y)) {
#    vv[i] <- rho(aa[i])
#  }
  f <- sum(vv)
  f
}
######################################
#' computes Tukey's bisquare rho-function
#'
#' @param u A real number
#' @return A real number
#' @export
rho <- function(x,cc)
{
  Mchi(x, cc, "bisquare")
}

#' computes Tukey's bisquare psi-function
#'
#' @param u A real number
#' @return A real number
#' @export

psi <- function(x,cc){
  Mpsi(x, cc = cc,"bisquare")
}

###############################
#' estimates the paramaters in the 4-parameter logistic regression function
#'
#' @param Dat A list with two elements: x is the vector of doses and y the vector of responses
#' @param tt A function to transform the responses
#' @param cc A calibrating constant
#' @return A real number
#' @export
estimateDRM <- function(Dat, tt = id, cc = 3.44, ccpsi= 3.14, ccscale =1.5) {
  concentra <- unique(Dat$x)
  ll <- length(concentra)
  s <- m<- 0
  for (i in 1:ll) {
    s[i] <- sd(tt(Dat$y[Dat$x == concentra[i]]))
    m[i] <- mad(tt(Dat$y[Dat$x == concentra[i]]))
  }
  sds <- madi <- 0
  for (i in 1:ll) {
    sds <- sds + s[i] * (Dat$x == concentra[i])
    madi <- madi + m[i] * (Dat$x == concentra[i])
  }
curva <-
  drm(
    tt(Dat$y) ~ Dat$x,
    data = Dat,
    weights = 1 / sds,
    fct = LL.4()
  )
  b <- curva$fit$par[1]
  c <- curva$fit$par[2]
  d <- curva$fit$par[3]
  e <- curva$fit$par[4]
  tita <- c(d, b, e, c)
  lower = c(
    tita1 = 0.001,
    tita2 = b - 5 * abs(b),
    tita3 = 0.0001,
    tita4 = c - 5 * abs(c)
  )
  upper = c(
    tita1 = d + 5 * abs(d),
    tita2 = b + 5 * abs(b),
    tita3 = e + 5 * abs(e),
    tita4 = c + 5 * abs(c)
  )
  salida.s <-
    robustbase:::nlrob.MM(
      y ~ funparas(x, tita1, tita2, tita3, tita4),
      data = Dat,
      ctrl = nlrob.control("MM", psi = "bisquare", init="S", fnscale = NULL,
                           tuning.chi.scale =  ccscale,
                           tuning.psi.M     = ccpsi),
      lower = lower,
      upper = upper
    )
  s.estimador <- salida.s$init$par
  theta0inis <- s.estimador
  res.inis <- (Dat$y - freg3(Dat$x, theta0inis))/madi 
  sigma.mad.inis <- mad(res.inis)
  nobs <- length(Dat$x)
  DL <- 1.0
  a <- rep(0, nobs)
  sigma0 <- sigma.mad.inis
  while (DL > 0.001) {
    for (r in 1:nobs) {
      a[r] <- funw0(res.inis[r] / sigma0) * (res.inis[r] ^ 2)
    }
    sigma <- sqrt((2 / nobs) * sum(a, na.rm = T))
    DL <- min(abs(sigma - sigma0) / sigma0, abs(sigma - sigma0))
    sigma0 <- sigma
  }
  sigma.m.inis  <- sigma
  ce0gmhr = cc
  salida.GM.w1.m <-
    nlminb(
      start = theta0inis,
      obj = gmest,
      x1 = Dat$x,
      y = Dat$y,
      ww = madi,
      sigma.est = sigma.m.inis,
      ce = ce0gmhr,
      lower = c(0, -Inf, 0.0001, 0)
    )
  salida.GM.w1.m$par
  ans <- list()
  ans$classic <- tita
  ans$robustinitial <- s.estimador
  ans$coefficients <- salida.GM.w1.m$par
  ans$dev <- madi
  ans$residuals <-
    (res <- (Dat$y - freg3(Dat$x, salida.GM.w1.m$par))) / madi
  ans$rweights <- robustbase::Mpsi(ans$residuals,cc=cc,psi="bisquare")/
    (Dat$y - freg3(Dat$x, salida.GM.w1.m$par))
  ans
}


if(FALSE){
  library(drc)
  theta<-c(0.32,   4.24, 112.78,   0.01 )
  doses<-c(0,  48,  96, 120, 240, 360, 480, 960)
  sigmas<-c(0.091028621, 0.088337152, 0.075790617, 0.083957533, 0.037103477,
            0.005550276, 0.005237775, 0.005616060)
  ns<-c(144, 144, 144, 144, 144, 144, 144, 144) 
  dat <- generateDRM_C0(doses,theta,sigmas,ns)
  system.time(estim_all1 <- estimateDRM(dat))
  estim_all1$coefficients
  estim_all1$robustinitial
  estim_all1$classic
}

    

#' Computes EDp for any p between 0 and 1
#'
#' @param Dat A list with two elements: x is the vector of doses and y the vector of responses
#' @param est Estimators of the parameters of the regression function
#' @param per The desired p, default is 0.5
#' @param ff The regression function, default is the four parameter regression function
#' @return A real number
#' @export
EDS <- function(Dat, est, per = 0.5, ff = freg3) {
  values <- ff(Dat$x, est)
  fmax <- max(values)
  fmin <- min(values)
  yper <- (fmax - fmin) * per + fmin
  faux <- function(x) {
    ff(x, est) - yper
  }
  ans <- uniroot(faux, c(min(Dat$x), max(Dat$x)))
  ans$root
}


#' Identity function
#'
#' @param x A real number
#' @return A real number
#' @export
id <- function(x)
  x

#' Computes concentration indices and isoboles for p=0.1,0.2,...0.9
#'
#' @param Dat1 A list with two elements: x is the vector of doses of substance 1 and y the vector of responses
#' @param Dat2 A list with two elements: x is the vector of doses of substance 2 and y the vector of responses
#' @param Datmix1 A list with two elements: x is the vector of doses of substance 1 in the mixture of substance 1 and 2and y the vector of responses
#' @param Datmix2 A list with two elements: x is the vector of doses of substance 2 in the mixture of substance 1 and 2and y the vector of responses
#' @param est1 Estimators of the parameters of the regression function for \code{Dat1}
#' @param est2 Estimators of the parameters of the regression function for \code{Dat2}
#' @param estmix1 Estimators of the parameters of the regression function for \code{Datmix1}
#' @param estmix2 Estimators of the parameters of the regression function for \code{Datmix2}
#' @param plot "ED" or "CI". If "ED" isobole plots are produced, if "CI", cooncentration index plots are produced.
#' @param ff The regression function, default is the four parameter regression function
#' @return Concentration indices and concentration index plot (if plot="CI") or isoboles (if plot="ED") for p=0.1,0.2,...0.9
#' @export
interaction_analysis <-
  function(Dat1,
           Dat2,
           Datmix1,
           Datmix2,
           est1,
           est2,
           estmix1,
           estmix2,
           plot="ED",
           ff = freg3) {
    edsrob1 <- edsrob2 <- edsrob3 <- edsrob4 <- NULL
    for (i in 1:9) {
      edsrob1[i] <- (EDS(Dat1, est1, i / 10, ff))
      edsrob2[i] <- (EDS(Dat2, est2, i / 10, ff))
      edsrob3[i] <- (EDS(Datmix1, estmix1, i / 10, ff))
      edsrob4[i] <- (EDS(Datmix2, estmix2, i / 10, ff))
    }
    if(plot=="ED"){
    par(mfrow = c(3, 3))
    for (i in 9:1) {
      plot(
        c(0, edsrob1[i]),
        c(edsrob2[i], 0),
        type = "l",
        xlab = "march",
        ylab = "ciper",
        xlim = c(0, max(edsrob1[i], edsrob3[i]) + 1),
        ylim = c(0, max(edsrob2[i], edsrob4[i]) + 1),
        main = paste("ED", i / 10)
      )
      points(edsrob3[i], edsrob4[i])
    }
   par(mfrow = c(1, 1))
   }
   ps<-seq(0.1,0.9,0.1)
   CI<- edsrob3/edsrob1+ edsrob4/edsrob2
   if(plot=="CI"){
   plot(ps,CI,xlab="Effect Level",ylab="Combination Index")
   lines(ps,CI)
   abline(1,0)
   }
   ans<-data.frame(as.matrix(CI))
   colnames(ans)<-"CI"
   rownames(ans)<-seq(0.1,0.9,0.1)
   ans
  }
