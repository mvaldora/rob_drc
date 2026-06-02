
######################################
#' computes the function P(y,theta), that is the loss function evaluated in one observation of the form (x1,x2,...xr,y), that is observation 1 for all doses.
#'
#' @param theta A vector of parameters for the four-parameter-logistic function
#' @param theta A vector of doses
#' @param y A vector of responses for all doses
#' @param sigmas A vector of the same length as x containing the standard deviations of y for each doses
#' @param cc The calibrating constant of the rho function
#' @return A real number
#' @export

Rho_all_doses <- function(theta, x, y, sigmass, cc){
  k<-length(unique(x))
  aa <- (y - freg3(x, theta)) / sigmass
  sum(rho(aa, cc))
}



#' computes derivative of the bisquare rho-function i.e. the bisquare psi function *6/c^2
#' @param u A real number
#' @return A real number
#' @export


rhoprime <- function(x, cc) {
  (6 * psi(x, cc) / cc ^ 2) * (abs(x) < cc)
}



#' computes the second derivative of the bisquare rho-function
#'
#' @param u A real number
#' @return A 4X4 matrix containing the second derivatives of the 4-paameter logistic regression function
#' @export

rhoprimeprime <- function(x, cc) {
  6 * (1 - (x / cc) ^ 2) * (1 - 5 * (x / cc) ^ 2) * (abs(x / cc) < 1) / cc ^
    2
}


#' computes derivative of the bisquare psi-function
#'
#' @param u A real number
#' @return A 4X4 matrix containing the second derivatives of the 4-paameter logistic regression function
#' @export

psiprime <- function(x, cc) {
  u <- x / cc
  (1 - u ^ 2) * (1 - 5 * u ^ 2) * (abs(u) < 1)
}

#' computes the gradient of the 4-parameter logistic regression function with respect to tita
#'
#' @param z A vector of doses of size n.
#' @param tita A vector of parameters
#' @return A nX4 matrix containing the gradient of the 4-parameter logistic regression function evaluated in each doese x
#' @export


dfdtita <- function(z, tita) {
  zord <- sort(z)
  if (zord[1] == 0) {
    hh <- (zord[1] + zord[2]) / 2
    z[z == 0] <- hh
  }
  tita1 <- tita[1]
  tita2 <- tita[2]
  tita3 <- tita[3]
  tita4 <- tita[4]
  denmenos1 <- exp(tita2 * (log(z) - log(tita3)))
  #denmenos1 <- (z/tita3)^tita2
  dtita1 <- 1 / (denmenos1 + 1)
  termdtita2y3 <- (tita1 - tita4) * (dtita1) ^ 2 * denmenos1
  dtita2 <- -termdtita2y3 * log(z / tita3)
  dtita3 <- termdtita2y3 * tita2 / tita3
  dtita4 <- 1 - dtita1
  cbind(dtita1, dtita2, dtita3, dtita4)
}

#' computes the score function of the robust estimator for dose response models
#'
#' @param x A vector of doses of size k.
#' @param y A vector of responses of size k.
#' @param sigmas A vector whose i-th entry is the standard deviation of the responses for the i-th dose.
#' @param tita A vector of parameters
#' @return A vector of length equal to 4, the length of tita.
#' @export

Psi_all_doses <- function(x, y, sigmass, tita, cc) {
  k <- length(unique(x))
  if (length(y) != length(x))
    print("error")
  ans <- matrix(NA, nrow = length(x), ncol = 4)
  rhoprimeev <-rhoprime((y - freg3(x, tita)) / sigmass, cc)
  dfdtitaev <- dfdtita(x, tita)
  ans[, 1] <- -rhoprimeev / sigmass * dfdtitaev[, 1]
  ans[, 2] <- -rhoprimeev / sigmass * dfdtitaev[, 2]
  ans[, 3] <- -rhoprimeev / sigmass * dfdtitaev[, 3]
  ans[, 4] <- -rhoprimeev / sigmass * dfdtitaev[, 4]
  colSums(ans)
}


#' computes the matrix of second derivatives of the 4-parameter logistic regression function with respect to tita
#' @param z A real number indicating a dose
#' @param tita A vector of parameters
#' @return A 4X4 matrix containing the second derivatives of the 4-paameter logistic regression function
#' @export

df2dtita2 <- function(z, tita) {
  tita1 <- tita[1]
  tita2 <- tita[2]
  tita3 <- tita[3]
  tita4 <- tita[4]
  denmenos1 <- exp(tita2 * (log(z) - log(tita3)))
  dtita1 <- 1 / (denmenos1 + 1)
  deriv_seg <- matrix(NA, 4, 4)
  deriv_seg[1, 1] <- deriv_seg[4, 4] <-
    deriv_seg[1, 4] <- deriv_seg[4, 1] <- 0
  deriv_seg[1, 2] <-
    deriv_seg[2, 1] <- -dtita1 ^ 2 * denmenos1 * log(z /
                                                       tita3)
  deriv_seg[1, 3] <-
    deriv_seg[3, 1] <- dtita1 ^ 2 * denmenos1 * tita2 / tita3
  deriv_seg[3, 3] <-
    (tita1 - tita4) * dtita1 * denmenos1 * tita2 / tita3 * (2 * deriv_seg[1, 3] - dtita1 * (tita2 +
                                                                                              1) / tita3)
  deriv_seg[2, 2] <-
    -(tita1 - tita4) * log(z / tita3) * denmenos1 * dtita1 * (2 * deriv_seg[1, 2] + dtita1 * log(z / tita3))
  deriv_seg[2, 3] <-
    deriv_seg[3, 2] <-
    -(tita1 - tita4) * denmenos1 * dtita1 * (2 * deriv_seg[1, 3] * log(z / tita3) -
                                               dtita1 *  log(z / tita3) * tita2 /
                                               tita3 -
                                               dtita1 / tita3)
  deriv_seg[2, 4] <-
    deriv_seg[4, 2] <- denmenos1 * dtita1 ^ 2 * log(z / tita3)
  deriv_seg[3, 4] <-
    deriv_seg[4, 3] <- -denmenos1 * dtita1 ^ 2 * tita2 / tita3
  deriv_seg
}


#' computes the jacobian matrix of the score function of the robust estimators for dose response models psi_all_doses
#'
#' @param z A real number indicating the dose
#' @param tita A vector of parameters
#' @return A 4X4 matrix containing the second derivatives of the 4-paameter logistic regression function
#' @export

Jpsi_all_doses <- function(x, y, sigmass, tita, cc) {
  k <- length(unique(x))
  ans <- matrix(NA, nrow = k, ncol = 4)
  rhoprimeev<- rhoprime((y - freg3(x, tita)) / sigmass, cc)
  rhoprimeprimeev <- rhoprimeprime((y - freg3(x, tita)) / sigmass, cc)
  dfdtitaev_fun <- function(t) {
    dfdtita(t, tita)
  }
  df2dtita2ev_fun <- function(t) {
    df2dtita2(t, tita)
  }
  ans <- 0
  for (i in 1:length(x)) {
    ans <- ans +
      rhoprimeprimeev[i] * matrix(dfdtitaev_fun(x[i])) %*% dfdtitaev_fun(x[i])/sigmass[i]^2 -
      rhoprimeev[i] * df2dtita2ev_fun(x[i])/sigmass[i]
  }
  ans
}



  #' computes the asymptotic covariance matrix of of the robust estimator for dose response models 
  #'
  #' @param z A real number indicating the dose
  #' @param tita A vector of parameters
  #' @return A 4X4 matrix containing 
  
  varest3<-function(dat, cc =3.44, tita, sigmass){
    concentra<- unique(dat$x)
    k<-length(concentra)
    termfact1 <- termfact2 <- 0
    lista_rtas_por_dosis<-list()
    n<-0
    for(i in 1:k){
      lista_rtas_por_dosis[[i]]<-dat$y[dat$x==concentra[i]]
      n[i]<-length(dat$y[dat$x==concentra[i]])
    }
    nmin<-min(n)
      for(i in 1:nmin){
      dat1 <- list()
      dat1$x <- unique(dat$x) 
      dat1$y <-0
      for(j in 1:k){dat1$y[j] <- lista_rtas_por_dosis[[j]][i]}
      sigmas <- rep(0,k)
      for( j in 1:k){
        sigmas[j]<-sigmass[which(dat1$x==dat1$x[j])]
      }
      termfact1 <- termfact1 + Jpsi_all_doses(dat1$x, dat1$y, sigmas, tita, cc)
      termfact2 <- termfact2 + Psi_all_doses(dat1$x, dat1$y, sigmas, tita, cc) %*% 
        t(Psi_all_doses(dat1$x, dat1$y, sigmas, tita, cc))
    }
    solve(termfact1/i)%*%(termfact2/i)%*%t(solve(termfact1/i))
  }
  

  
#' computes confidence intervals for the 4-paramter-logistic dose response models 
#'
#' @param z A real number indicating the dose
#' @param tita A vector of parameters
#' @return A 4X4 matrix containing 
#' @export

  
  confintDRC <- function(dat, tt = id, cc = 3.44, alpha=0.05){
    ajuste_robDRC<- estimateDRM(dat, tt = id, cc = 3.44)
    titahat<-ajuste_robDRC$coefficients
    sigmashat<-ajuste_robDRC$dev
    n <- length(dat$x)/length(unique(dat$x))
    asvar <- varest3(dat, cc, titahat, sigmashat)/n
    intervals <- 
    c(titahat, sqrt(diag(asvar)), titahat-qnorm((1-alpha/2))*sqrt(diag(asvar)), titahat+qnorm((1-alpha/2))*sqrt(diag(asvar)))
    matrix(intervals,nrow=4)
    } 
  
  
  
  
  