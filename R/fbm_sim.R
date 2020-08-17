## ---------------------------------------------------------------------------------------------------------------
##
## fbm_sim.R
##
## This script contains the primary function for simulatiing discrete paths of Fractional Brownian Motion.
## Included are several different simulation methods (Hosking, Cholesky, Davies-Harte) to provide the user
## with different options for their specific needs.
##
##
## Author: Justin Yu, M.S. Financial Engineering, Stevens Institute of Technology
##
## Date Created: 2020-08-13
##
## Email: 732jhy@gmail.com
##
## ---------------------------------------------------------------------------------------------------------------
##
## Notes:
##        This code was inspired by the following paper:
##
##        Dieker, Ton. Simulation of fractional Brownian motion. M.Sc. theses,
##        University of Twente, Amsterdam, The Netherlands (2004)
##
##        This paper can be found at the following web address:
##        http://www.columbia.edu/~ad3217/fbm/thesis.pdf
## ---------------------------------------------------------------------------------------------------------------


# Fractional Brownian Motion (Exact) Methods:
#
# 1. Davies-Harte Method --- Complexity: O(nlog(n))
# 2. Hosking Method --- Complexity: O(N^2)
# 3. Cholesky Method --- Complexity: O(N^3)


#' @title Fractional Brownian Motion
#'
#' @description Returns a discrete simulated realization of fractional Brownian motion of size N.
#'
#' @param N Number of simulated points
#' @param H Hurst parameter
#' @param t Length of the simulated series
#' @param method The simulation method to be used
#'
#' @return A single discrete path of fractional Brownian motion by simulating fractional Gaussian noise (see \code{\link{FGN}}).
#' @export
#'
#' @examples
#' \dontrun{
#' FBM(1000, 0.25, 1.0, 'davies-harte')
#' }
#' \dontrun{
#' FBM(1000, 0.25, 1.0, 'hosking')
#' }
#' \dontrun{
#' FBM(1000, 0.25, 1.0, 'cholesky')
#' }
FBM <- function(N, H, t=1, method='davies-harte') {

  # Hurst parameter check
  if((H < 0) | (H > 1)) {
    stop("Hurst parameter must be between 0 and 1.")
  }

  # autocovariance function
  gamma <- function(k, h) {
    return(0.5*( (abs(k-1)^(2*h)) - (2*(abs(k)^(2*h))) + (abs(k+1)^(2*H))))
  }

  # Davies Harte method:
  if(method == 'davies-harte') {

    # Get eigenvalues
    g = c()
    for(k in 0:(N-1)){g <- c(g, gamma(k,H))}
    r = c(g, 0, rev(g)[1:(N-1)])
    j = seq(0, ((2*N)-1))
    K = (2*N)-1
    i = complex(real=0, imaginary=1)
    lk = rev(fft(r*exp(2*pi*i*K*j*(1/(2*N)))))

    # Generate random variables
    Vj <- cbind(seq(0,0,length.out=2*N), seq(0,0,length.out=2*N))
    Vj[1,1] <- rnorm(1)
    Vj[N+1, 1] <- rnorm(1)
    Vj1 <- rnorm(N-1)
    Vj2 <- rnorm(N-1)
    Vj[2:N,1] <- Vj1
    Vj[2:N,2] <- Vj2
    Vj[(N+2):(2*N),1] <- rev(Vj1)
    Vj[(N+2):(2*N),2] <- rev(Vj2)

    # Compute Z (fractional Gaussian Noise)
    wk = seq(0,0,length.out=2*N)
    wk[1] <- sqrt(lk[1]/(2*N))*Vj[1,1]
    wk[2:N] <- sqrt(lk[2:N]/(4*N))*(Vj[2:N,1] + i*Vj[2:N,2])
    wk[N+1] <- sqrt(lk[N+1]/(2*N))*Vj[N+1,1]
    wk[(N+2):(2*N)] <- sqrt(lk[(N+2):(2*N)]/(4*N))*(Vj[(N+2):(2*N),1] - i*Vj[(N+2):(2*N),2])
    Z = fft(wk)
    fGn = Z[1:N]
    fBm = cumsum(fGn)*(N^(-H))
    return(Re((t^H)*fBm))


  } else if(method == 'hosking') {

    # Starting values for the recursion
    X = c(rnorm(1)) # vector containing fractional Gaussian Noise
    mu = c(gamma(1,H)*X[1])
    sigsq = c(1-(gamma(1,H)^2))
    tau = c(gamma(1,H)^2)
    d = c(gamma(1,H))

    for(n in 2:N){
      Fmat = apply(diag(n),2,rev)
      cn = c()
      for(k in 0:(n-1)){cn <- c(cn, gamma(k+1,H))}

      # Generate sig(n+1)^2
      s = sigsq[n-1] - ((gamma(n,H) - tau[n-1])^2)/sigsq[n-1]
      sigsq = c(sigsq, s)

      # Generate d(n+1)
      phi = (gamma(n,H) - tau[n-1])/sigsq[n-1]
      d = d - phi*rev(d)
      d = c(d, phi)

      # mu(n+1) and tau(n+1)
      Xnplus1 = mu[n-1] + sigsq[n-1]*rnorm(1)
      X = c(X, Xnplus1)
      mu = c(mu, d %*% rev(X))
      tau = c(tau, cn %*% Fmat %*% d)
    }

    fBm = cumsum(X)*(N^(-H))
    return((t^H)*fBm)


  } else if(method == 'cholesky') {

    # Starting values for recursion
    L = matrix(0, N,N)
    X = seq(0,0,length.out=N)
    V = rnorm(N)

    L[1,1] = 1.0
    X[1] = V[1]

    L[2,1] = gamma(1,H)
    L[2,2] = sqrt(1 - (L[2,1]^2))
    X[2] = sum(L[2,1:2])

    # Generate Cholesky matrix by row
    for(i in 3:N) {

      # First row value
      L[i,1] = gamma(i-1,H)

      # Middle values
      for(j in 2:(i-1)){
        L[i,j] = (1/L[j,j])*(gamma(i-j,H) - (L[i,1:(j-1)]%*%(L[j,1:(j-1)])))
      }

      # Last value
      L[i,i] = sqrt(1 - sum(L[i,1:(i-1)]^2))

      # Use row to compute X(n+1)
      X[i] = L[i,1:i] %*% V[1:i]
    }

    fBm = cumsum(X)*(N^(-H))
    return((t^H)*fBm)

  } else {stop("Incorrect method")}
}



# Testing each method
#
# # Davies-Harte method
# h02 = FBM(200,0.2,1,method='davies-harte')
# plot(seq(1,200),h02, type='l')
# head(h02)
#
# h05 = FBM(200,0.5,1,method='davies-harte')
# plot(seq(1,200),h05, type='l')
# head(h05)
#
# h08 = FBM(200,0.8,1,method='davies-harte')
# plot(seq(1,200),h08, type='l')
# head(h08)
#
#
# # Hosking's method
# h02 = FBM(200,0.2,1,method='hosking')
# plot(seq(1,200),h02, type='l')
# head(h02)
#
# h05 = FBM(200,0.5,1,method='hosking')
# plot(seq(1,200),h05, type='l')
# head(h05)
#
# h08 = FBM(200,0.8,1,method='hosking')
# plot(seq(1,200),h08, type='l')
# head(h08)
#
#
# # Cholesky's method
# h02 = FBM(200,0.2,1,method='cholesky')
# plot(seq(1,200),h02, type='l')
# head(h02)
#
# h05 = FBM(200,0.5,1,method='cholesky')
# plot(seq(1,200),h05, type='l')
# head(h05)
#
# h08 = FBM(200,0.8,1,method='cholesky')
# plot(seq(1,200),h08, type='l')
# head(h08)














