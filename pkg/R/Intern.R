grossA <- function (x, y, sigma, n) {
### this can never be infinity!
  return(n / ( (n-1) * sigma^2 * (1 - x -y)^2) ) 
}

psi <- function(x, k, betapara){
    qbeta(k- pbeta(x, betapara, betapara), betapara, betapara)
  }

integrandOC.1 <- function(t, n, deltaL, deltaU, x) {
  wt <- sqrt(t*(n-1))
  dchisq(t,n-1) * ( pnorm ( -deltaU + (2*x -1) * wt ) -
                   pnorm ( -deltaL + wt ) )
}

integralOC.1 <- function (x, mu, sigma, n){
  dU <-  sqrt(n)* ( mu-1) / sigma # deltaU( mu, sigma, n) 
  dL <- sqrt(n) *( mu+1 ) / sigma # deltaL( mu, sigma, n) 
  integrate(
            integrandOC.1,
            0, grossA(x, 0,sigma, n),
            n=n, deltaL=dL, deltaU=dU, x=x, subdivisions=2500)$value
}


OC.1 <- function(mu, sigma,  n, k) {
 betap <- (n-2)/2
 integralOC.1(qbeta(k, betap, betap ), mu, sigma,  n)
}


OC.2 <- function(mu, sigma,  n, k, resolution=2501){
### ToDo: rewrite using integrate.
  betap <- (n-2)/2
  xgitter <- seq(0,qbeta(k, betap, betap), length.out=resolution )[-resolution]
  resolution <- resolution-1
  werte <- rep(NA, resolution)
  gridwidth <- (xgitter[resolution] - xgitter[1])/resolution 
  for (i in 1:resolution) {werte[i] <-  innerOC.2(xgitter [i], mu, sigma, n, k)  }
  sum(werte)*gridwidth
}


innerOC.2 <- function(x, mu, sigma, n, k) {
  dL <- sqrt(n) *( mu+1 ) / sigma # deltaL( mu, sigma, n) 
  result <- integrate(ininOC.2,0,grossA( psi(x, k, (n-2)/2 ), x, sigma,  n ),subdivisions=2500, x=x, deltaL=dL,  n=n)$value
  result
}

ininOC.2 <- function(t, x, deltaL, n) {
  wt <- sqrt(t*(n-1))
  dchisq(t,n-1)* wt * dnorm(-deltaL + (1-2*x) * wt)
}






mup <- function(mu,plevel,sigma){ pnorm((-1-mu)/sigma) + pnorm((mu-1)/sigma) - plevel }
### diese Funktion ist nur eine Hilfsfunktion für die numerische Nullstellensuche
### über den Parameter mu gegeben plevel und sigma. Funktioniert.

sigmap0 <- function(p){ -1/qnorm(p/2) }
### Dies ist das maximal interessante sigma zu einem gegebenen p-level,
### korrespondiert zu einem mu=0

integrandOCML <- function(t, n, mu, sigma, k){
  t[1] <- min( 1e-6, (t[1]+t[2]) /2)
  musigma <- sigma*sqrt(t/(n-1))
   
  lokalmu <- rep(NA, length(musigma))
  for ( i in 1:length(musigma)) {
    lokalmu[i] <- uniroot(mup, interval=c(1e-6, 1 ) , sigma=musigma[i], plevel=k )$root
  }
  dchisq(t,df=n-1) * (pnorm( sqrt(n)*(lokalmu - mu )/sigma )-
             pnorm( sqrt(n)*(-lokalmu -mu )/sigma))        
}

iOCML <- function(t, n, mu, sigma, k){
  musigma <- sigma*sqrt(t/(n-1))
  lokalmu <- rep(NA, length(musigma))
  for ( i in 1:length(musigma)) {
    lokalmu[i] <- uniroot(mup, interval=c(1e-6, 1 ) , sigma=musigma[i], plevel=k )$root
  }
  gn <- dchisq(t,df=n-1)
  gn * pnorm( sqrt(n)*(lokalmu - mu )/sigma ) -  gn * pnorm( sqrt(n)*(-lokalmu -mu )/sigma )
}
