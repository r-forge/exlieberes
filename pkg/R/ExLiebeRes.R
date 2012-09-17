### Following a complete rewrite with L=-1, U=1 fixed (speedup and hope for numerical improvements.)
## DEBUG <- FALSE
## This variable is used for debugging the package. Disabled for publication.

### mup <- function(mu,plevel,sigma){pnorm((-1-mu)/sigma) + pnorm((mu-1)/sigma) - plevel}
### sigmap0 <- function(p){-1/qnorm(p/2) } 

grossA <- function (x, y, sigma, n) {
### if (DEBUG)  cat(x, " ", y, "  ", sigma, " ", n, "\n")
  if ( x + y < 1 ) return(n / ( (n-1) * sigma^2 * (1 - x -y)^2) ) else
  return(Inf) 
}

psi <- function(x, k, betapara){
    qbeta(k- pbeta(x, betapara, betapara), betapara, betapara)
  }


OC <- function(mu, sigma,n, k) {
  part1 <- OC.1(mu, sigma,  n, k)
                                        #  cat("Part1 ", part1,"\n")
  part2 <- OC.2(mu, sigma, n, k )
                                        #  cat("Part2 ", part2,"\n")
  part1 + 2*part2
}

deltaU <- function(mu, sigma, n) { sqrt(n)* ( mu-1) / sigma}
deltaL <- function(mu, sigma, n) { sqrt(n) *( mu+1 ) / sigma}

integrandOC.1 <- function(t, n, deltaL, deltaU, x) {
  wt <- sqrt(t*(n-1))
  dchisq(t,n-1) * ( pnorm ( -deltaU + (2*x -1) * wt ) -
                   pnorm ( -deltaL + wt ) )
}

integralOC.1 <- function (x, mu, sigma, n){
  dU <- deltaU( mu, sigma, n)
  dL <- deltaL( mu, sigma, n)
  integrate(
            integrandOC.1,
            0, grossA(x,0,sigma, n),
            n=n, deltaL=dL, deltaU=dU, x=x, subdivisions=2500)$value
}


OC.1 <- function(mu, sigma,  n, k) {
 betap <- (n-2)/2
 integralOC.1(qbeta(k, betap, betap ), mu, sigma,  n)
}


OC.2 <- function(mu, sigma,  n, k, resolution=2501){
  betap <- (n-2)/2
  xgitter <- seq(0,qbeta(k, betap, betap), length.out=resolution )[-resolution]
  resolution <- resolution-1
#  cat("Int Obergrenze B^{-1}", xgitter[resolution], "\n")
  werte <- rep(NA, resolution)
  gridwidth <- (xgitter[resolution] - xgitter[1])/resolution 
  for (i in 1:resolution) {werte[i] <-  innerOC.2(xgitter [i], mu, sigma, n, k)  }
 # plot(werte)
  sum(werte)*gridwidth
}


innerOC.2 <- function(x, mu, sigma, n, k) {
  dL <- deltaL( mu, sigma, n)
### if (DEBUG) cat("x ", x, "mu ", mu, "sigma ", sigma, "n ", n, "k ", k, "\n")  
  result <- integrate(ininOC.2,0,grossA( psi(x, k, (n-2)/2 ), x, sigma,  n ),subdivisions=2500, x=x, deltaL=dL,  n=n)$value
###  cat("Inneres Integral = ", result, "\n")
  result
}

ininOC.2 <- function(t, x, deltaL, n) {
  wt <- sqrt(t*(n-1))
  dchisq(t,n-1)* wt * dnorm(-deltaL + (1-2*x) * wt)
}

LRP <- function(p1, p2, alpha, beta)
  {
    ### Es wird das minimale n zurückgegeben, für das eine Lösung exisitert.
    ### Näherungsplan
#cat(p1, p2, alpha, beta,l,"\n")
    n1 <- ceiling(((qnorm(1-alpha) - qnorm(beta)) / (qnorm(p2)-qnorm(p1)) )^2 )
#cat("n1 ", n1, "\n")
    l1n <- qnorm(1-alpha) + sqrt(n1)*qnorm(p1)
    l2n <- qnorm(beta) + sqrt(n1)*qnorm(p2)
    l1 <- (l1n+l2n)/2
#   cat("l1 ", l1, "\n")
    
    n0 <- floor( n1*(1 + l1^2/(2*n1)) )
    l <- l1 * sqrt(1 + l1^2/(2*n1))
    

   cat("Näherungsplan n0=",n0, " l0=", l, "\n" )
    return(list(ntilde=n0-1, ltilde=l,ktilde=pbeta(0.5 + l/(2*(n0-1)), (n0-2)/2, (n0-2)/2) ))
### man braucht nicht den exakten Näherungsplan!
    n<- n0 - 2
    
    while ( n < n0+3){
      cat("n= ", n, "\n")
      cat("(i)  ", pt(l,n-1,sqrt(n)*qnorm(p1)), ">= 1-alpha", 1-alpha, "\n" )
      cat("(ii) ", pt(l,n-1,sqrt(n)*qnorm(p2)), "<  beta", beta , "\n")
      if ( pt(l,n-1,sqrt(n)*qnorm(p1)) >= (1-alpha) & pt(l,n-1,sqrt(n)*qnorm(p2)) < beta){
        cat("min n =", n,"\n")
        plan <- list(ntilde=n, ltilde=l, ktilde=pbeta(0.5 + l/(2*(n-1)), (n-2)/2, (n-2)/2))
        return(plan)
      }
      n <- n+1
    }
  }

### z.B: L(n = 34 l=-11.41 p1=0.01 p2=0.06 alpha, beta = 0.1) ist ein Plan

## OC.Isop <- function( p, n, k, resolution=100)
##   {
## ###    cat("OC.Isop( p=", p, " n=", n, " k=", k, ")\n")
##     sigma.max <- sigmap0(p)
##     resolution <- 100
##     sigmas<- seq(1e-6, sigma.max*(1-1e-5), length.out=resolution)
##     mus <- rep(NA, resolution)
##     integralwert<- mus
    
##     for (i in 3:resolution){
##       mus[i] <- uniroot(mup, interval=c(0,1e6), tol=1e-15, sigma=sigmas[i], plevel=p)$root
##     }
##     for (i in 3:resolution){
##        integralwert[i] <- OC(mus [i], sigmas[i], n, k )

##     }
##     return(integralwert)
##   }

OC.Isop <- function( p, n, k, resolution=100)
  {
###    cat("OC.Isop( p=", p, " n=", n, " k=", k, ")\n")
    sigma.max <- sigmap0(p)
    resolution <- 100
    sigmas<- seq(1e-6, sigma.max*(1-1e-5), length.out=resolution)
    mus <- rep(NA, resolution)
    integralwert<- mus
    
    for (i in 1:resolution){
      mus[i] <- uniroot(mup, interval=c(0,4*sigma-max), tol=1e-15, sigma=sigmas[i], plevel=p)$root
    }
    for (i in 1:resolution){
       integralwert[i] <- OC(mus [i], sigmas[i], n, k )

    }
    return(integralwert)
  }


OC.Isop.quick <- function( p, n, k)
  {
###    cat("OC.Isop( p=", p, " n=", n, " k=", k, ")\n")
    sigma.max <- sigmap0(p)*(1-1e-5)
##    resolution <- 100
##    sigmas<- seq(1e-6, sigma.max*(1-1e-5), length.out=resolution)
##    mus <- rep(NA, resolution)
##   integralwert<- mus

    mus <- uniroot(mup, interval=c(0,10), tol=1e-15, sigma=sigma.max, plevel=p)$root
    
##    for (i in 3:resolution){
##      mus[i] <- uniroot(mup, interval=c(0,1e6), tol=1e-15, sigma=sigmas[i], plevel=p)$root
##    }
##    for (i in 3:resolution){
##       integralwert[i] <- OC(mus [i], sigmas[i], n, k )
##
##    }
    return(OC(mus, sigma.max,n,k))
  }


#### ok, das liefe durch wie oben.
#### Automatisieren
#### Zu \alpha \beta p_1 p2 finde optimalen Plan

fOptPlan <- function( p1, p2, alpha, beta, given.plan=NULL)
  {
    if (is.null(given.plan))
      {
        startplan <- LRP(p1, p2, alpha, beta)
        nmin <- startplan$ntilde
        kopt <- startplan$ktilde
      }
    else
      {
        nmin <- given.plan$ntilde
        kopt <- given.plan$ktilde
      }
    kdiff <- 2e-5
    PlanGefunden <- FALSE
    direction <- "unknown"
    
    while(! PlanGefunden) {
      cat("Teste Plan mit p1=" , p1, " p2= ", p2, " n= ", nmin, " und k=", kopt, "\n")
      result1 <- OC.Isop(p1, nmin, kopt)
      plot(result1[result1 > 0.1])
      Lminp1 <- min(result1[result1 > 0.1], na.rm=TRUE)
      cat("Lminp1 =", Lminp1)
      result2 <- OC.Isop(p2, nmin, kopt)
      plot(result2[result2>0.02])
      Lmaxp2 <- max(result2, na.rm=TRUE)
      cat(" Lmaxp2 =", Lmaxp2, "\n")
      if ( (Lminp1 >= 1-alpha) & (Lmaxp2 <= beta))
        {### Plan gefunden!
          PlanGefunden <- TRUE 
        }
      else if ( !(Lminp1 >= 1-alpha ) & !(Lmaxp2 < beta) )
        {### Es gibt keine Lösung
          nmin <- nmin+1
          direction<-"unknown"
          cat("nächstes n: ", nmin, "\n")
        }
      else
        { ### es muss am k gerüttelt werden.
          if (Lminp1 < 1-alpha) {
            if (direction == "unknown") direction <- "up"
            if (direction == "up") kopt <- kopt + kdiff else {
              kdiff <- kdiff/2 ; kopt <- kopt + kdiff}
            cat("neues kopt: " , kopt, "\n")
          } else
          {
            if (direction == "unknown") direction="down"
            if (direction == "down") kopt <- kopt - kdiff else {
              kdiff <- kdiff/2 ; kopt <- kopt - kdiff}
            cat("neues kopt: " , kopt, "\n")
          }
        }
    }
  }

fOptPlan.quick <- function( p1, p2, alpha, beta)
  {
    startplan <- LRP(p1, p2, alpha, beta)
    nmin <- startplan$ntilde
    kopt <- startplan$ktilde
    kdiff <- 2e-4
    PlanGefunden <- FALSE
    direction <- "unknown"
    
    while(! PlanGefunden) {
      cat("Teste Plan mit p1=" , p1, " p2= ", p2, " n= ", nmin, " und k=", kopt, "\n")
      Lminp1 <- OC.Isop.quick(p1, nmin, kopt)
      cat("Lminp1 =", Lminp1)
      Lmaxp2 <- OC.Isop.quick(p2, nmin, kopt)
      cat(" Lmaxp2 =", Lmaxp2, "\n")
      if ( (Lminp1 >= 1-alpha) & (Lmaxp2 <= beta))
        {### Plan gefunden!
          PlanGefunden <- TRUE 
        }
      else if ( !(Lminp1 >= 1-alpha ) & !(Lmaxp2 < beta) )
        {### Es gibt keine Lösung
          nmin <- nmin+1
          direction<-"unknown"
          cat("nächstes n: ", nmin, "\n")
        }
      else
        { ### es muss am k gerüttelt werden.
          if (Lminp1 < 1-alpha) {
            if (direction == "unknown") direction <- "up"
            if (direction == "up") kopt <- kopt + kdiff else {
              kdiff <- kdiff/2 ; kopt <- kopt + kdiff
            }
            cat("neues kopt: " , kopt, "\n")
          } else
          {
            if (direction == "unknown") direction <- "down"
            if (direction == "down") kopt <- kopt - kdiff else {
              kdiff <- kdiff/2 ; kopt <- kopt - kdiff
            }
            cat("neues kopt: " , kopt, "\n")
          }
        }
    }
  }


fOptPlan2 <- function( p1, p2, alpha, beta, kdiff = 2e-5, verbose = FALSE, show.plots=TRUE)
  {
    startplan <- LRP(p1, p2, alpha, beta)
    n <- startplan$ntilde
    kopt <- startplan$ktilde
    PlanGefunden <- FALSE
    direction <- "unknown"
    
    while(! PlanGefunden) {
      if (verbose) cat("Teste Plan mit p1=" , p1, " p2= ", p2, " n= ", nmin, " und k=", kopt, "\n")
      result1 <- OC.Isop(p1, nmin, kopt)
      if (show.plots) plot(result1[result1 > 0.1])
      Lminp1 <- min(result1[result1 > 0.1], na.rm=TRUE)
      if (verbose) cat("Lminp1 =", Lminp1)
      result2 <- OC.Isop(p2, nmin, kopt)
      if (show.plots) plot(result2[result2>0.02])
      Lmaxp2 <- max(result2, na.rm=TRUE)
      if (verbose) cat(" Lmaxp2 =", Lmaxp2, "\n")
      if ( (Lminp1 >= 1-alpha) & (Lmaxp2 <= beta))
        {### Plan gefunden!
          PlanGefunden <- TRUE
          LRPlan <- list(n=n, Lminp1=Lminp1, Lmaxp2=Lmaxp2)
        }
      else if ( !(Lminp1 >= 1-alpha ) & !(Lmaxp2 < beta) )
        {### Es gibt keine Lösung
          n <- n+1
          direction<-"unknown"
          if (verbose) cat("nächstes n: ", nmin, "\n")
        }
      else
        { ### es muss am k gerüttelt werden.
          if (Lminp1 < 1-alpha) {
            if (direction == "up") kopt <- kopt + kdiff else {kdiff <- kdiff/2 ; kopt <- kopt + kdiff; direction <- "up"}
            if (verbose) cat("neues kopt: " , kopt, "\n")
          } else
          {
            if (direction == "down") kopt <- kopt - kdiff else {kdiff <- kdiff/2 ; kopt <- kopt - kdiff; direction <- "down"}
            cat("neues kopt: " , kopt, "\n")
          }
        }
    }

    return(LRPlan)
    
  }



mup <- function(mu,plevel,sigma){pnorm((-1-mu)/sigma) + pnorm((mu-1)/sigma) - plevel}
### diese Funktion ist nur eine Hilfsfunktion für die numerische Nullstellensuche
### über den Parameter mu gegeben plevel und sigma. Funktioniert.

sigmap0 <- function(p){-1/qnorm(p/2) }
### Dies ist das maximal interessante sigma zu einem gegebenen p-level,
### korrespondiert zu einem mu=0

OCML <- function(mu, sigma, n, k) {
  B <- (n-1)/(sigma * qnorm(k/2))^2
  cat("Parameter für integrate: B ", B, " n ", n, " mu", mu, " sigma ", sigma, " k ", k,  "\n")
 result <- quadgk(i2OCML, 0 , B, tol=1e-5, n=n, mu=mu , sigma=sigma, k=k)
###  result <- integrate(integrandOCML, 0 , B, n=n, mu=mu , sigma=sigma, k=k, subdivisions=100)
###  cat("B ", B, " subdivisions: " , result$subdivisions," Wert : ", result$value, "\n")
  cat("B ", B, " Wert : ", result, "\n")
###  result$value
  ### für integrate
  result
}
### Das ist L_3 aus dem Paper. Die OC für den ML Ansatz.
### Die Probleme treten bei der Integration auf.

integrandOCML <- function(t, n, mu, sigma, k){
  t[1] <- min( 1e-6, (t[1]+t[2]) /2)
  cat("intergrand OCML aufgerufen mit t = ", t, " länge ", length(t), "\n")
  cat(" n = ", n)
  cat(" mu = ", mu)
  cat(" sigma = ", sigma)
  cat(" k = ", k, "\n")
  musigma <- sigma*sqrt(t/(n-1))
   
  lokalmu <- rep(NA, length(musigma))
  for ( i in 1:length(musigma)) {
    cat("Parameter für uniroot sigma", musigma[i], "plevel ",k, "\n")
    lokalmu[i] <- uniroot(mup, interval=c(1e-6, 1 ) , sigma=musigma[i], plevel=k )$root
 cat("lokalmu [ ",i," ]= ", lokalmu[i], "\r")
  }
    dchisq(t,df=n-1) * (pnorm( sqrt(n)*(lokalmu - mu )/sigma )-
       pnorm( sqrt(n)*(-lokalmu -mu )/sigma))        
  }

i2OCML <- function(t,n,mu,sigma,k){
  t[1] <- min( 1e-6, (t[1]+t[2]) /2)
  ## cat("intergrand OCML aufgerufen mit t = ", t, " länge ", length(t), "\n")
  ## cat(" n = ", n)
  ## cat(" mu = ", mu)
  ## cat(" sigma = ", sigma)
  ## cat(" k = ", k, "\n")
  musigma <- sigma*sqrt(t/(n-1))
   
  lokalmu <- rep(NA, length(musigma))
  for ( i in 1:length(musigma)) {
    ## cat("Parameter für uniroot sigma", musigma[i], "plevel ",k, "\n")
    lokalmu[i] <- uniroot(mup, interval=c(1e-6, 1 ) , sigma=musigma[i], plevel=k )$root
    ###cat("lokalmu [ ",i," ]= ", lokalmu[i], "\r")
  }
    dchisq(t,df=n-1) * (pnorm( sqrt(n)*(lokalmu - mu )/sigma )-
       pnorm( sqrt(n)*(-lokalmu -mu )/sigma))        
  }


OC.Isop.ML <- function( p, n, k, resolution=100)
  {
###    cat("OC.Isop( p=", p, " n=", n, " k=", k, ")\n")
    sigma.max <- sigmap0(p)
    resolution <- 100
    sigmas<- seq(1e-6, sigma.max*(1-1e-5), length.out=resolution)
    mus <- rep(NA, resolution)
    integralwert<- mus
    
    for (i in 1:resolution){
      mus[i] <- uniroot(mup, interval=c(0,1), tol=1e-15, sigma=sigmas[i], plevel=p)$root
    }
   for (i in 1:resolution){
      cat( "Aufruf in Oc.Isop.ML mit Parametern: mus[i] ", mus[i], "sigma[i] ", sigmas[i], " ",i," n: ", n,"k: ",k ,"\n")
       integralwert[i] <- OCML(mus [i], sigmas[i], n, k )
    }
    return(integralwert)
  }

