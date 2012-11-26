OC <- function ( mu, sigma, n, k, resolution=100, type="LR")
  {
    if (!( (type=="LR") | (type == "BSK") ) ) stop("Type must be either \"LR\" or \"BSK\"!")
    if (type == "LR") return (OC.1( mu, sigma, n, k ) + 2 * OC.2( mu, sigma, n, k ) )
    if (type == "BSK") {
       B <- (n-1)/(sigma * qnorm(k/2))^2
       return (integrate(iOCML, 0 , B, n=n, mu=mu , sigma=sigma, k=k, subdivisions=resolution)$value)
    }
  }

OCIsop <- function( p, n, k, resolution=100, type="LR"){
  if (!( (type=="LR") | (type == "BSK") ) ) stop("Type must be either \"LR\" or \"BSK\"!")
  sigma.max <- sigmap0(p)
  sigmas<- seq(0,sigma.max,length.out=resolution+2)[2:(resolution+1)]
  mus <- rep(NA, resolution)
  integralwert<- mus
  
  for (i in 1:resolution){
    mus[i] <- uniroot(mup, interval=c(0,1), tol=1e-15, sigma=sigmas[i], plevel=p)$root
  }
  for (i in 1:resolution){
    integralwert[i] <- OC(mus [i], sigmas[i], n, k, type=type )
  }
  return(integralwert)
}

singleSidedApproxPlan <- function(p1, p2, alpha, beta, verbosity.level=0)
  {
### Einseitiger Naeherungsplan; ist eine gute Approximation für
### einen Startplan, um nach BSK oder LR zu suchen
### Es wird das minimale n zurückgegeben, für das eine Loesung exisitert.

    n1 <- ceiling(((qnorm(1-alpha) - qnorm(beta)) / (qnorm(p2)-qnorm(p1)) )^2 )
    l1n <- qnorm(1-alpha) + sqrt(n1)*qnorm(p1)
    l2n <- qnorm(beta) + sqrt(n1)*qnorm(p2)
    l1 <- (l1n+l2n)/2
    
    n0 <- floor( n1*(1 + l1^2/(2*n1)) )
    l <- l1 * sqrt(1 + l1^2/(2*n1))

   if (verbosity.level >= 1) cat("Normal approx. of single sided LR plan:  n_0 =",n0, " l_0 =", l, "\n" )
   return(list(ntilde=n0, ltilde=l,ktilde=pbeta(0.5 + l/(2*(n0-1)), (n0-2)/2, (n0-2)/2) ))
  }

singleSidedLRPlan <- function(p1, p2, alpha, beta, verbosity.level=0)
  {
    startplan <- singleSidedApproxPlan( p1, p2, alpha, beta, verbosity.level=verbosity.level)
    nmin <- startplan$ntilde - 1
    lopt <- startplan$ltilde

    ldiff <- 0.02 * abs(lopt)
        PlanGefunden <- FALSE
    direction <- "unknown"
    
    while(! PlanGefunden) {
      if (verbosity.level >= 2) cat("Testing plan  p1 = " , p1, " p2 = ", p2, " n = ", nmin, " l =", lopt, "\n")
      result1 <- pt(lopt, nmin-1, sqrt(nmin)*qnorm(p1) )
      if (verbosity.level >= 2) cat("condition (i) ", result1)

      result2 <- pt(lopt, nmin-1, sqrt(nmin)*qnorm(p2) )
      if (verbosity.level >= 2) cat(" condition (ii)", result2, "\n")

      if ( (result1 >= 1-alpha) & (result2 <= beta))
        {### Plan gefunden!
          PlanGefunden <- TRUE
         
          plan <- list(ntilde=nmin, ltilde=lopt, ktilde=pbeta(0.5 + lopt/(2*nmin -2), (nmin-2)/2, (nmin-2)/2 ))
        }
      else if ( !(result1 >= 1-alpha ) & !(result2 < beta) )
        {### Es gibt keine Loesung
          nmin <- nmin+1
          direction<-"unknown"
          if (verbosity.level >= 2) cat("increase n now: ", nmin, "\n")
        }
      else
        { ### es muss am k geruettelt werden.
          if (result1 < 1-alpha) {
            if (direction == "unknown") direction <- "up"
            if (direction == "up") lopt <- lopt + ldiff else {
              direction <- "up" ; ldiff <- ldiff/2 ; lopt <- lopt + ldiff}
            if (verbosity.level >= 2) cat("new l_opt: " , lopt, "\n")
          } else
          {
            if (direction == "unknown") direction="down"
            if (direction == "down") lopt <- lopt - ldiff else {
              direction <- "down" ; ldiff <- ldiff/2 ; lopt <- lopt - ldiff}
            if (verbosity.level >= 2) cat("new l_opt: " , lopt, "\n")
          }
        }
    }
    if (verbosity.level >= 1) cat("Single sided LR plan:  n_tilde =",plan$ntilde, " l_tilde =", plan$ltilde," k_tilde = ", plan$ktilde, "\n" )
    return(plan)
  }

findOptPlan <- function( p1, p2, alpha, beta, type="LR", use.quickApproxPlan=TRUE, given.plan=NULL, verbosity.level=0, show.plots=FALSE) {
  if (!( (type=="LR") | (type == "BSK") ) ) stop("Type must be either \"LR\" or \"BSK\"!")

  if (use.quickApproxPlan) {
    given.plan <- quickApproxPlan( p1, p2, alpha, beta, type=type, verbosity.level=verbosity.level)
  }

if (is.null(given.plan)) {
  startplan <- singleSidedLRPlan(p1, p2, alpha, beta, verbosity.level=verbosity.level)
  nmin <- startplan$ntilde
  if (type == "LR") kopt <- startplan$ktilde
  if (type == "BSK") kopt <- pnorm(startplan$ltilde/sqrt(startplan$ntilde))
}
else {
  nmin <- given.plan$ntilde
  kopt <- given.plan$ktilde
}
kdiff <- 0.01*kopt
PlanGefunden <- FALSE
direction <- "unknown"
    
while(! PlanGefunden){
  if (verbosity.level >= 2) cat("Testing plan  p1 = " , p1, " p2 = ", p2, " n = ", nmin, " k =", kopt, "\n")
  result1 <- OCIsop(p1, nmin, kopt, type=type)
  if (show.plots) plot(result1[result1 > 0.1])
  Lminp1 <- min(result1[result1 > 0.1], na.rm=TRUE)
  if (verbosity.level >= 2) cat("L_min(p1) =", Lminp1)
  result2 <- OCIsop(p2, nmin, kopt, type=type)
  if (show.plots) plot(result2[result2>0.02])
  Lmaxp2 <- max(result2, na.rm=TRUE)
  if (verbosity.level >= 2) cat(" L_max(p2) =", Lmaxp2, "\n")
  
  if ( (Lminp1 >= 1-alpha) & (Lmaxp2 <= beta)){ ### Plan gefunden!
    PlanGefunden <- TRUE
    plan <- list(ntilde=nmin, ltilde=NA, ktilde=kopt)
  }
  else if ( !(Lminp1 >= 1-alpha ) & !(Lmaxp2 < beta) ) { ### Es gibt keine Loesung
    nmin <- nmin+1
    direction<-"unknown"
    if (verbosity.level >= 2) cat("naechstes n: ", nmin, "\n")
  }
  else { ### es muss am k geruettelt werden.
    if (Lminp1 < 1-alpha) {
      if (direction == "unknown") direction <- "up"
      if (direction == "up") kopt <- kopt + kdiff else {
        direction <- "up" ; kdiff <- kdiff/2 ; kopt <- kopt + kdiff
      }
      if (verbosity.level >= 2) cat("new k_opt: " , kopt, "\n")
    } else {
      if (direction == "unknown") direction="down"
      if (direction == "down") kopt <- kopt - kdiff else {
        direction <- "down" ; kdiff <- kdiff/2 ; kopt <- kopt - kdiff
      }
      if (verbosity.level >= 2) cat("new k_opt: " , kopt, "\n")
    }
  }
}
if (verbosity.level >= 1) cat("Optimal plan of type ",type , " p_1 = ", p1, " p_2 = ", p2, " alpha = ", alpha, " beta = ", beta, ": n = ", plan$ntilde, " k = ", plan$ktilde, "\n")
invisible(plan)
}

quickApproxPlan <- function( p1, p2, alpha, beta, type="LR", verbosity.level=0){

  if (!( (type=="LR") | (type == "BSK") ) ) stop("Type must be either \"LR\" or \"BSK\"!")
  
  startplan <- singleSidedLRPlan(p1, p2, alpha, beta, verbosity.level=verbosity.level)
  nmin <- startplan$ntilde
  if (type == "LR") { kopt <- startplan$ktilde }
  if (type == "BSK") { kopt <- pnorm(startplan$ltilde/sqrt(startplan$ntilde)) }
  
  kdiff <- 0.02 * kopt
  
  PlanGefunden <- FALSE
  direction <- "unknown"
    
  while(! PlanGefunden) {
    if (verbosity.level >= 2) cat("Testing plan p1=" , p1, " p2 = ", p2, " n = ", nmin, " k =", kopt, "\n")
    Lminp1 <- quickOCIsop(p1, nmin, kopt, type=type)
    if (verbosity.level >= 2) cat("L_min(p1) =", Lminp1)
    Lmaxp2 <- quickOCIsop(p2, nmin, kopt, type=type)
    if (verbosity.level >= 2) cat(" L_max(p2) =", Lmaxp2, "\n")
    if ( (Lminp1 >= 1-alpha) & (Lmaxp2 <= beta))
      {### Plan gefunden!
        PlanGefunden <- TRUE
        plan <- list(ntilde=nmin, ltilde=NA, ktilde=kopt)
      }
    else if ( !(Lminp1 >= 1-alpha ) & !(Lmaxp2 < beta) )
      {### Es gibt keine Loesung
        nmin <- nmin+1
        direction<-"unknown"
        if (verbosity.level >= 2) cat("increase n now: ", nmin, "\n")
      }
    else
      { ### es muss am k geruettelt werden.
        if (Lminp1 < 1-alpha) {
          if (direction == "unknown") direction <- "up"
          if (direction == "up") kopt <- kopt + kdiff else {
            direction <- "up" ; kdiff <- kdiff/2 ; kopt <- kopt + kdiff
          }
          if (verbosity.level >= 2) cat("new k_opt: " , kopt, "\n")
        } else
          {
            if (direction == "unknown") direction <- "down"
            if (direction == "down") kopt <- kopt - kdiff else {
              direction <- "down" ; kdiff <- kdiff/2 ; kopt <- kopt - kdiff
            }
            if (verbosity.level >= 2) cat("new k_opt: " , kopt, "\n")
          }
      }
  }
  if (verbosity.level >= 1) cat("Fast approximative plan for p_1 = ", p1, " p_2 = ", p2, " alpha = ", alpha, " beta = ", beta, ": n = ", plan$ntilde, " k = ", plan$ktilde, "\n")
  return(plan)
}

quickOCIsop <- function (p , n, k, type="LR")
  {
    sigma.max <- sigmap0(p)*(1-1e-5)
    mus <- uniroot(mup, interval=c(0,10), tol=1e-15, sigma=sigma.max, plevel=p)$root
    return(OC(mus, sigma.max,n,k, type=type))
  
  }
