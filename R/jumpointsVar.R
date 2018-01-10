#' @importFrom grDevices grey
#' @importFrom graphics abline par plot points
#' @importFrom stats Gamma coef dgamma glm.fit quantile resid sd smooth.spline
#' @import MASS lars

#' @export
jumpointsVar <- function(y, x=NULL, y.res=FALSE, k=min(30, round(length(y)/10)), print.level=0, plot.it=FALSE,
                         psi=NULL, round=TRUE, control=fit.control(), selection=sel.control()){
    # require(MASS)
    # require(lars)
    # require(cumSeg)

    call <- match.call()
    fam <- Gamma(link="log")
    input <- list(y=y, x=x, y.res=y.res, k=k, psi=psi, round=round,
                  control=control, selection=selection, fam=fam)
    n <- length(y)

    reg <- y.start <- NULL
    if(!y.res){
      y.start <- y
      x <- 1:n
      reg <- smooth.spline(x, y)
      mu <- reg$y
      h <- reg$lev
      if(control$scale.res){
        y <- ((y-mu)/sd(y-mu))^2/(1-h) + 1 #
      }else{
        y <- (y-mu)^2/(1-h) + 1 #(y-mu)^2+(var(y)*(h))+1
      }
    }
    if(is.null(x)){
      x <- 1:n
      Y <- cumsum(y)
    }else{
      if(length(x) != n) stop("Lengths of x and y differ")
      y <- y[order(x)]
      x <- sort(x)
      diffx <- c(x[1], diff(x))
      Y <- cumsum(y*diffx)
    }
    if(is.null(psi)) psi <- quantile(x, prob=seq(0, 1, l=k+2)[-c(1, k+2)], names=FALSE)
    k <- length(psi)
    Z <- matrix(rep(x, k), nrow=n)
    PSI <- matrix(rep(psi, rep(n, k)), ncol=k)

    it.max <- old.it.max <- control$it.max
    if(it.max == 0) U <- pmax((Z - PSI), 0)

    obj <- seg.lm.fit0.gamma(y=Y, Z=Z, PSI=PSI, control=control, round=round, print.level=print.level)
    if(is.null(obj$psi)) stop("No estimated breakpoint in obj")

    # obj$Y <- Y
    # obj$y <- y
    # obj$y.start <- y.start

    # if(!is.null(obj$psi)){
    #   obj$fitted.values <- drop(abs(obj$V) %*% obj$beta.c)
    #   if("firstSlope" %in% names(coef(obj))) obj$fitted.values <- obj$fitted.values + obj$coef["firstSlope"]
    #   obj$id.group <- -rowSums(obj$V)
    # }

    # if(!is.null(obj$psi)){
    #   fitted.values <- drop(abs(obj$V) %*% obj$beta.c)
    #   if("firstSlope" %in% names(coef(obj))) fitted.values <- fitted.values + obj$coef["firstSlope"]
    #   id.group <- -rowSums(obj$V)
    # }else{
    #   stop("No estimated breakpoint in obj")
    # }

    # obj$est.means <- cumsum(c(obj$coef["firstSlope"], obj$beta.c))
    est.means0 <- cumsum(c(obj$coef["firstSlope"], obj$beta.c))
    psi0 <- obj$psi
    # est.means0 <- obj$est.means

    display1 <- if(print.level==2) TRUE else FALSE
    if(print.level == 1 | print.level == 2) cat("\n")
    edf.psi <- selection$edf.psi
    type <- selection$type
    S <- selection$S
    tipoAlg <- selection$alg
    Cn <- eval(parse(text=selection$Cn))
    if(control$scale.res) Cn <- Cn/sd(resid(reg))

    olars <- lars(abs(obj$V), y=log(y), type=tipoAlg, normalize=FALSE, intercept=TRUE, trace=display1)
    id.var.entry <- (1:ncol(obj$V))[order(olars$entry)]
    edf <- if(edf.psi) (olars$df-1)*2+1 else olars$df
    RSS <- olars$RSS
    nRSS <- length(RSS)
    LL <- vector(length=nRSS)
    xv <- matrix(1, n, 1)
    mod1 <- glm.fit(xv, y, family=fam)
    eta1 <- (xv) %*% (mod1$coefficients)
    fit <- fam$linkinv(eta1)
    phi1 <- 2
    LL[1] <- -2*sum(dgamma(y, 1/phi1, scale=fit*phi1, log=TRUE))
    for(i in 1:(nRSS-1)){
      xv <- cbind(1, obj$V[, id.var.entry[1:i]])
      mod1 <- glm.fit(xv, y, family=fam)
      eta1 <- (xv) %*% (mod1$coefficients)
      fit <- fam$linkinv(eta1)
      # phi1 <- 2
      LL[i+1] <- -2* sum(dgamma(y, 1/phi1, scale=fit*phi1, log=TRUE))
    }

    crit <- switch(type,
                   bic = (LL + log(n)*edf*Cn),
                   # mdl = (n*log(RSS/(n - edf)) + Cn*pen.MDL(id.var.entry, as.numeric(table(-rowSums(obj$V))))),
                   rss = (RSS))

    min.r <- switch(type,
                    bic = which.min(crit),
                    # mdl = which.min(crit),
                    rss = max.rss(crit, S, n))

    if(plot.it){
      plot(1:nRSS, crit, type="b", pch=19, xlab="Number of breakpoints", ylab=input$selection$type, frame.plot=FALSE)
      points(min.r, crit[min.r], col=2, pch=19)
    }

    id <- sort(c(0, id.var.entry)[1:min.r])
    if(length(id) <= 1){
      out <- list(y=y, est.means=mean(y), id.var.entry=id.var.entry, n.psi=0, criterion=crit,
                input=input, call=call)
      class(out) <- "jumpointsVar"
      return(out)
    }

    id <- id[-1]
    psi1 <- obj$psi[id]
    k <- length(psi1)

    Z <- matrix(rep(x, k), nrow=n)
    PSI <- matrix(rep(psi1, rep(n, k)), ncol=k)
    obj <- seg.lm.fit0.gamma(y=Y, Z=Z, PSI=PSI, control=control,#fit.control(toll=1e-04, it.max=1, stop.if.error=FALSE),
                             round=round, print.level=print.level)
    if(!is.null(obj$est.means)){
      obj$n.psi <- 0
      obj$psi0 <- psi0
      obj$psi1 <- psi1
      obj$criterion <- crit
      return(obj)
    }

    fitted.v <- drop(abs(obj$V) %*% obj$beta.c)
    if("firstSlope" %in% names(coef(obj))) fitted.v <- fitted.v + obj$coef["firstSlope"]
    id.group <- -rowSums(obj$V)
    est.means <- cumsum(c(obj$coef["firstSlope"], obj$beta.c))
    out <- list(psi=psi1, est.means=est.means, n.psi=length(psi1), id.var.entry=id.var.entry,
                psi.order.entry=psi0[id.var.entry], psi0=psi0, est.means0=est.means0, criterion=crit,
                fitted.values=fitted.v, input=input, call=call)
    class(out) <- "jumpointsVar"
    return(out)
}

#' @export
print.jumpointsVar <- function(x, digits=max(3L, getOption("digits") - 3L), ...){
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nNumber of breakpoints:", x$n.psi, "/", x$input$k, "\n", sep="")
  cat("IDs of breakpoints:", as.vector(x$psi), "\n")
  # cat("Estimated means:\n")#, signif(round(x$est.means, digits), digits), "\n")
  # print.default(signif(round(x$est.means, digits), digits), print.gap=2L, quote=FALSE)

  cat("\n")
  invisible(x)
}

#' @export
plot.jumpointsVar <- function(x, ...){
  psi <- x$psi
  y <- x$input$y
  est.means <- x$est.means

  par(mfrow=c(1, 1))
  plot(y, type="l", col=grey(0.7), ylab="y", xlab="segmented variable", ...)
  abline(v=psi, col=2, lty=2)

  return(invisible(x))
}

seg.lm.fit0.gamma <- function(y, Z, PSI, control, round=FALSE, print.level=0){
  fam <- Gamma(link="log")
  it.max <- old.it.max <- control$it.max
  toll <- control$toll
  visual <- print.level
  last <- control$last
  stop.if.error <- control$stop.if.error
  h <- min(abs(control$h), 1)
  if(h < 1) it.max <- it.max + round(it.max/2)
  it <- 1
  epsilon <- 10
  k <- ncol(PSI)
  psi.values <- NULL
  H <- 1
  psi <- PSI[1, ]
  XREG <- cbind(Z[, 1])
  obj <- list(residuals=rep(10, 3))
  while(abs(epsilon) > toll){
    U <- pmax((Z - PSI), 0)
    V <- ifelse((Z > PSI), -1, 0)
    X <- cbind(XREG, U, V)
    dev.old <- sum(obj$residuals^2)
    rownames(X) <- NULL
    if(ncol(V) == 1){
      colnames(X)[ncol(XREG):ncol(X)] <- c("firstSlope","U", "V")
    }else{
      colnames(X) <- rev(c(paste("V", ncol(V):1, sep = ""),
                           paste("U",ncol(U):1, sep = ""),
                           rep("firstSlope", ncol(X)-ncol(U)-ncol(V))))
    }
    obj <- glm.fit(x = X, y = y, family=fam, control=list(maxit=10000))
    dev.new <- obj$deviance
    if(visual == 1 | visual == 2){
      if(it == 1) cat(0, "", formatC(dev.old, 3, format="f"),"", "(No breakpoint(s))", "\n")
      spp <- if(it < 10) "" else NULL
      cat(it, spp, "", formatC(dev.new, 3, format="f"), "---", ncol(V), "breakpoints", "\n")
    }
    epsilon <- (dev.new - dev.old)/dev.old
    obj$epsilon <- epsilon
    it <- it + 1
    obj$it <- it
    class(obj) <- c("segmented", class(obj))
    if(k == 1){
      beta.c <- (coef(obj))["U"]
      gamma.c <- (coef(obj))["V"]
    }else{
      beta.c <- (coef(obj))[paste("U", 1:ncol(U), sep = "")]
      gamma.c <- (coef(obj))[paste("V", 1:ncol(V), sep = "")]
    }
    if(it > it.max) break
    psi.values[[length(psi.values) + 1]] <- psi.old <- psi
    if(it >= old.it.max && h < 1) H <- h
    psi <- round(psi.old + H * gamma.c/beta.c, 0)
    PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
    a <- apply((Z <= PSI), 2, all)
    b <- apply((Z >= PSI), 2, all)
    if(stop.if.error){
      if(sum(a + b) != 0 || is.na(sum(a + b))) stop("(Some) estimated psi out of its range")
    }else{
      id.psi.ok <- !is.na((a+b) <= 0) & (a+b) <= 0
      Z <- Z[, id.psi.ok, drop=FALSE]
      psi <- psi[id.psi.ok]
      PSI <- PSI[, id.psi.ok, drop=FALSE]
    }
    if(ncol(PSI) <= 0){
      warning("No breakpoint estimated", call. = FALSE)
      obj <- glm.fit(x=XREG, y=y, family=fam, control=list(maxit=10000))
      obj$fitted.values <- rep((obj$coef), length(y))
      obj$est.means <- (obj$coef)
      return(obj)
    }
  }
  if(round){
    psi <- round(psi, 0)
    PSI <- matrix(rep(psi, rep(nrow(Z), ncol(Z))), ncol = ncol(Z))
    V <- ifelse((Z > PSI), -1, 0)
  }

  obj$psi <- sort(psi)
  obj$beta.c <- beta.c[order(psi)]
  obj$gamma.c <- gamma.c[order(psi)]
  obj$epsilon <- epsilon
  obj$V<- V[, order(psi)]
  obj$psi <- obj$psi[!is.na(obj$beta.c)]
  obj$V <- as.matrix(as.matrix(obj$V)[, !is.na(obj$beta.c)])
  obj$beta.c <- obj$beta.c[!is.na(obj$beta.c)]
  return(obj)
}

max.rss <- function(RSS, S, n){
  var.mod <- RSS/n
  ll <- -(log(2*pi*var.mod)+1)*n/2
  new.r <- ((ll[length(ll)]-ll[-1])/(ll[length(ll)]-ll[2]))*(length(ll)-1) + 1
  diff2 <- diff(new.r, diff=2) > S
  if(!any(diff2)) return(0)
  maxll <- max(which(diff2)) + 1
  return(maxll)
}

#' @export
sel.control <- function(type=c("bic", "rss"), S=1,
                        Cn="2*log(log(n))", alg=c("lasso", "stepwise"), edf.psi=TRUE){
  # type=c("bic", "mdl"
  type <- match.arg(type)
  alg <- match.arg(alg)
  return(list(type=type, S=S, Cn=Cn, edf.psi=edf.psi, alg=alg))
}

#' @export
fit.control <- function(toll=0.001, it.max=10, last=TRUE, scale.res=FALSE,
                        maxit.glm=30, h=1, stop.if.error=FALSE){
  return(list(toll=toll, it.max=it.max, last=last, scale.res=FALSE, maxit.glm=maxit.glm, h=h,
              stop.if.error=stop.if.error))
}

# pen.MDL <- function(id, n){
#   # restituisce un vettore (di dim=length(id)) che rappresenta la penalità 2*sum\log n_j per ogni partizione (active set)
#   # length(id) è il num (max) di breakpoints ed n è il vettore delle numerosità della partizione.
#
#   do.m <- function(id, n.col){
#     blockdiag <- function(...){
#       args <- list(...)
#       nc <- sapply(args,ncol)
#       cumnc <- cumsum(nc)
#       ##  nr <- sapply(args,nrow)
#       ## NR <- sum(nr)
#       NC <- sum(nc)
#       rowfun <- function(m,zbefore,zafter) {
#         cbind(matrix(0, ncol=zbefore, nrow=nrow(m)), m, matrix(0, ncol=zafter, nrow=nrow(m)))
#       }
#       ret <- rowfun(args[[1]], 0, NC-ncol(args[[1]]))
#       for(i in 2:length(args)){
#         ret <- rbind(ret, rowfun(args[[i]], cumnc[i-1], NC-cumnc[i]))
#       }
#       ret
#     } #end blockgiag
#     id <- sort(id) #sort(unlist(id))
#     if(length(id) == 1){
#       m <- t(rep(1, id))
#     }else{
#       m <- do.call(blockdiag, lapply(c(id[1], diff(id)), function(xx) t(rep(1, xx))))
#     }
#     m <- blockdiag(m, t(rep(1, n.col-ncol(m))))
#     m
#   } #end do.m
#
#   #inizio codici veri
#   if(length(n) != (length(id)+1)) stop("Error in 'id' or 'n'")
#   A <- matrix(rev(id), length(id), length(id), byrow=FALSE)
#   A[col(A) > row(A)] <- NA
#   r <- rev(apply(A, 2, function(x) x[!is.na(x)]))
#   lista.m <- lapply(r, do.m, n.col=length(n))
#   #sapply(lista.m,function(xx)drop(xx%*%n))
#
#   ris <- sapply(lista.m, function(xx) 2*sum(log(drop(xx %*% n))))
#   ris <- c(2*sum(log(n)), ris)
#   ris
# }
