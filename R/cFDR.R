#' cFDR: a package for estimating conditional false discovery rates
#'
#' The main functions are cFDR and FDR
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom spatstat owin union.owin area
#' @importFrom Rcpp evalCpp
#' @useDynLib cFDR
#' @docType package
#' @name cFDR
NULL


##' Calculate FDR given a vector of Z scores
##' 
##' @title FDR
##' @param Z vector of Z scores
##' @export
##' @return vector of local FDR
##' @author Chris Wallace
FDR <- function(Z) {
    p=2*pnorm(abs(Z),lower.tail=FALSE)
    f <- ecdf(p)
    p/f(p)
}
##' Assess power and FDR control for simulated data where it is known whether each Z score corresponds to the null or alternative hypothesis
##'
##' .. content for \details{} ..
##' @title assess power and FDR control of FDR procedure
##' @param Z vector of Z scores
##' @param alt logical vector of equal length with Z, with TRUE indicating alternative hypotheses
##' @param alpha FDR threshold for assessment
##' @export
##' @return named vector of alpha, power, FDR
##' @author Chris Wallace
FDR.assess <- function(Z,alt,alpha) {
    pp <- FDR(Z)
    reject <- pp<=alpha
    power=if(any(alt)) { mean(reject[ alt ]) } else { 0 }
    fdr=mean(!alt[reject])
    ret <- c(alpha=alpha,power=power,fdr=fdr)
    names(ret) <- c("alpha","power","fdr")
    return(ret)
 
}

##' Estimate overal FDR for the cFDR procedure at a given alpha
##'
##' @title cFDR.tune
##' @param alpha.current alpha to assess
##' @param p p values corresponding to test statistics
##' @param pc p values corresponding to conditional test statistics
##' @param pp local FDR calculated by FDR function for p
##' @param alpha.target optional (default 0). If set, will return the
##'     absolute difference between the estimated FDR and alpha.  Used
##'     if you want to optimise alpha to a pre-specified overall FDR
##'     control
##' @param method method for estimating overall FDR in a set of
##'     overlapping rectangles given known FDR control in each.
##'     default is liley.area, which is the published method - the
##'     others are purely experimental and shouldn't be used.
##' @export
##' @return absolute difference between estimated overall FDR control
##'     and alpha.target
##' @author Chris Wallace
cFDR.tune <- function(alpha.current=alpha, p,pc,pp,alpha.target=0,
                      method=c("liley.area","constant.density","liley.number","unadjusted")) {
    method <- match.arg(method)
#    message(method,"\t",alpha.current)
    reject <- pp<=alpha.current #& p<=p.max
    nL <- sum(reject)
    if(nL==0)
        return(0)
    allM <- cbind(p,pc,A=p*pc)[reject,,drop=FALSE]
    biggestM <- which.max(allM[,3])
    x <- allM[biggestM,1]
    y <- allM[biggestM,2]
    nM <- sum(p<x & pc < y)
    alpha.obs <- switch(method,
                        unadjusted=alpha.current,
                        liley.area=alpha.current * areaLoverM(allM),
                        liley.number=alpha.current * nL/nM,
                        constant.density = alpha.current * nM * areaLoverM(allM) / nL,
                        conservative=1-(1-alpha.current)*nM/nL)
    abs(alpha.obs - alpha.target)
}

areaLoverM <- function(allM) {
    biggestM <- which.max(allM[,3])
    areaM <- allM[biggestM,3]
    x <- allM[biggestM,1]
    y <- allM[biggestM,2]
    allM <- allM[allM[,1]>=x | allM[,2] >=y,,drop=FALSE]
    polys <- lapply(1:nrow(allM), function(i) { owin(c(0, allM[i,1]), c(0, allM[i,2])) })
    merged.poly <- do.call("union.owin",polys)
    areaL <- area(merged.poly)
    if(areaL < areaM) # numerical inaccuracies
        areaL <- areaM
    areaL/areaM
}

##' Calculate cFDR given vectors of *independent* test and conditional Z scores
##'
##' @title cFDR
##' @inheritParams FDR.assess
##' @inheritParams cFDR.assess
##' @param eps tolerance for FDR control to alpha
##' @param alpha vector of thresholds. overall cFDR will be assessed
##'     if all tests with cFDR<alpha were rejected
##' @param do.optimise logical, default FALSE. if true, and a single
##'     value of alpha is supplied, the rejection threshold will be
##'     optimised to control the overall FDR at alpha.
##' @param ... other arguments passed to cFDR.tune
##' @return a named list with three components.  alpha is the
##'     threshold for rejection, est.cFDR is the estimated FDR at this
##'     threshold, cFDR is a vector of conditional FDRs in rectangles
##'     defined the corresponding absolute Z and Zc values.
##' @author Chris Wallace
cFDR <- function(Z,Zc,eps=0.01,alpha=seq(0.01,0.1,by=0.005),do.optimise=FALSE,...) {
    p=2*pnorm(abs(Z),lower.tail=FALSE)
    pc <- 2*pnorm(abs(Zc),lower.tail=FALSE)
    pp <- p / cecdf(p,pc)
    summary(pp)
    if(do.optimise) {
        if(length(alpha)>1)
            stop("can only optimise to a single alpha")
        o <- optimise(cFDR.tune, p=p, pc=pc, pp=pp, ..., alpha.target=alpha, interval=alpha*c(0.00001,1),tol=eps)
        return(list(alpha=o$minimum, est.cFDR=alpha, cFDR=pp))
    }
    est.cFDR <- sapply(alpha, cFDR.tune, p=p, pc=pc, pp=pp,...)
    list(alpha=alpha, est.cFDR=est.cFDR, cFDR=pp)
}

##' Calculate cFDR given vectors of *independent* test and conditional Z scores
##'
##' @title cFDR
##' @inheritParams FDR.assess
##' @inheritParams cFDR.assess
##' @param eps tolerance for FDR control to alpha
##' @param alpha vector of thresholds. overall cFDR will be assessed
##'     if all tests with cFDR<alpha were rejected
##' @param do.optimise logical, default FALSE. if true, and a single
##'     value of alpha is supplied, the rejection threshold will be
##'     optimised to control the overall FDR at alpha.
##' @param ... other arguments passed to cFDR.tune
##' @return a named list with three components.  alpha is the
##'     threshold for rejection, est.cFDR is the estimated FDR at this
##'     threshold, cFDR is a vector of conditional FDRs in rectangles
##'     defined the corresponding absolute Z and Zc values.
##' @author Chris Wallace
cFDR2 <- function(Z1,Z2,Zc1,Zc2,eps=0.01,alpha=seq(0.01,0.1,by=0.005),do.optimise=FALSE,...) {
    p1=2*pnorm(abs(Z1),lower.tail=FALSE)
    pc1 <- 2*pnorm(abs(Zc1),lower.tail=FALSE)
    p2=2*pnorm(abs(Z2),lower.tail=FALSE)
    pc2 <- 2*pnorm(abs(Zc2),lower.tail=FALSE)
    ## pp2 <- p2 / cecdf(p2,pc2)
    return(p1 * p2 / cecdf2(p1,p2,pc1,pc2))
    ## if(do.optimise) { ## NB this is a shortcut - should optimise product of cFDRs to alpha, not each to sqrt(alpha)
    ##     if(length(alpha)>1)
    ##         stop("can only optimise to a single alpha")
    ##     o1 <- optimise(cFDR.tune, p=p1, pc=pc1, pp=pp1, ..., alpha.target=sqrt(alpha), interval=sqrt(alpha)*c(0.001,1),tol=eps)
    ##     o2 <- optimise(cFDR.tune, p=p2, pc=pc2, pp=pp2, ..., alpha.target=sqrt(alpha), interval=sqrt(alpha)*c(0.001,1),tol=eps)
    ##     return(list(alpha=o1$minimum*o2$minimum, est.cFDR=alpha, cFDR=pp1*pp2))
    ## }
    ## est.cFDR1 <- sapply(sqrt(alpha), cFDR.tune, p=p1, pc=pc1, pp=pp1,...)
    ## est.cFDR2 <- sapply(sqrt(alpha), cFDR.tune, p=p2, pc=pc2, pp=pp2,...)
    ## list(alpha=alpha, est.cFDR=est.cFDR1*est.cFDR2, cFDR=pp1*pp2)
}

##' Assess power and FDR control for the *independent* cFDR procedure for simulated
##' data where it is known whether each Z score corresponds to the
##' null or alternative hypothesis
##'
##' @title cFDR
##' @inheritParams FDR.assess
##' @param Zc vector of Z scores for the conditional tests
##' @param ... other arguments passed to cFDR
##' @export
##' @return named vector of requested alpha, estimated alpha, power, fdr
##' @author Chris Wallace
cFDR.assess <- function(Z,Zc,alt,alpha,Z2=NULL,Zc2=NULL,...) {
    result <- if(is.null(Z2)) {
                  cFDR(Z,Zc,alpha=alpha,do.optimise=TRUE,...)
              } else {
                  cFDR2(Z1=Z,Zc1=Zc,Z2=Z2,Zc2=Zc2,alpha=alpha,do.optimise=FALSE,...)
              }
    reject <- result$cFDR < result$alpha 
    power=if(any(alt)) { mean(reject[ alt ]) } else { 0 }
    fdr=mean(!alt[reject])
    c(alpha.input=alpha,alpha=result$alpha,power=power,fdr=fdr)
}

##' Calculate cFDR given vectors of *dependent* test and conditional Z scores
##'
##' @title depcFDR
##' @inheritParams FDR.assess
##' @param rho correlation between Z and Zc under the null hypothesis
##' @param ... other arguments passed to cFDR.tune
##' @return a named list with three components.  alpha is the
##'     threshold for rejection, est.cFDR is the estimated FDR at this
##'     threshold, cFDR is a vector of conditional FDRs in rectangles
##'     defined the corresponding absolute Z and Zc values.
##' @author Chris Wallace
depcFDR <- function(Z,Zc,rho,mu=0,muc,eps=0.01,p.max=0.1,alpha=seq(0.01,0.1,by=0.005),do.optimise=FALSE,...) {
    p=2*pnorm(abs(Z),lower.tail=FALSE)
    pc <- 2*pnorm(abs(Zc),lower.tail=FALSE)
    upper <- rbind(-abs(Z),-abs(Zc))
    lower <- -upper
    sigma <- matrix(c(1,rho,rho,1),2)
    skip <- p>p.max
    ## microbenchmark(numer <- 2 * sapply(1:ncol(upper), function(i) { pmvnorm(upper=upper[,i]) + pmvnorm(lower=lower[,i]) }),
    ##                {
    denom <- numer <- rep(1,length(p))
    if(length(mu)==1 & length(muc)==1) {
        MU <- matrix(c(mu,muc),ncol=length(p),nrow=2)
    } else {
        MU <- rbind(mu,muc)
    }
    numer[!skip] <- sapply(which(!skip), function(i) {
        pmvnorm(upper=upper[,i],sigma=sigma,mean=MU[,i]) +
        pmvnorm(lower=lower[,i],sigma=sigma,mean=MU[,i]) +
        pmvnorm(lower=c(-Inf,lower[2,i]), upper=c(upper[1,i],Inf),sigma=sigma,mean=MU[,i]) +
        pmvnorm(lower=c(lower[1,i],-Inf), upper=c(Inf,upper[2,i]),sigma=sigma,mean=MU[,i]) })## },
    denom[!skip] <- sapply(which(!skip), function(i) {
        pnorm(upper[2,i],mean=MU[2,i]) +
        pnorm(lower[2,i],mean=MU[2,i],lower.tail=FALSE) })
    ## times=10)
    pp <- numer / (denom * cecdf(p,pc))
    ## if(do.optimise) {
    ##     if(length(alpha)>1)
    ##         stop("can only optimise to a single alpha")
    ##     o <- optimise(cFDR.tune, p=p, pc=pc, pp=pp, ..., alpha.target=alpha, interval=alpha*c(0.001,1),tol=eps)
    ##     return(list(alpha=o$minimum, est.cFDR=alpha, cFDR=pp))
    ## }
    ## est.cFDR <- sapply(alpha, cFDR.tune, p=p, pc=pc, pp=pp,...)
    list(alpha=alpha, ## est.cFDR=est.cFDR, 
         cFDR=pp)
}
