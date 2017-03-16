#' cFDR: a package for estimating conditional false discovery rates
#'
#' The main functions are cFDR and FDR
#' 
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
FDR <- function(Z,alt,alpha,p.max) {
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
##' @param alpha.target optional (default 0). If set, will return the absolute difference
##'     between the estimated FDR and alpha.  Used if you want to
##'     optimise alpha to a pre-specified overall FDR control
##' @param method method for estimating overall FDR in a set of
##'     overlapping rectangles given known FDR control in each
##' @export
##' @return absolute difference between estimated overall FDR control and alpha.target
##' @author Chris Wallace
cFDR.tune <- function(alpha.current=alpha, p,pc,pp,alpha.target=0,
                      method=c("liley.number","constant.density","liley.area")) {
    method <- match.arg(method)
    message(alpha.current)
    reject <- pp<=alpha.current #& p<=p.max
    nL <- sum(reject)
    if(nL==0)
        return(0)
    allM <- cbind(p,pc,A=p*pc)[reject,,drop=FALSE]
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
    nM <- sum(p<x & pc < y)
    alpha.obs <- switch(method,
                        liley.area=alpha.current * areaL/areaM,
                        liley.number=alpha.current * nL/nM,
                        constant.density = alpha.current * nM * areaL / ( nL * areaM),
                        conservative=1-(1-alpha.current)*nM/nL)
    abs(alpha.obs - alpha.target)
}

##' Calculate cFDR given a vector of test and conditional Z scores
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
    if(do.optimise) {
        if(length(alpha)>1)
            stop("can only optimise to a single alpha")
        o <- optimise(cFDR.tune, p=p, pc=pc, pp=pp, ..., alpha.target=alpha, interval=alpha*c(0.001,1),tol=eps)
        return(list(alpha=o$minimum, est.cFDR=alpha, cFDR=pp))
    }
    est.cFDR <- sapply(alpha, cFDR.tune, p=p, pc=pc, pp=pp,...)
    list(alpha=alpha, est.cFDR=est.cFDR, cFDR=pp)
}

##' Assess power and FDR control for the cFDR procedure for simulated
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
cFDR.assess <- function(Z,Zc,alt,alpha,...) {
    result <- cFDR(Z,Zc,alpha=alpha,do.optimise=TRUE)
    reject <- result$cFDR < result$alpha #& p < p.max
    power=if(any(alt)) { mean(reject[ alt ]) } else { 0 }
    fdr=mean(!alt[reject])
    c(alpha.input=alpha,alpha=result$alpha,power=power,fdr=fdr)
}
