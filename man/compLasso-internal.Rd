\name{compLasso-internal}
\title{Internal compLasso functions}
\alias{nnls}
\alias{cv.elnet}
\alias{cvcompute}
\alias{error.bars}
\alias{plotCoef}
\alias{nonzeroCoef}
\alias{getmin}
\description{Internal functions used by compLasso.}
\usage{
nnls(x,y,SMALL=1e-8)
cv.elnet(outlist, lambda, x, y, weights, foldid, type.measure, keep = FALSE)
cvcompute(mat, weights, foldid, nlams) 
error.bars(x, upper, lower, width = 0.02, ...) 
plotCoef(beta, norm, lambda, df, label=FALSE, xvar=c("norm","lambda"),
xlab=iname, ylab="Coefficients", ...) 
nonzeroCoef(beta, bystep = FALSE) 
getmin(lambda, cvm, cvsd) 
}
\author{Nadine Hussami and Robert Tibshirani}
\details{These functions are not to be called by the user and all but nnls are taken from the glmnet package.}
\keyword{internal}




