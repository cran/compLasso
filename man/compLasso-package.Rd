\name{compLasso-package}
\alias{compLasso-package}
\docType{package}
\title{
Implements the component lasso regression method.
}
\description{
Functions for implementing the component lasso method for regression.
}
\details{
\tabular{ll}{
Package: \tab compLasso\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-10-04\cr
License: \tab GPL (>=2)\cr
}

The main function is "compLasso", which takes a data matrix x as input and a response vector y. It estimates the regression coefficients for a range of tuning parameter values which can also be supplied to the function. The package also contains functions to perform prediction and cross-validation. 

}
\author{
Nadine Hussami, Robert Tibshirani, and Jerome Friedman.

Maintainer: Nadine Hussami <nadinehu@stanford.edu>
}
\references{
Hussami and Tibshirani (2013) A Component Lasso.
}
\keyword{ package }
\examples{
set.seed(1000)

n=200
p=50
k=5
clus=rep(1:k,p/k)

x=matrix(rnorm(n*p),ncol=p)
b0=rep(0,p)
b0[clus<5]=3
y=x\%*\%b0+.1*rnorm(n)

a=compLasso(x,y,clus,nnls=TRUE,nlambda=10)
predict(a,newx=x[1:10,],s=c(0.1,0.02))

foldid=c(rep(2,40),rep(1,40),rep(4,40),rep(5,40),rep(3,40))
a=cv.compLasso(x, y, clus, parm=1, foldid=foldid, type.measure = "mse", nfolds = 5,nnls=TRUE) 
}
