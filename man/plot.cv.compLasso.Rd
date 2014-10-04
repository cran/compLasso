\name{plot.cv.compLasso}
\alias{plot.cv.compLasso}
\title{
Plot the component lasso cross-validation curve returned by cv.compLasso.
}
\description{
This function plots the cross-validation curve versus the lambda values, along with the upper/lower corresponding standard deviation curves.
}
\usage{
\method{plot}{cv.compLasso}(x, sign.lambda = 1, ...)
}
\arguments{
  \item{x}{
Fitted "compLasso" object.
}
  \item{sign.lambda}{
Plots the curve against the log(lambda) values or its negative (if sign.lambda=-1).}
  \item{\dots}{
Other graphical parameters to plot.
}
}
\references{Hussami, N and Tibshirani, R (2013) A Component Lasso.}
\author{FORTRAN code by Jerry Friedman. R interface by Nadine Hussami and Robert Tibshirani}
\note{The FORTRAN code that this function links to was kindly written
  and provided by Jerry Friedman. This function is taken from the glmnet package.
}
\seealso{ cv.compLasso}
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

foldid=c(rep(2,40),rep(1,40),rep(4,40),rep(5,40),rep(3,40))
a=cv.compLasso(x, y, clus, parm=1, foldid=foldid, type.measure = "mse", nfolds = 5) 
print(a$cvm)

plot(a)
}
