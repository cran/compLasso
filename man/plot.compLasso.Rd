\name{plot.compLasso}
\alias{plot.compLasso}
\title{
Plot the results of the component lasso method.
}
\description{
A function to plot the coefficient paths of the fitted component lasso object.
}
\usage{
\method{plot}{compLasso}(x, xvar = c("norm", "lambda"), label = FALSE, ...)
}
\arguments{
  \item{x}{
Fitted "compLasso" object.
}
  \item{xvar}{
Represents what is on the x-axis of the plot. "norm" plots the coefficient paths against the L1 norm of the coefficients and "lambda" against the log-lambda sequence of values of the fitted object.
}
  \item{label}{
If TRUE, this labels the curves with variable sequence numbers.
}
  \item{\dots}{
Other graphical parameters to plot.
}
}
\details{
This function plots the coefficient paths of the component lasso versus the L1 norm of the coefficient vectors or the log(lambda) values. 

}
\references{Hussami, N and Tibshirani, R (2013) A Component Lasso.}
\author{FORTRAN code by Jerry Friedman. R interface by Nadine Hussami and Robert Tibshirani}
\note{The FORTRAN code that this function links to was kindly written
  and provided by Jerry Friedman. This function is taken from the glmnet package.
}
\seealso{ compLasso, cv.compLasso}
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

a=compLasso(x,y,clus,nnls=TRUE)
plot(a,xvar="norm",label=TRUE)

}
