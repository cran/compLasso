\name{cv.compLasso}
\alias{cv.compLasso}
\title{
Evaluate the cross validation error of the component lasso.
}
\description{
Performs k-fold cross-validation for the component lasso and returns a vector of cross-validated errors and an optimal value of lambda.
}
\usage{
cv.compLasso(x, y, clus, lambda = NULL, type.measure = c("mse", "mae"),
weights, nfolds = 10,
foldid, keep = FALSE, ...)
}
\arguments{
  \item{x}{
Input matrix containing an observation per row.
}
  \item{y}{
Response variable.
}
  \item{clus}{
A vector of cluster labels.
}
  \item{lambda}{
Supplied lambda values.
}
  \item{type.measure}{
Mse type, mean squared error ("mse") or mean absolute error ("mae").
}
  \item{weights}{
A vector of observation weights.
}
  \item{nfolds}{
Number of folds for cross validation.
}
  \item{foldid}{
As in glmnet, an optional vector of values between 1 and nfold identifying what fold each observation is in.
}
  \item{keep}{
Return matrix of predicted responses.
}
  \item{\dots}{Other arguments that can be passed to compLasso.
}
}
\value{
  \item{lambda}{Range of lambda values.}
  \item{cvm}{Vector of mean-cross validation error per lambda value.}
  \item{cvsd}{Corresponding estimate of the standard error of cvm.}
  \item{lamin}{Optimal lambda value.}
}
\references{Hussami, N and Tibshirani, R (2013) A Component Lasso.}
\author{FORTRAN code by Jerry Friedman. R interface by Nadine Hussami and Robert Tibshirani}
\note{The FORTRAN code that this function links to was kindly written
  and provided by Jerry Friedman. This function is based on the cv.glmnet function.
}\seealso{ plot.cv.compLasso, compLasso, predict.compLasso}
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
a=cv.compLasso(x, y, clus, foldid=foldid, type.measure = "mse", nfolds = 5) 
print(a$cvm)
}

