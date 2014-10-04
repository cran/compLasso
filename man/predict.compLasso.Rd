\name{predict.compLasso}
\alias{predict.compLasso}
\title{
Make predictions for a component lasso object.
}
\description{
The function predicts fitted values and coefficients from a fitted component lasso object.
}
\usage{
\method{predict}{compLasso}(object, newx, s = NULL, type = c("response", "coefficients",
"nonzero"), exact = FALSE, ...)
}
\arguments{
  \item{object}{
Fitted "compLasso"" object.
}
  \item{newx}{
A matrix of new x values at which response predictions are to be made. This argument is only required for the type "response".
}
  \item{s}{
Lambda values at which the predictions are to be made.
}
  \item{type}{
The type "response" gives the fitted values. The type "coefficients" returns the coefficient vectors at the new desired s values. The type "nonzero" returns the indices of the non zero coefficients at the different penalty parameter values.
}
  \item{exact}{
If TRUE and there are new s values that were not in the fitted object, the object is refit before making any new predictions. If FALSE, the function uses linear interpolation to predict for values of s that were not in the original fitted object.
}
  \item{\dots}{
Other arguments that can be passed to compLasso.
}
}
\value{
 \item{nbeta}{Coefficients at the required s values if the type is set to "coefficients".}
  \item{nfit}{Predicted response values if the type is set to "response".}
}
\references{
Hussami, N and Tibshirani, R (2013) A Component Lasso.
}
\author{
FORTRAN code by Jerry Friedman. R interface by Nadine Hussami and Robert Tibshirani.
}
\note{
The FORTRAN code that this function links to was kindly written
  and provided by Jerry Friedman. This function is based on the predict.glmnet function.
}

\seealso{
compLasso, cv.compLasso
}
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
print(a)

predict(a,newx=x[1:10,],s=c(0.1,0.02))

}