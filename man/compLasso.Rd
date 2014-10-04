\name{compLasso}
\alias{compLasso}
\title{
Fit a linear model with the component lasso method.
}
\description{
Fits a linear model via the component lasso method given an initial clustering of the variables and an elastic net penalty to use when fitting each component. 
	 
Suppose you expect the covariance matrix of your data to be block diagonal. Then the component lasso method exploits that structure by running a separate regression on each cluster of correlated variables before combining the results into a final coefficient vector. This can be done by running compLasso(x,y,clus).
}
\usage{
compLasso(x, y, clus, weights, parm = 1, jd = rep(0, ncol(x)), 
penalty.factor = rep(1, ncol(x)), lambda = NULL,
lambda.min.ratio = ifelse(nobs < ncol(x), 0.01, 1e-04), standardize = TRUE, 
intercept = TRUE, thresh = 1e-07, nlambda = 100, dfmax = ncol(x) + 1, 
pmax = min(dfmax * 2 + 20, ncol(x)), exclude = NULL, lower.limits = NULL,
upper.limits = NULL, maxit = 1e+05, return.coef = TRUE, nnls = TRUE)
}
\arguments{
  \item{x}{
Input matrix containing an observation per row.
}
  \item{y}{
Response variable.
}
  \item{clus}{
A vector of cluster labels identifying the cluster each variable is in.
}
  \item{weights}{
A vector of observation weights.
}
  \item{parm}{
Elastic Net parameter (1=lasso,0=ridge).
}
  \item{jd}{
Predictor variable deletion flag (jd(1) = 0  => use all variables, jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1) )
}
  \item{penalty.factor}{
Relative penalties for each predictor variable

       penalty.factor(j) = 0 => jth variable unpenalized.
}
  \item{lambda}{
Supplied lambda values.
}
  \item{lambda.min.ratio}{
Smallest value of lambda as a fraction of the smallest lambda for which the returned coefficient vector is zero.
}
  \item{standardize}{
Flag to standardize the x data.
}
  \item{intercept}{
Flag to use an intercept when fitting the model.
}
  \item{thresh}{
Convergence threshold.
}
  \item{nlambda}{
Number of lambda values to use.
}
  \item{dfmax}{
Maximum number of variables allowed to enter largest model.
}
  \item{pmax}{
Maximum number of variables allowed to enter all models along the path.
}
  \item{exclude}{
Which predictors to leave out.
}
  \item{lower.limits}{
Vector of lower range limits for coefficients.
}
  \item{upper.limits}{
Vector of upper range limits for coefficients.
}
  \item{maxit}{
Maximum number of iterations.
}
  \item{return.coef}{
Should the coefficients be returned?
}
  \item{nnls}{
Shoud non-negative least squares be appplied to the final coefficient estimates?
}
}
\value{
  \item{beta}{Matrix of coefficients.}
  \item{a0}{Intercept seqeunce.}
  \item{lambda}{List of lambda values.}
  \item{rss}{R squared values for each solution.}
  \item{npasses}{Actual number of passes over the data for all lambda values.}
}
\references{
Hussami, N and Tibshirani, R (2013) A Component Lasso.
}
\author{
FORTRAN code by Jerry Friedman. R interface by Nadine Hussami and Robert Tibshirani.
}
\note{
The FORTRAN code that this function links to was kindly written
  and provided by Jerry Friedman.
}
\seealso{
cv.compLasso, predict.compLasso
}
\examples{
set.seed(1000)

# Suppose we want to run the component lasso method as in Section 3.2 
# of Hussami and Tibshirani (2013), the algorithm steps are:
# 1. Estimate the connected components by running hierarchical clustering on the
# sample covariance matrix
# 2. Find \beta by running a separate regression on each component
# 3. Re-scale \beta by using a non-negative least squares fit of the response on the responses
# predicted by each component


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

}
