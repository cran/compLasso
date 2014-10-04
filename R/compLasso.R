compLasso=
  function(x,y,clus,weights,parm=1,jd=rep(0,ncol(x)),
           penalty.factor=rep(1,ncol(x)), lambda=NULL, lambda.min.ratio=ifelse(nobs<ncol(x),1e-2,1e-4),
           standardize=TRUE,intercept=TRUE,thresh=1e-7,nlambda=100,
           dfmax=ncol(x)+1,pmax=min(dfmax*2+20,ncol(x)), exclude=NULL, lower.limits=NULL,
           upper.limits=NULL,
           maxit=100000,return.coef=TRUE, nnls=TRUE)
  {
    # inputs x,y,
    # clus- cluster labels
    # obs weights (default= all 1)
    # parm= EN parameter (1=lasso)
    #  jd(jd(1)+1) = predictor variable deletion flag
    #      jd(1) = 0  => use all variables
    #      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
    # penalty.factor  relative penalties for each predictor variable
    #     penalty.factor(j) = 0 => jth variable unpenalized
    # lambda- supplied lambda values
    #  nlambda- number of lambda values to use
    #  dfmax=maximum number of variables allowed to enter largest model
    #        (stopping criterion)
    # pmax= maximum number of variables allowed to enter all models
    #        along path
    # exclude= which predictors to leave out
    # lower.limits and upper.limits= vector of range limits for coefs
    # maxit= max # of iters
    # return.coef= should coefs be returned?
    # nnls= should nonneg LS be applied to final estimates?
    
    
    BIG= 9.9e+35

    nobs=nrow(x)
    nvars=ncol(x)
    nlam=as.integer(nlambda)

    if(is.null(lower.limits)){ lower.limits=rep(-Inf,nvars)}
    if(is.null(upper.limits)){ upper.limits=rep(Inf,nvars)}
    lower.limits[lower.limits==-Inf]=-BIG
    upper.limits[upper.limits==Inf]=BIG
    cl=rbind(lower.limits,upper.limits)



    if(is.null(lambda)){
      if(lambda.min.ratio>=1)stop("lambda.min.ratio should be less than 1")
      flmin=as.double(lambda.min.ratio)
      ulam=double(nlam)
    }
    else{
      flmin=as.double(1)
      if(any(lambda<0))stop("lambdas should be non-negative")
      ulam=as.double(rev(sort(lambda)))
      nlam=as.integer(length(lambda))
    }

    if(!missing(exclude)){
      jd=match(exclude,seq(nvars),0)
      if(!all(jd>0))stop("Some excluded variables out of range")
      jd=as.integer(c(length(jd),jd))
    }
    else jd=as.integer(0)

    isd=as.integer(standardize)
    intr=as.integer(FALSE)
    if(intercept){
      xxx=x
      yyy=y
      x=scale(x,T,F)
      y=y-mean(y)
    }


    if(missing(weights)) weights=rep(1,nobs)


    mode(parm)="double"
    mode(nobs)="integer"
    mode(nvars)="integer"
    mode(x)="double"
    mode(y)="double"
    mode(clus)="integer"
    w=as.double(weights)
    mode(jd)="integer"
    mode(cl)="double"
    ne=as.integer(dfmax)
    nx=as.integer(pmax)
    mode(parm)="double"
    mode(nlam)="integer"
    mode(maxit)="integer"
    thr=as.double(thresh)
    vp=as.double(penalty.factor)
   

    junk=.Fortran("clustelNet",
                  parm,
                  nobs,
                  nvars,
                  x,
                  y,
                  w,
                  clus,
                  jd,
                  vp,
                  cl,
                  ne,
                  nx,
                  nlam,
                  flmin,
                  ulam,
                  thr,
                  isd,
                  intr,
                  maxit,
                  lmu=integer(1),
                  a0=double(nlam),
                  ca=double(nx*nlam),
                  ia=integer(nx),
                  nin=integer(nlam),
                  rsq=double(nlam),
                  alm=double(nlam),
                  nlp=integer(1),
                  jerr=integer(1)
    )

    lmu=junk$lmu
    beta=matrix(NA,nrow=lmu,ncol=ncol(x))
    lamlist=junk$alm[1:lmu]
    
    if(!intercept & !standardize)
     lamlist[1]=max(abs(t(x)%*%y))/(length(y)*parm)
    if(intercept & !standardize)
     lamlist[1]=max(abs(t(x)%*%(y-mean(y))))/(length(y)*parm)
    if( standardize){
      for(i in 1:ncol(x)){
        vc=sum(x[,i]^2)/nrow(x)-(sum(x[,i])/nrow(x))^2
        x[,i]=x[,i]/sqrt(vc)
      }
     if(!intercept)   
      lamlist[1]=max(abs(t(x)%*%y))/(length(y)*parm)
     if(intercept)
      lamlist[1]=max(abs(t(x)%*%(y-mean(y))))/(length(y)*parm)
    }

    if(return.coef){
      nin=junk$nin
      ca=matrix(junk$ca,nrow=nx,ncol=nlam)
      ia=junk$ia
      ni=ncol(x)

      mode(lmu)="integer"
      mode(ca)="double"
      mode(ia)="integer"
      mode(nin)="integer"
      mode(ni)="integer"
      mode(nx)="integer"

      junk2=.Fortran("solns",
                     ni,
                     nx,
                     lmu,
                     ca,
                     ia,
                     nin,
                     beta=double(ni*lmu))

      beta=t(matrix(junk2$beta,ncol=ni,byrow=T))
      dimnames(beta)=list(paste("V",as.character(1:ncol(x)),sep=""),
                          paste("s",as.character(1:ncol(beta)),sep=""))
    }

    if(nnls){
      nc=length(table(clus))
      nlam=length(lamlist)
      xx=matrix(NA,nrow=length(y),ncol=nc)
      bb0=beta
      
      for(k in 1:nlam){
        if( abs(max(beta[,k])) > 1e-8 ){
        for(ii in 1:nc){
          xx[,ii]=x[,clus==ii,drop=F]%*%beta[clus==ii,k,drop=F]
        }
        aaa=nnls(xx,y)
        for(ii in 1:nc){
          beta[clus==ii,k]=beta[clus==ii,k,drop=F]*aaa[ii]
        }
       }
      }
    }

    a0=rep(0,lmu)
    if(intercept){
      a0=rep(mean(yyy), lmu)-colMeans(xxx, na.rm = FALSE, dims = 1)%*%beta
    }

    result <- list(a0=a0, beta=beta,lambda=lamlist,rss=junk$rsq[seq(lmu)],
                   npasses=junk$nlp,jerr=junk$jerr)
    class(result) <- "compLasso"
    result
}


nnls=
  function(x,y,SMALL=1e-8){
  o=colSums(abs(x))>SMALL
  pp=sum(o)
  Dmat=t(x[,o])%*%x[,o]
  dvec=t(x[,o])%*%y
  Amat=diag(pp)
  bvec=rep(0,pp)
  a=solve.QP(Dmat, dvec, Amat, bvec)
  bb=rep(0,ncol(x))
  bb[o]=a$sol
  return(bb)
}

# this function is based on the cv.glmnet function
cv.compLasso=
  function (x, y, clus, lambda = NULL, type.measure = c("mse","mae"), weights, nfolds = 10, foldid,
             keep = FALSE, ...)
  {
    
    # inputs x,y,
    # clus- cluster labels
    # obs weights (default= all 1)
    # type.measure mse type, mean squared error or mean absolute error
    # nfolds number of folds for cross validation
    # foldid an optional vector of values between 1 and nfold identifying what fold each observation is in.
    # keep to return the matrix of predicted responses
    # rest of component lasso arguments

    if (missing(type.measure))
      type.measure = "default"
    else type.measure = match.arg(type.measure)

    if (!is.null(lambda) && length(lambda) < 2)
      stop("Need more than one value of lambda for cv.compLasso")

    N = nrow(x)

    if (missing(weights))
      weights = rep(1, N)
    else weights = as.double(weights)

    y = drop(y)
    
    compLasso.object = compLasso(x, y, clus, lambda=lambda, ...)# if lambda null generates its own sequence 
    lambda = compLasso.object$lambda
    
    nz = sapply(predict(compLasso.object, type = "nonzero"),
                length)

    if (missing(foldid))
      foldid = sample(rep(seq(nfolds), length = N))
    else nfolds = max(foldid)
    if (nfolds < 3)
      stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist = as.list(seq(nfolds))


  for (i in seq(nfolds)) {
    which = foldid == i
    y_sub = y[!which]
    
    outlist[[i]] = compLasso(x[!which, , drop = FALSE],
                             y_sub, clus, lambda = lambda, weights=weights[!which],
                             ...)
  }
 

cvstuff = cv.elnet(outlist, lambda, x, y, weights, foldid, type.measure, keep)
cvm = cvstuff$cvm
cvsd = cvstuff$cvsd
cvname = cvstuff$name

out = list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvup = cvm +
             cvsd, cvlo = cvm - cvsd, nzero = nz, name = cvname, compLasso.fit = compLasso.object)

if (keep)
  out = c(out, list(fit.preval = cvstuff$fit.preval, foldid = foldid))
lamin = getmin(lambda, cvm, cvsd)
obj = c(out, as.list(lamin))
class(obj) = "cv.compLasso"
obj

}

# this function is based on the glmnet.predict function
predict.compLasso=
  function(object, newx, s = NULL, type = c("response", "coefficients", "nonzero"), exact = FALSE,
           ...)
  {
    type = match.arg(type)
    if (missing(newx)) {
      if (!match(type, c("coefficients", "nonzero"), FALSE))
        stop("You need to supply a value for 'newx'")
    }

    if (exact && (!is.null(s))) {
      lambda = object$lambda
      which = match(s, lambda, FALSE)
      if (!all(which > 0)) {
        lambda = unique(rev(sort(c(s, lambda))))
        object = update(object, lambda = lambda)
      }
    }

    a0 = matrix(object$a0, nrow=1, ncol=length(object$a0), byrow=TRUE)
    rownames(a0) = "(Intercept)"
    nbeta = rbind2(a0, object$beta)

    if (!is.null(s)) {
      vnames = dimnames(nbeta)[[1]]
      dimnames(nbeta) = list(NULL, NULL)
      lambda = object$lambda
      lamlist = lambda.interp(lambda, s)
      nbeta = nbeta[, lamlist$left, drop = FALSE] %*% Diagonal(x = lamlist$frac) +
        nbeta[, lamlist$right, drop = FALSE] %*% Diagonal(x = 1 -
                                                            lamlist$frac)
      dimnames(nbeta) = list(vnames, paste(seq(along = s)))
    }

    if (type == "coefficients")
      return(nbeta)
    if (type == "nonzero")
      return(nonzeroCoef(nbeta[-1, , drop = FALSE], bystep = TRUE))
    if (inherits(newx, "sparseMatrix"))
      newx = as(newx, "dgCMatrix")
    nfit = as.matrix(cbind2(1, newx) %*% nbeta)

    nfit
  }


plot.compLasso=
  function (x, xvar = c("norm", "lambda"), label = FALSE,
            ...)
  {
    xvar = match.arg(xvar)
    df = apply(x$beta, 2, function(c)sum(abs(c)>1e-8))
    plotCoef(x$beta, lambda = x$lambda, df=df,
             label = label, xvar = xvar, ...)
  }


# function taken from the glmnet package
plot.cv.compLasso=
  function (x, sign.lambda = 1, ...)
  {
    cvobj = x
    xlab = "log(Lambda)"
    if (sign.lambda < 0)
      xlab = paste("-", xlab, sep = "")
    plot.args = list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm,
                     ylim = range(cvobj$cvup, cvobj$cvlo), xlab = xlab, ylab = cvobj$name,
                     type = "n")
    new.args = list(...)
    if (length(new.args))
      plot.args[names(new.args)] = new.args
    do.call("plot", plot.args)
    error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvup, cvobj$cvlo,
               width = 0.01, col = "darkgrey")
    points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20,
           col = "red")
    axis(side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz),
         tick = FALSE, line = 0)
    abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
    abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
    invisible()
}

# The following functions are all internal functions taken from the glmnet package

cv.elnet=
  function(outlist, lambda, x, y, weights, foldid, type.measure,
            keep = FALSE)
  {

    typenames = c( mse = "Mean-Squared Error",
                   mae = "Mean Absolute Error")
    if (type.measure == "default")
      type.measure = "mse"
    if (!match(type.measure, c("mse", "mae"), FALSE)) {
      warning("Only 'mse', or 'mae'  available for Gaussian models; 'mse' used")
      type.measure = "mse"
    }

    predmat = matrix(NA, length(y), length(lambda))
    nfolds = max(foldid)
    nlams = double(nfolds)

    for (i in seq(nfolds)) {
      which = foldid == i
      fitobj = outlist[[i]]
      preds = predict(fitobj, x[which, , drop = FALSE])
      nlami = length(outlist[[i]]$lambda)
      predmat[which, seq(nlami)] = preds
      nlams[i] = nlami
    }

    N = length(y) - apply(is.na(predmat), 2, sum)

    cvraw = switch(type.measure, mse = (y - predmat)^2, mae = abs(y - predmat))

    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                      w = weights, na.rm = TRUE)/(N - 1))
    out = list(cvm = cvm, cvsd = cvsd, name = typenames[type.measure])
    if (keep)
      out$fit.preval = predmat
    out


  }

cvcompute=
  function (mat, weights, foldid, nlams)
  {
    wisum = tapply(weights, foldid, sum)
    nfolds = max(foldid)
    outmat = matrix(NA, nfolds, ncol(mat))
    good = matrix(0, nfolds, ncol(mat))
    mat[is.infinite(mat)] = NA

    for (i in seq(nfolds)) {
      mati = mat[foldid == i, , drop = FALSE]
      wi = weights[foldid == i]
      outmat[i, ] = apply(mati, 2, weighted.mean, w = wi, na.rm = TRUE)
      good[i, seq(nlams[i])] = 1
    }

    N = apply(good, 2, sum)
    list(cvraw = outmat, weights = wisum, N = N)
  }


error.bars=
  function (x, upper, lower, width = 0.02, ...)
  {
    xlim <- range(x)
    barw <- diff(xlim) * width
    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
    range(upper, lower)
  }


plotCoef=
  function (beta, norm, lambda, df, label = FALSE, xvar = c("norm","lambda"),xlab = iname, ylab = "Coefficients", ...)
{
  which = nonzeroCoef(beta)
  nwhich = length(which)
  switch(nwhich + 1, `0` = {
    warning("No plot produced since all coefficients zero")
    return()
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  beta = as.matrix(beta[which, , drop = FALSE])
  xvar = match.arg(xvar)
  switch(xvar, norm = {
    index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(lambda)
    iname = "Log Lambda"
    approx.f = 0
  })
  dotlist = list(...)
  type = dotlist$type
  if (is.null(type))
    matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
            type = "l", ...)
  else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab,
               ...)
  atdf = pretty(index)
  prettydf = approx(x = index, y = df, xout = atdf, rule = 2,
                    method = "constant", f = approx.f)$y
  axis(3, at = atdf, labels = prettydf, tcl = NA)
  if (label) {
    nnz = length(which)
    xpos = max(index)
    pos = 4
    if (xvar == "lambda") {
      xpos = min(index)
      pos = 2
    }
    xpos = rep(xpos, nnz)
    ypos = beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
}

nonzeroCoef=
  function(beta, bystep = FALSE)
{
    nr = nrow(beta)
    if (nr == 1) {
      if (bystep)
        apply(beta, 2, function(x) if (abs(x) > 0)
          1
          else NULL)
      else {
        if (any(abs(beta) > 0))
          1
        else NULL
      }
    }
    else {
      beta = abs(beta) > 0
      which = seq(nr)
      ones = rep(1, ncol(beta))
      nz = as.vector((beta %*% ones) > 0)
      which = which[nz]
      if (bystep) {
        if (length(which) > 0) {
          beta = as.matrix(beta[which, , drop = FALSE])
          nzel = function(x, which) if (any(x))
            which[x]
          else NULL
          which = apply(beta, 2, nzel, which)
          if (!is.list(which))
            which = data.frame(which)
          which
        }
        else {
          dn = dimnames(beta)[[2]]
          which = vector("list", length(dn))
          names(which) = dn
          which
        }
      }
      else which
    }
}


getmin=
  function (lambda, cvm, cvsd)
{
  cvmin = min(cvm, na.rm = TRUE)
  idmin = cvm <= cvmin
  lambda.min = max(lambda[idmin], na.rm = TRUE)
  idmin = match(lambda.min, lambda)
  semin = (cvm + cvsd)[idmin]
  idmin = cvm <= semin
  lambda.1se = max(lambda[idmin], na.rm = TRUE)
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

lambda.interp=function (lambda, s)
{
    if (length(lambda) == 1) {
        nums = length(s)
        left = rep(1, nums)
        right = left
        sfrac = rep(1, nums)
    }
    else {
        s[s > max(lambda)] = max(lambda)
        s[s < min(lambda)] = min(lambda)
        k = length(lambda)
        sfrac <- (lambda[1] - s)/(lambda[1] - lambda[k])
        lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
        coord <- approx(lambda, seq(lambda), sfrac)$y
        left <- floor(coord)
        right <- ceiling(coord)
        sfrac = (sfrac - lambda[right])/(lambda[left] - lambda[right])
        sfrac[left == right] = 1
    }
    list(left = left, right = right, frac = sfrac)
}


