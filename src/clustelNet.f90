c
c                          ClusterNet (12/28/13)
c
c
c              Clustered elastic net with squared-error loss
c
c
c call clustelnet(parm,no,ni,x,y,w,ke,jd,vp,cl,ne,nx,nlam,flmin,ulam,thr,isd,
c            intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
c
c input:
c
c   parm = penalty member index (0 <= parm <= 1)
c        = 0.0 => ridge
c        = 1.0 => lasso
c   no = number of observations
c   ni = number of predictor variables
c   x(no,ni) = predictor data matrix flat file (overwritten)
c   y(no) = response vector
c   w(no)= observation weights (overwritten)
c   ke(ni) = predictor variable cluster encoder
c      ke(j)=k => variable j in cluster k
c   jd(jd(1)+1) = predictor variable deletion flag
c      jd(1) = 0  => use all variables
c      jd(1) != 0 => do not use variables jd(2)...jd(jd(1)+1)
c   vp(ni) = relative penalties for each predictor variable
c      vp(j) = 0 => jth variable unpenalized
c   cl(2,ni) = interval constraints on coefficient values (overwritten)
c      cl(1,j) = lower bound for jth coefficient value (<= 0.0)
c      cl(2,j) = upper bound for jth coefficient value (>= 0.0)
c   ne = maximum number of variables allowed to enter largest model
c        (stopping criterion)
c   nx = maximum number of variables allowed to enter all models
c        along path (memory allocation, nx > ne).
c   nlam = (maximum) number of lamda values
c   flmin = user control of lamda values (>=0)
c      flmin < 1.0 => minimum lamda = flmin*(largest lamda value)
c      flmin >= 1.0 => use supplied lamda values (see below)
c   ulam(nlam) = user supplied lamda values (ignored if flmin < 1.0)
c   thr = convergence threshold for each lamda solution.
c      iterations stop when the maximum reduction in the criterion value
c      as a result of each parameter update over a single pass
c      is less than thr times the null criterion value.
c      (suggested value, thr=1.0e-5)
c   isd = predictor variable standarization flag:
c      isd = 0 => regression on original predictor variables
c      isd = 1 => regression on standardized predictor variables
c      Note: output solutions always reference original
c            variables locations and scales.
c   intr = intercept flag
c      intr = 0/1 => don't/do include intercept in model
c   maxit = maximum allowed number of passes over the data for all lambda
c      values (suggested values, maxit = 100000)
c
c output:
c
c   lmu = actual number of lamda values (solutions)
c   a0(lmu) = intercept values for each solution
c   ca(nx,lmu) = compressed coefficient values for each solution
c   ia(nx) = pointers to compressed coefficients
c   nin(lmu) = number of compressed coefficients for each solution
c   rsq(lmu) = R**2 values for each solution
c   alm(lmu) = lamda values corresponding to each solution
c   nlp = actual number of passes over the data for all lamda values
c   jerr = error flag:
c      jerr  = 0 => no error
c      jerr > 0 => fatal error - no output returned
c         jerr < 7777 => memory allocation error
c         jerr = 7777 => all used predictors have zero variance
c         jerr = 10000 => maxval(vp) <= 0.0
C      jerr < 0 => non fatal error - partial output:
c         Solutions for larger lamdas (1:(k-1)) returned.
c         jerr = -k => convergence for kth lamda value not reached
c            after maxit (see above) iterations.
c         jerr = -10000-k => number of non zero coefficients along path
c            exceeds nx (see above) at kth lamda value.
c
c
c
c least-squares utility routines:
c
c
c uncompress coefficient vectors for all solutions:
c
c call solns(ni,nx,lmu,ca,ia,nin,b)
c
c input:
c
c    ni,nx = input to clustelnet
c    lmu,ca,ia,nin = output from clustelnet
c
c output:
c
c    b(ni,lmu) = all elnet returned solutions in uncompressed format
c
c
c uncompress coefficient vector for particular solution:
c
c call uncomp(ni,ca,ia,nin,a)
c
c input:
c
c    ni = total number of predictor variables
c    ca(nx) = compressed coefficient values for the solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for the solution
c
c output:
c
c    a(ni) =  uncompressed coefficient vector
c             referencing original variables
c
c
c evaluate linear model from compressed coefficients and
c uncompressed predictor matrix:
c
c call modval(a0,ca,ia,nin,n,x,f);
c
c input:
c
c    a0 = intercept
c    ca(nx) = compressed coefficient values for a solution
c    ia(nx) = pointers to compressed coefficients
c    nin = number of compressed coefficients for solution
c    n = number of predictor vectors (observations)
c    x(n,ni) = full (uncompressed) predictor matrix
c
c output:
c
c    f(n) = model predictions
c
c
c
c
c                Changing internal parameter values
c
c
c call chg_fract_dev(fdev)
c   fdev = minimum fractional change in deviance for stopping path
c      default = 1.0e-5
c
c call chg_dev_max(devmax)
c   devmax = maximum fraction of explained deviance for stopping path
c      default = 0.999
c
c call chg_min_flmin(eps)
c   eps = minimum value of flmin (see above). default= 1.0e-6
c
c call chg_big(big)
c   big = large floating point number. default = 9.9e35
c
c call chg_min_lambdas(mnlam)
c   mnlam = minimum number of path points (lambda values) allowed
c      default = 5
c
c
c             Obtain current internal parameter values
cg
c call get_int_parms(fdev,eps,big,mnlam,devmax)
c
c
c             
      subroutine get_int_parms(sml,eps,big,mnlam,rsqmax)                    250
      data sml0,eps0,big0,mnlam0,rsqmax0  /1.0e-5,1.0e-6,9.9e35,5,0.999/    252
      sml=sml0                                                              252
      eps=eps0                                                              252
      big=big0                                                              252
      mnlam=mnlam0                                                          252
      rsqmax=rsqmax0                                                        253
      return                                                                254
      entry chg_fract_dev(arg)                                              254
      sml0=arg                                                              254
      return                                                                255
      entry chg_dev_max(arg)                                                255
      rsqmax0=arg                                                           255
      return                                                                256
      entry chg_min_flmin(arg)                                              256
      eps0=arg                                                              256
      return                                                                257
      entry chg_big(arg)                                                    257
      big0=arg                                                              257
      return                                                                258
      entry chg_min_lambdas(irg)                                            258
      mnlam0=irg                                                            258
      return                                                                259
      end                                                                   260
      subroutine clustelnet  (parm,no,ni,x,y,w,ke,jd,vp,cl,ne,nx,nlam,fl    263 
     *min,ulam,thr,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real x(no,ni),y(no),w(no),vp(ni),ca(nx,nlam),cl(2,ni)                 264
      real ulam(nlam),a0(nlam),rsq(nlam),alm(nlam)                          265
      integer jd(*),ia(nx),nin(nlam),ke(ni)                                 266
      real, dimension (:), allocatable :: vq;                                   
      if(maxval(vp) .gt. 0.0)goto 10021                                     269
      jerr=10000                                                            269
      return                                                                269
10021 continue                                                              270
      allocate(vq(1:ni),stat=jerr)                                          270
      if(jerr.ne.0) return                                                  271
      vq=max(0.0,vp)                                                        271
      vq=vq*ni/sum(vq)                                                      272
      call elnetn (parm,no,ni,x,y,w,ke,jd,vq,cl,ne,nx,nlam,flmin,ulam,th    275 
     *r,isd,intr,maxit,  lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      deallocate(vq)                                                        276
      return                                                                277
      end                                                                   278
      subroutine elnetn (parm,no,ni,x,y,w,ke,jd,vp,cl,ne,nx,nlam,flmin,u    280 
     *lam,thr,  isd,intr,maxit,lmu,a0,ca,ia,nin,rsq,alm,nlp,jerr)
      real vp(ni),x(no,ni),y(no),w(no),ulam(nlam),cl(2,ni)                  281
      real ca(nx,nlam),a0(nlam),rsq(nlam),alm(nlam)                         282
      integer jd(*),ia(nx),nin(nlam),ke(ni)                                 283
      real, dimension (:), allocatable :: xm,xs,xv,vlam,fe                      
      integer, dimension (:), allocatable :: ju,kp,kq                           
      allocate(xm(1:ni),stat=jerr)                                          288
      allocate(xs(1:ni),stat=ierr)                                          288
      jerr=jerr+ierr                                                        289
      allocate(ju(1:ni),stat=ierr)                                          289
      jerr=jerr+ierr                                                        290
      allocate(xv(1:ni),stat=ierr)                                          290
      jerr=jerr+ierr                                                        291
      allocate(vlam(1:nlam),stat=ierr)                                      291
      jerr=jerr+ierr                                                        292
      allocate(kq(1:ni),stat=ierr)                                          292
      jerr=jerr+ierr                                                        293
      allocate(kp(1:(ni+1)),stat=ierr)                                      293
      jerr=jerr+ierr                                                        294
      allocate(fe(1:ni),stat=ierr)                                          294
      jerr=jerr+ierr                                                        295
      if(jerr.ne.0) return                                                  296
      call chkvars(no,ni,x,ju)                                              297
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0                                  298
      if(maxval(ju) .gt. 0)goto 10041                                       298
      jerr=7777                                                             298
      return                                                                298
10041 continue                                                              299
      call standard1(no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)           300
      if(jerr.ne.0) return                                                  301
      cl=cl/ys                                                              301
      if(isd .le. 0)goto 10061                                              301
10070 do 10071 j=1,ni                                                       301
      cl(:,j)=cl(:,j)*xs(j)                                                 301
10071 continue                                                              301
10072 continue                                                              301
10061 continue                                                              302
      if(flmin.ge.1.0) vlam=ulam/ys                                         303
10080 do 10081 j=1,ni                                                       303
      fe(j)=ke(j)                                                           303
      kq(j)=j                                                               303
10081 continue                                                              303
10082 continue                                                              303
      call psort7(fe,kq,1,ni)                                               304
      deallocate(fe)                                                        305
      call spcol(ni,ke,kq,kl,kp)                                            305
      kk=kl-1                                                               306
      call elnet2(parm,ni,kk,kp,kq,ju,vp,cl,y,no,ne,nx,x,nlam,flmin,vlam    308 
     *,thr,  maxit,xv,lmu,ca,ia,nin,rsq,alm,nlp,jerr)
      if(jerr.gt.0) return                                                  309
      deallocate(kp,kq,vlam,ju)                                             310
10090 do 10091 k=1,lmu                                                      310
      alm(k)=ys*alm(k)                                                      310
      nk=nin(k)                                                             311
10100 do 10101 l=1,nk                                                       311
      ca(l,k)=ys*ca(l,k)/xs(ia(l))                                          311
10101 continue                                                              311
10102 continue                                                              311
      a0(k)=0.0                                                             312
      if(intr.ne.0) a0(k)=ym-dot_product(ca(1:nk,k),xm(ia(1:nk)))           313
10091 continue                                                              314
10092 continue                                                              314
      deallocate(xm,xs,xv)                                                  315
      return                                                                316
      end                                                                   317
      subroutine standard1 (no,ni,x,y,w,isd,intr,ju,xm,xs,ym,ys,xv,jerr)    318
      real x(no,ni),y(no),w(no),xm(ni),xs(ni),xv(ni)                        318
      integer ju(ni)                                                        319
      real, dimension (:), allocatable :: v                                     
      allocate(v(1:no),stat=jerr)                                           322
      if(jerr.ne.0) return                                                  323
      w=w/sum(w)                                                            323
      v=sqrt(w)                                                             324
      if(intr .ne. 0)goto 10121                                             324
      ym=0.0                                                                324
      y=v*y                                                                 325
      ys=sqrt(dot_product(y,y)-dot_product(v,y)**2)                         325
      y=y/ys                                                                326
10130 do 10131 j=1,ni                                                       326
      if(ju(j).eq.0)goto 10131                                              326
      xm(j)=0.0                                                             326
      x(:,j)=v*x(:,j)                                                       327
      xv(j)=dot_product(x(:,j),x(:,j))                                      328
      if(isd .eq. 0)goto 10151                                              328
      xbq=dot_product(v,x(:,j))**2                                          328
      vc=xv(j)-xbq                                                          329
      xs(j)=sqrt(vc)                                                        329
      x(:,j)=x(:,j)/xs(j)                                                   329
      xv(j)=1.0+xbq/vc                                                      330
      goto 10161                                                            331
10151 continue                                                              331
      xs(j)=1.0                                                             331
10161 continue                                                              332
10141 continue                                                              332
10131 continue                                                              333
10132 continue                                                              333
      go to 10170                                                           334
10121 continue                                                              335
10180 do 10181 j=1,ni                                                       335
      if(ju(j).eq.0)goto 10181                                              336
      xm(j)=dot_product(w,x(:,j))                                           336
      x(:,j)=v*(x(:,j)-xm(j))                                               337
      xv(j)=dot_product(x(:,j),x(:,j))                                      337
      if(isd.gt.0) xs(j)=sqrt(xv(j))                                        338
10181 continue                                                              339
10182 continue                                                              339
      if(isd .ne. 0)goto 10201                                              339
      xs=1.0                                                                339
      goto 10211                                                            340
10201 continue                                                              340
10220 do 10221 j=1,ni                                                       340
      if(ju(j).eq.0)goto 10221                                              340
      x(:,j)=x(:,j)/xs(j)                                                   340
10221 continue                                                              341
10222 continue                                                              341
      xv=1.0                                                                342
10211 continue                                                              343
10191 continue                                                              343
      ym=dot_product(w,y)                                                   343
      y=v*(y-ym)                                                            343
      ys=sqrt(dot_product(y,y))                                             343
      y=y/ys                                                                344
10170 continue                                                              344
      deallocate(v)                                                         345
      return                                                                346
      end                                                                   347
      subroutine spcol(n,jx,m,j,ia)                                         348
      integer jx(n),m(n),ia(n+1)                                            349
      j=1                                                                   349
      jx0=jx(m(j))                                                          349
      ia(j)=1                                                               350
10230 do 10231 i=2,n                                                        351
      if(jx(m(i)) .eq. jx0)goto 10251                                       351
      j=j+1                                                                 351
      ia(j)=i                                                               351
10251 continue                                                              352
      jx0=jx(m(i))                                                          353
10231 continue                                                              354
10232 continue                                                              354
      j=j+1                                                                 354
      ia(j)=n+1                                                             355
      return                                                                356
      end                                                                   357
      subroutine elnet2(beta,ni,kk,kp,kq,ju,vp,cl,y,no,ne,nx,x,nlam,flmi    359 
     *n,  ulam,thr,maxit,xv,lmu,ao,ia,kin,rsqo,almo,nlp,jerr)
      real vp(ni),y(no),x(no,ni),ulam(nlam),ao(nx,nlam),rsqo(nlam),almo(    360 
     *nlam),xv(ni)
      real cl(2,ni)                                                         361
      integer ju(ni),ia(nx),kin(nlam),kp(kk+1),kq(ni)                       362
      real, dimension (:), allocatable :: a,g                                   
      integer, dimension (:), allocatable :: mm,ix,kch                          
      real, dimension(:,:), allocatable :: r                                    
      allocate(r(1:no,1:kk),stat=jerr)                                          
      call get_int_parms(sml,eps,big,mnlam,rsqmax)                          369
      allocate(a(1:ni),stat=ierr)                                           369
      jerr=jerr+ierr                                                        370
      allocate(mm(1:ni),stat=ierr)                                          370
      jerr=jerr+ierr                                                        371
      allocate(g(1:ni),stat=ierr)                                           371
      jerr=jerr+ierr                                                        372
      allocate(ix(1:ni),stat=ierr)                                          372
      jerr=jerr+ierr                                                        373
      allocate(kch(1:kk),stat=ierr)                                         373
      jerr=jerr+ierr                                                        374
      if(jerr.ne.0) return                                                  375
      bta=beta                                                              375
      omb=1.0-bta                                                           375
      ix=0                                                                  376
      if(flmin .ge. 1.0)goto 10271                                          376
      eqs=max(eps,flmin)                                                    376
      alf=eqs**(1.0/(nlam-1))                                               376
10271 continue                                                              377
      rsq=0.0                                                               377
      a=0.0                                                                 377
      mm=0                                                                  377
      nlp=0                                                                 377
      nin=nlp                                                               377
      mnl=min(mnlam,nlam)                                                   377
      alm=0.0                                                               377
      kch=0                                                                 378
10280 do 10281 j=1,ni                                                       378
      if(ju(j).eq.0)goto 10281                                              378
      g(j)=dot_product(y,x(:,j))                                            378
10281 continue                                                              379
10282 continue                                                              379
10290 do 10291 kc=1,kk                                                      379
      r(:,kc)=y                                                             379
10291 continue                                                              379
10292 continue                                                              379
      fkk=kk                                                                380
10300 do 10301 m=1,nlam                                                     380
      alm0=alm                                                              381
      if(flmin .lt. 1.0)goto 10321                                          381
      alm=ulam(m)                                                           381
      goto 10311                                                            382
10321 if(m .le. 2)goto 10331                                                382
      alm=alm*alf                                                           382
      goto 10311                                                            383
10331 if(m .ne. 1)goto 10341                                                383
      alm=big                                                               383
      goto 10351                                                            384
10341 continue                                                              384
      alm0=0.0                                                              385
10360 do 10361 j=1,ni                                                       385
      if(ju(j).eq.0)goto 10361                                              386
      if(vp(j).gt.0.0) alm0=max(alm0,abs(g(j))/vp(j))                       387
10361 continue                                                              388
10362 continue                                                              388
      alm0=alm0/max(bta,1.0e-3)                                             388
      alm=alf*alm0                                                          389
10351 continue                                                              390
10311 continue                                                              390
      dem=alm*omb                                                           390
      ab=alm*bta                                                            390
      rsq0=rsq                                                              391
      tlam=bta*(2.0*alm-alm0)                                               392
10370 do 10371 k=1,ni                                                       392
      if(ix(k).eq.1)goto 10371                                              392
      if(ju(k).eq.0)goto 10371                                              393
      if(abs(g(k)).gt.tlam*vp(k)) ix(k)=1                                   394
10371 continue                                                              395
10372 continue                                                              395
10380 do 10381 kc=1,kk                                                      396
10390 continue                                                              396
10391 continue                                                              396
10400 continue                                                              397
      nlp=nlp+1                                                             397
      dlx=0.0                                                               398
10410 do 10411 l=kp(kc),kp(kc+1)-1                                          398
      k=kq(l)                                                               398
      if(ix(k).eq.0)goto 10411                                              399
      if(kch(kc).gt.0) g(k)=dot_product(r(:,kc),x(:,k))                     400
      ak=a(k)                                                               400
      u=g(k)+ak*xv(k)                                                       400
      v=abs(u)-vp(k)*ab                                                     400
      a(k)=0.0                                                              402
      if(v.gt.0.0) a(k)=max(cl(1,k),min(cl(2,k),sign(v,u)/(xv(k)+vp(k)*d    403 
     *em)))
      if(a(k).eq.ak)goto 10411                                              404
      if(mm(k) .ne. 0)goto 10431                                            404
      nin=nin+1                                                             404
      if(nin.gt.nx)goto 10412                                               405
      mm(k)=nin                                                             405
      ia(nin)=k                                                             406
10431 continue                                                              407
      del=a(k)-ak                                                           407
      rsq=rsq+del*(2.0*g(k)-del*xv(k))/fkk                                  408
      r(:,kc)=r(:,kc)-del*x(:,k)                                            408
      kch(kc)=1                                                             409
      dlx=max(xv(k)*del**2,dlx)                                             411
10411 continue                                                              412
10412 continue                                                              412
      if(nin.gt.nx)goto 10392                                               413
      if(dlx .ge. thr)goto 10451                                            413
      ixx=0                                                                 414
10460 do 10461 l=kp(kc),kp(kc+1)-1                                          414
      k=kq(l)                                                               414
      if(ix(k).eq.1)goto 10461                                              415
      if(ju(k).eq.0)goto 10461                                              416
      if(kch(kc).gt.0) g(k)=dot_product(r(:,kc),x(:,k))                     417
      if(abs(g(k)) .le. ab*vp(k))goto 10481                                 417
      ix(k)=1                                                               417
      ixx=1                                                                 417
10481 continue                                                              418
10461 continue                                                              419
10462 continue                                                              419
      if(ixx.eq.1) go to 10400                                              420
      goto 10392                                                            421
10451 continue                                                              422
      if(nlp .le. maxit)goto 10501                                          422
      jerr=-m                                                               422
      return                                                                422
10501 continue                                                              423
      goto 10391                                                            424
10392 continue                                                              424
10381 continue                                                              425
10382 continue                                                              425
      if(nin .le. nx)goto 10521                                             425
      jerr=-10000-m                                                         425
      goto 10302                                                            425
10521 continue                                                              426
      if(nin.gt.0) ao(1:nin,m)=a(ia(1:nin))                                 426
      kin(m)=nin                                                            427
      rsqo(m)=rsq                                                           427
      almo(m)=alm                                                           427
      lmu=m                                                                 428
      if(m.lt.mnl)goto 10301                                                428
      if(flmin.ge.1.0)goto 10301                                            429
      me=0                                                                  429
10530 do 10531 j=1,nin                                                      429
      if(ao(j,m).ne.0.0) me=me+1                                            429
10531 continue                                                              429
10532 continue                                                              429
      if(me.gt.ne)goto 10302                                                430
      if(rsq-rsq0.lt.sml*rsq)goto 10302                                     430
      if(rsq.gt.rsqmax)goto 10302                                           431
10301 continue                                                              432
10302 continue                                                              432
      deallocate(a,mm,g,ix)                                                 433
      return                                                                434
      end                                                                   435
      subroutine chkvars(no,ni,x,ju)                                        436
      real x(no,ni)                                                         436
      integer ju(ni)                                                        437
10540 do 10541 j=1,ni                                                       437
      ju(j)=0                                                               437
      t=x(1,j)                                                              438
10550 do 10551 i=2,no                                                       438
      if(x(i,j).eq.t)goto 10551                                             438
      ju(j)=1                                                               438
      goto 10552                                                            438
10551 continue                                                              439
10552 continue                                                              439
10541 continue                                                              440
10542 continue                                                              440
      return                                                                441
      end                                                                   442
      subroutine solns(ni,nx,lmu,a,ia,nin,b)                                443
      real a(nx,lmu),b(ni,lmu)                                              443
      integer ia(nx),nin(lmu)                                               444
10560 do 10561 lam=1,lmu                                                    444
      call uncomp(ni,a(:,lam),ia,nin(lam),b(:,lam))                         444
10561 continue                                                              445
10562 continue                                                              445
      return                                                                446
      end                                                                   447
      subroutine uncomp(ni,ca,ia,nin,a)                                     448
      real ca(*),a(ni)                                                      448
      integer ia(*)                                                         449
      a=0.0                                                                 449
      if(nin.gt.0) a(ia(1:nin))=ca(1:nin)                                   450
      return                                                                451
      end                                                                   452
      subroutine modval(a0,ca,ia,nin,n,x,f)                                 453
      real ca(nin),x(n,*),f(n)                                              453
      integer ia(nin)                                                       454
      f=a0                                                                  454
      if(nin.le.0) return                                                   455
10570 do 10571 i=1,n                                                        455
      f(i)=f(i)+dot_product(ca(1:nin),x(i,ia(1:nin)))                       455
10571 continue                                                              456
10572 continue                                                              456
      return                                                                457
      end                                                                   458
      subroutine psort7 (v,a,ii,jj)                                             
c                                                                               
c     puts into a the permutation vector which sorts v into                     
c     increasing order. the array v is not modified.                            
c     only elements from ii to jj are considered.                               
c     arrays iu(k) and il(k) permit sorting up to 2**(k+1)-1 elements           
c                                                                               
c     this is a modification of cacm algorithm #347 by r. c. singleton,         
c     which is a modified hoare quicksort.                                      
c                                                                               
      dimension a(jj),v(jj),iu(20),il(20)                                       
      integer t,tt                                                              
      integer a                                                                 
      real v                                                                    
      m=1                                                                       
      i=ii                                                                      
      j=jj                                                                      
 10   if (i.ge.j) go to 80                                                      
 20   k=i                                                                       
      ij=(j+i)/2                                                                
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 30                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
 30   l=j                                                                       
      if (v(a(j)).ge.vt) go to 50                                               
      a(ij)=a(j)                                                                
      a(j)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 50                                               
      a(ij)=a(i)                                                                
      a(i)=t                                                                    
      t=a(ij)                                                                   
      vt=v(t)                                                                   
      go to 50                                                                  
 40   a(l)=a(k)                                                                 
      a(k)=tt                                                                   
 50   l=l-1                                                                     
      if (v(a(l)).gt.vt) go to 50                                               
      tt=a(l)                                                                   
      vtt=v(tt)                                                                 
 60   k=k+1                                                                     
      if (v(a(k)).lt.vt) go to 60                                               
      if (k.le.l) go to 40                                                      
      if (l-i.le.j-k) go to 70                                                  
      il(m)=i                                                                   
      iu(m)=l                                                                   
      i=k                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 70   il(m)=k                                                                   
      iu(m)=j                                                                   
      j=l                                                                       
      m=m+1                                                                     
      go to 90                                                                  
 80   m=m-1                                                                     
      if (m.eq.0) return                                                        
      i=il(m)                                                                   
      j=iu(m)                                                                   
 90   if (j-i.gt.10) go to 20                                                   
      if (i.eq.ii) go to 10                                                     
      i=i-1                                                                     
 100  i=i+1                                                                     
      if (i.eq.j) go to 80                                                      
      t=a(i+1)                                                                  
      vt=v(t)                                                                   
      if (v(a(i)).le.vt) go to 100                                              
      k=i                                                                       
 110  a(k+1)=a(k)                                                               
      k=k-1                                                                     
      if (vt.lt.v(a(k))) go to 110                                              
      a(k+1)=t                                                                  
      go to 100                                                                 
      end                                                                       
