# Kernel and kernel matrix
kstand <- function(x,ker='normal',knorm="sq")
{
  # Densities such as the integral of the square of the density is one if knorm is sq 
  # or such that the sd is one if knorm is sd.
  
  if (ker=='normal') 
  { 
    s <- switch(knorm, sd = 1, sq = 1/(2*sqrt(pi)))
    aux <- exp(- x^2 /(2*s^2)) / (sqrt(2*pi)*s)
  }
  if (ker=='triangle')
  {
    s <- switch(knorm, sd = sqrt(6), sq = (2/3))
    aux <- (s-abs(x))
    aux <- aux * (aux >0)/s^2
  }
  if (ker=='laplace')
  {
    s <- switch(knorm, sd = 1/sqrt(2), sq = 1/4)
    aux <- exp(-abs(x)/s)/(2*s)
  }
  aux
}

# 2D Kernel and kernel matrix
kstand2d <- function(x, y, ker='normal', knorm="sq") {
  # Densities such as the integral of the square of the density is one if knorm is sq 
  # or such that the sd is one if knorm is sd.
  
  # Compute the 2D kernel value as the product of 1D kernel values
  kx <- kstand(x, ker, knorm)
  ky <- kstand(y, ker, knorm)
  
  k2d <- kx * ky
  k2d
}

## computation of the Omega matrix

wmat <- function(x, h=1 ,ker='normal',knorm="sq",remove=FALSE)
{
  #   ICM and smoothing matrix 
  #   If no bandwidth is provided for the function, h=1 is used  (no smoothing)
  #   The principle diagonal is replaced with zeros if remove = TRUE.
  n <- length(x)
  mat <- sapply(x,function(e) x - e*pracma::ones(n,1))
  wmat <-  kstand(mat / h,ker=ker,knorm=knorm)/h; # kernel smoothing 
  # principle diagonal of the matrix replaced with zeros
  if (remove==TRUE) wmat <-  wmat - diag(diag(wmat))
  wmat
}

# Omega matrix for multidimensional instruments

wmat2 <- function(w1, w2, h = 1, ker = 'normal', knorm = "sq", remove = FALSE) {
  # ICM and smoothing matrix
  # If no bandwidth is provided for the function, h=1 is used (no smoothing)
  # The principle diagonal is replaced with zeros if remove = TRUE.
  n <- length(w1)
  matrix <- cbind(w1, w2)
  for(i in 1:n){
    for (j in 1:n) {
      if(i==1 && j==1) {
          mat <- matrix[i, ]-matrix[j,]
        }
        else{
          mat <- cbind(mat, matrix[i, ]-matrix[j,])
        }
    }
  }
  
  mat <- matrix(mat, nrow=2*n, ncol=n)
  
  wmat <- matrix(0, n, n)
  
  for (i in seq(1, 2*n, by=2)){
    for(j in 1:n){
      wmat[(i+1)/2, j] <- kstand2d(mat[i, j], mat[i+1, j])
    }
  }
  
  wmat
}

# Optimized wmat2 function
wmat3 <- function(w1, w2, h = 1, ker = 'normal', knorm = "sq", remove = FALSE) {
  n <- length(w1)
  mat1 <- outer(w1, w1, "-")
  mat2 <- outer(w2, w2, "-")
  
  # Combine the differences into a single 3-dimensional array
  diffs <- array(c(mat1, mat2), dim = c(n, n, 2))
  
  # Apply the kernel function to each pair of differences
  wmat <- apply(diffs, c(1, 2), function(d) kstand2d(d[1], d[2], ker, knorm))
  
  # Optionally, remove the diagonal elements
  if (remove) {
    diag(wmat) <- 0
  }
  
  wmat
}


tpsmat <- function(knots)
{
  k <- length(knots)
  tmat <- rbind(pracma::ones(1,k),t(knots))
  emat <- sapply(knots, function(x) abs(x-knots)^3)/12
  
  list(tmat=tmat,emat=emat)
}

tpsmat2 <- function(knots1, knots2)
{
  k <- length(knots1)
  
  tmat1 <- rbind(pracma::ones(1,k),t(knots1))
  tmat1 <- rbind(tmat1, knots2)
  tmat1 <- cbind(tmat1, tmat1)
  
  tmat2 <- rbind(pracma::ones(1,k),t(knots1))
  tmat2 <- rbind(tmat2, rep(0, k))
  aux <- rbind(matrix(0, 2, k), knots2)
  tmat2 <- cbind(tmat2, aux)
  
  emat1 <- sapply(knots1, function(x) abs(x-knots1)^3)/12
  emat2 <- sapply(knots2, function(x) abs(x-knots2)^3)/12
  
  emat1 <- rbind(emat1, emat1)
  emat2 <- rbind(emat2, emat2)
  
  emat <- cbind(emat1, emat2)
  
  list(tmat1=tmat1,tmat2=tmat2, emat=emat)
}

tpseval <- function(t,knots=knots,a=a,d=d)
{
  cub <- (abs(t-knots)^3)/12
  c(1,t)%*%a + cub%*%d
}

tpseval21 <- function(t,knots=knots,a=a,d=d)
{
  cub <- (abs(t-knots)^3)/12
  c(1,t)%*%a + cub%*%d
}

tpseval22 <- function(t,knots=knots,a=a,d=d)
{
  cub <- (abs(t-knots)^3)/12
  t*a + cub%*%d
}

tpseval.prime <- function(t,knots=knots,a=a,d=d){ ## This is the function returning the first derivative of ghat.iv
  ## This is written in the same logic as the tps function above
  sqsign <- (3/12)*((t-knots)^2)*ifelse(t-knots>0,1,-1)
  a[2]+sqsign%*%d
}

bmat <- function(l,emat,tmat,W=pracma::eye(n)/n,n)
{
  rbind(cbind(W%*%emat + l*pracma::eye(n),W%*%t(tmat)),cbind(tmat,pracma::zeros(2)))
}

bmat2 <- function(l,emat,tmat1, tmat2, W=pracma::eye(n)/n, n)
{
  # nulmat <- matrix(0, n, n)
  # W1 <- cbind(W, nulmat)
  # W2 <- cbind(nulmat, W)
  # W <- rbind(W1, W2)
  
  rbind(cbind(W%*%emat + l*pracma::eye(2*n),W%*%t(tmat1)),cbind(tmat2,pracma::zeros(3)))
}
## 

fun <- function(x,case=1)  
{
  switch(case,  
         x,
         x^2/sqrt(2),
         sqrt(3*sqrt(3))*x*exp(-(x^2)/2),
         sqrt(3)*x*sin(pi*x/2),
         4*exp(-abs(x)) , 
         log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3),
         sqrt(2)*x*cos(pi*x),
         (log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3) - 0.6*x+ (x^3)*2)/8
  )
}

fun2 <- function(z1, z2,case)  
{
  switch (case, 
          z1,
          z1+z2^2/sqrt(2), #2
          sqrt(3*sqrt(3))*z1*exp(-(z1^2)/2) + z2^2/sqrt(2), #3
          (log(abs(z1-1)+1)*ifelse(z1>1,1,-1)*sqrt(10/3) - 0.6*z1+ (z1^3)*2)/8 + sqrt(3*sqrt(3))*z2*exp(-(z2^2)/2)#4
  )
  
  # sqrt(3*sqrt(3))*x*exp(-(x^2)/2),
  # sqrt(3)*x*sin(pi*x/2),
  # 4*exp(-abs(x)) , 
  # log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3),
  # sqrt(2)*x*cos(pi*x),
  # (log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3) - 0.6*x+ (x^3)*2)/8
  
}

fun21 <- function(z1,case=1)  
{
  switch (case, 
          z1,
          z1, #2
          sqrt(3*sqrt(3))*z1*exp(-(z1^2)/2), #3
          (log(abs(z1-1)+1)*ifelse(z1>1,1,-1)*sqrt(10/3) - 0.6*z1+ (z1^3)*2)/8 #4
  )
  
  # sqrt(3*sqrt(3))*x*exp(-(x^2)/2),
  # sqrt(3)*x*sin(pi*x/2),
  # 4*exp(-abs(x)) , 
  # log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3),
  # sqrt(2)*x*cos(pi*x),
  # (log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3) - 0.6*x+ (x^3)*2)/8
  
}

fun22 <- function(z2,case=1)  
{
  switch (case, 
          z2,
          z2^2/sqrt(2), #2
          z2^2/sqrt(2), #3
          sqrt(3*sqrt(3))*z2*exp(-(z2^2)/2) #4
  )
  
  # sqrt(3*sqrt(3))*x*exp(-(x^2)/2),
  # sqrt(3)*x*sin(pi*x/2),
  # 4*exp(-abs(x)) , 
  # log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3),
  # sqrt(2)*x*cos(pi*x),
  # (log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3) - 0.6*x+ (x^3)*2)/8
  
}

fun.prime <- function(x,case=1){
  switch(case,  
         1,
         x*sqrt(2),
         sqrt(3*sqrt(3))*(1-x^2)*exp(-(x^2)/2),
         sqrt(3)*(sin(pi*x/2)+(x*pi/2)*cos(pi*x/2)),
         4*exp(-abs(x)) , 
         log(abs(x-1)+1)*ifelse(x>1,1,-1)*sqrt(10/3),
         sqrt(2)*( cos(pi*x)-pi*x*sin(pi*x) )
  )
}


wmatp <- function(x, h=1 ,ker='normal',knorm="sq",remove=FALSE)
{
  #   ICM and smoothing matrix 
  #   If no bandwidth is provided for the function, h=1 is used  (no smoothing)
  #   The principle diagonal is replaced with zeros if remove = TRUE.
  n <- dim(x)[1]
  p <- dim(x)[2]
  wmat <- pracma::ones(n,n)
  for (i in 1:p)
  {
    mat <- sapply(x[,i],function(e) x[,i] - e*pracma::ones(n,1))
    wmat <-  wmat*kstand(mat / h,ker=ker,knorm=knorm)/h; # kernel smoothing   
  }
  # principle diagonal of the matrix replaced with zeros
  if (remove==TRUE) wmat <-  wmat - diag(diag(wmat))
  wmat
}

wmat.eval.fun <- function(xdata, xeval , h=1 ,ker='normal', knorm="sd", remove=FALSE)
{
  #   Kernel smoothing matrix having at the row index xeval and column index xdata 
  #   If no bandwidth is provided for the function, h=1 is used  (no smoothing)
  #   The principle diagonal is replaced with zeros if remove = TRUE.
  m <- length(xeval)
  mat <- sapply(xdata,function(e) xeval - e*pracma::ones(m,1))
  wmat <-  kstand(mat / h,ker=ker,knorm=knorm)/h; # kernel smoothing 
  # principle diagonal of the matrix replaced with zeros
  if (remove==TRUE) wmat <-  wmat - diag(diag(wmat))
  wmat
}

RSS <- function(lambda, n, r, A, Astar, AstarA){
  loo.mat <- pracma::ones(n,n) - pracma::eye(n)
  IN <- solve(  as.numeric(lambda)*diag(n) + AstarA  )
  phi.loo <- IN%*%( (Astar*loo.mat)%*%r )
  rss <- (A*loo.mat)%*%phi.loo -  r
  rss <- t(rss)%*%rss
  return(rss)
}

lambda.tikh.fun <- function(y,x,w,hx, hw, lambdamin, lambdamax){
  n <- length(y)
  Kx <- wmat(x=x, h=hx ,ker='normal',knorm="sq",remove=FALSE)
  Kx <- Kx/rowSums(Kx)
  Kw <- wmat(x=w, h=hw ,ker='normal',knorm="sq",remove=FALSE)
  Kw <- Kw/rowSums(Kw)
  r <- Kw%*%y
  KxKw <- Kx%*%Kw
  lambda.star <- pracma::fminbnd(function(e) RSS(e, n=n, r=r, A=Kw, Astar=Kx, AstarA=KxKw ), lambdamin, lambdamax)$xmin    
  ## phihat <- solve(  as.numeric(lambda.star)*diag(n) + KxKw  )%*%Kx%*%r
  Hatmat <- solve(  as.numeric(lambda.star)*diag(n) + KxKw  )%*%Kx%*%Kw
  phihat <- Hatmat%*%y
  df <- sum(diag(Hatmat))
  return(list(phihat=phihat , lambda.star=lambda.star, df=df))
}

phihat.eval.tikh.fun <- function(y,x,w,zeval,hx, hw, lambda, phihat.data){
  Kw <- wmat(x=w, h=hw ,ker='normal',knorm="sq",remove=FALSE)
  Kw <- Kw/rowSums(Kw) 
  Kxeval <- wmat.eval.fun(xdata=x, xeval=zeval , h=hx ,ker='normal', knorm="sd", remove=FALSE)
  Kxeval <- Kxeval/rowSums(Kxeval)
  phihat.eval <- (1/as.numeric(lambda))*(Kxeval%*%Kw%*%(y-phihat.data) )
  return(phihat.eval)
}


hermite.fun <- function(x,d){
  n <- length(x)
  her <- matrix(rep(NA,n*(d+1)), nrow=n,ncol=d+1)
  for(j in 0:d){
    herj <- EQL::hermite(x,j,prob=TRUE)  
    her[,j+1] <- herj
  }
  return(her)
}


ghat.gal.fun <- function(y,x,w, zeval){
  n <- length(y)
  ## tun <- n^(-1/2) #" This is used for tuning , necessary to improve the performance of the galerking estimator
  tun <- 0 #" This is used for tuning , necessary to improve the performance of the galerking estimator
  ## wnorm <- (w-min(w))/(max(w)-min(w)) ## maybe these should be in [-1,1] instead of [0,1]
  ## wnorm <- (2*(w-min(w))/(max(w)-min(w)))-1 ## normalization in [-1,1]
  ## xnorm <- (x-min(x))/(max(x)-min(x)) ## maybe these should be in [-1,1] instead of [0,1]
  ## xnorm <- ( 2*(x-min(x))/(max(x)-min(x)) ) - 1 ## Normalization in [-1,1] 
  nk <- 0
  crit <- -100
  while(crit<0){
    nk <- nk+1
    ## xx <- t(pracma::legendre(nk,xnorm)) ## n by nk matrix
    ## ww <- t(pracma::legendre(nk+1,wnorm)) ## n by (nk+1) matrix
    xx <- hermite.fun(x=x,d=nk)
    ww <- hermite.fun(x=w,d=nk+1)
    ## pmat <- t(xx)%*%( ww%*%solve(qr(t(ww)%*%ww, LAPACK=TRUE))%*%t(ww) )%*%xx
    ## pmat <- t(xx)%*%( ww%*%pracma::pinv(t(ww)%*%ww)%*%t(ww) )%*%xx
    pmat <- t(xx)%*%( ww%*%pracma::pinv(t(ww)%*%ww + tun*diag(nk+2))%*%t(ww) )%*%xx
    #pmat <- ( xx%*%pracma::pinv(t(xx)%*%xx + tun*diag(nk+1))%*%t(xx) )%*%( ww%*%pracma::pinv(t(ww)%*%ww + tun*diag(nk+2))%*%t(ww) ) ## tuned regularizaed version
    ## pmatx <- xx%*%solve(qr(t(xx)%*%xx, LAPACK=TRUE))%*%t(xx)
    ## pmatw <- ww%*%solve(qr(t(ww)%*%ww, LAPACK=TRUE))%*%t(ww)
    ## pmat <- pmatx%*%pmatw
    eigenvalues <- eigen(pmat)$values 
    rho2 <- 1/Re(eigenvalues[which.min(abs(eigenvalues))])
    crit <- rho2*(nk^3.5)/n - 1
  }
  ## qmat <- t(xx)%*%(  ww%*%solve( qr(t(ww)%*%ww, LAPACK=TRUE) )%*%t(ww)  )%*%xx
  ## qmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww)%*%t(ww)  )%*%xx
  qmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww+tun*diag(nk+2))%*%t(ww)  )%*%xx ## tuned regularizaed version
  ## nmat <- t(xx)%*%(  ww%*%solve( qr(t(ww)%*%ww, LAPACK=TRUE) )%*%t(ww)  )
  ## nmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww)%*%t(ww)  )
  nmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww+tun*diag(nk+2))%*%t(ww)  ) ## tuned regularizaed version
  nmaty <- nmat%*%y
  ## bhat_G <- solve( qr( qmat  , LAPACK=TRUE), nmaty )
  bhat_G <- solve( qr( qmat +tun*diag(nk+1)  , LAPACK=TRUE), nmaty ) ## tuned regularizaed version
  
  phi_hat_G <- xx%*%bhat_G
  Jhatfun <- rep(NA,nk)
  for(j in 1:nk){
    ## xx <- t(pracma::legendre(j, xnorm)) ## n by j matrix
    xx <- hermite.fun(x=x,d=j)
    ## ww <- t(pracma::legendre(j+1,wnorm)) ## n by (j+1) matrix
    ww <- hermite.fun(x=w,d=j+1)
    ## xxinv <- pracma::pinv(t(xx)%*%xx)
    xxinv <- pracma::pinv(t(xx)%*%xx+tun*diag(j+1)) ## tuned regularizaed version
    ## Amin1 <- ww%*%( solve( qr(t(ww)%*%(xx%*%xxinv%*%t(xx))%*%ww, LAPACK=TRUE))%*%t(ww) )
    Amin1 <- ww%*%( solve( qr(t(ww)%*%(xx%*%xxinv%*%t(xx))%*%ww+tun*diag(j+2), LAPACK=TRUE))%*%t(ww) ) ## tuned regularizaed version
    ## qmat <- t(xx)%*%(  ww%*%solve( qr(t(ww)%*%ww, LAPACK=TRUE) )%*%t(ww)  )%*%xx 
    ## qmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww)%*%t(ww)  )%*%xx
    qmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww+tun*diag(j+2))%*%t(ww)  )%*%xx ## tuned regularizaed version
    ## nmat <- t(xx)%*%(  ww%*%solve( qr(t(ww)%*%ww, LAPACK=TRUE) )%*%t(ww)  )
    ## nmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww)%*%t(ww)  )
    nmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww+tun*diag(j+2))%*%t(ww)  ) ## tuned regularizaed version
    nmaty <- nmat%*%y
    ## Jhatfun[j] <- (2/3)*(log(n)/(n^2))*sum( ((y-phi_hat_G)^2)*rowSums((Amin1%*%xx)^2)  ) - sum( (xx%*%solve(qr(qmat, LAPACK=TRUE), nmaty))^2 )
    Jhatfun[j] <- (2/3)*(log(n)/(n^2))*sum( ((y-phi_hat_G)^2)*rowSums((Amin1%*%xx)^2)  ) - sum( (xx%*%solve(qr(qmat+tun*diag(j+1), LAPACK=TRUE), nmaty))^2 ) ## tuned regularizaed version
  }
  nk <- which.min(Jhatfun)
  ## xx <- t(pracma::legendre(nk,xnorm)) ## n by nk matrix
  xx <- hermite.fun(x=x,d=nk)
  ## ww <- t(pracma::legendre(nk+1,wnorm))  ## n by (nk+1) matrix
  ww <- hermite.fun(x=w,d=nk+1)
  ## qmat <- t(xx)%*%(  ww%*%solve( qr(t(ww)%*%ww, LAPACK=TRUE) )%*%t(ww)  )%*%xx
  ## qmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww)%*%t(ww)  )%*%xx
  qmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww+tun*diag(nk+2))%*%t(ww)  )%*%xx ## tuned regularizaed version
  ## nmat <- t(xx)%*%(  ww%*%solve( qr(t(ww)%*%ww, LAPACK=TRUE) )%*%t(ww)  )
  ## nmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww)%*%t(ww)  )
  nmat <- t(xx)%*%(  ww%*%pracma::pinv(t(ww)%*%ww+tun*diag(nk+2))%*%t(ww)  ) ## tuned regularizaed version
  nmaty <- nmat%*%y
  ## bhat_G <- solve( qr( qmat  , LAPACK=TRUE), nmaty )
  bhat_G <- solve( qr( qmat +tun*diag(nk+1) , LAPACK=TRUE), nmaty ) ## tuned regularizaed version
  zevalnorm <- (zeval-min(zeval))/(max(zeval)-min(zeval)) ## Normalization in [0,1]
  ## zevalnorm <- ( 2*(zeval-min(zeval))/(max(zeval)-min(zeval)) ) -1 ## Normalization in [-1,1]
  ## phi_hat_G_zeval <- t(pracma::legendre(nk,zevalnorm))%*%bhat_G ## phihat at the evaluation points
  phi_hat_G_zeval <- hermite.fun(x=zeval,d=nk)%*%bhat_G ## phihat at the evaluation points
  ## phi_hat_G <- t(pracma::legendre(nk,xnorm))%*%bhat_G ## phihat at the data points
  phi_hat_G <- hermite.fun(x=x,d=nk)%*%bhat_G ## phihat at the data points
  Hatmat <- hermite.fun(x=x,d=nk)%*%solve( qr( qmat +tun*diag(nk+1) , LAPACK=TRUE))%*%nmat
  df <- sum(diag(Hatmat))
  return(list(ghat_gal_eval=phi_hat_G_zeval, ghat_gal =phi_hat_G , nk=nk, df=df))
}

fun.sam <- function(x,case=1)  
{
  switch(case,  
         x,
         sqrt(2)*(x^2),
         0.5*sin(1.5*pi*x))
}

## Functions to implement Horowitz' procedure. 


## AAA: Need to install the package "orthopolynom" needed for legendre.basis.fun below (that is in turn needed for Horowitz estimator)


power.mat.fun <- function(x,degree){ ## here x can be a vector of observations on x
  x <- as.vector(x)
  n <- length(x)
  xmat.pow <- matrix(rep(NA,(degree+1)*n), nrow=degree+1, ncol=n)
  for(s in 1:n){
    xvec.pow <- rep(NA,degree+1)
    for(j in 0:(degree)){
      xj <- as.numeric(x[s])^j
      xvec.pow[j+1] <- xj
    }
    xmat.pow[,s] <- xvec.pow
  }
  return(xmat.pow) 
}

legendre.basis.fun <- function(x,degree){
  
  leg.list <- orthopolynom::legendre.polynomials(degree, normalized=TRUE)
  leg <- unlist(leg.list)
  
  leg.mat <- matrix(rep(0,(degree+1)^2), nrow=(degree+1), ncol=(degree+1))
  m0 <- 1
  for(m in 0:(degree)){ ## Build the matrix of coefficients
    a <- rep(0,(degree+1))
    am <- leg[max(m0,1):(m0+m)]
    a[1:(length(am))] <- am ## *sqrt(m+.5) 
    leg.mat[m+1,] <- a
    m0 <- m0+m+1
  }
  
  power.mat <- power.mat.fun(x=x,degree=degree) ## (degree+1) times n matrix of basis function:: jeneric element \psi_j(X_s) (j=row, s=col)
  legendre.basis <- t(leg.mat%*%power.mat) ## Gives the matrix \cal{X} or \cal{W} of the notes
  return(legendre.basis)
}
## The one below run provides the Horowitz' originary estimator
ghat.horow.fun <- function(y,z,w,zeval){
  PW <- ecdf(w)
  PZ <- ecdf(z)
  z.norm <- PZ(z) ## normalize on [0,1]
  w.norm <- PW(w)
  zeval.norm <- PZ(zeval)
  
  n <- length(y)
  tun <- 0.001 ## try diferent values here, to see the performance
  nk <- 0
  crit <- -100
  
  while(crit<0){ ## This selects the J_n in Horowitz
    nk <- nk+1
    ## xx <- hermite.fun(x=x,d=nk) ## n X nk matrix
    ## ww <- hermite.fun(x=w,d=nk) ## n X nk matrix
    z.input <- 2*z.norm - 1 ## To shift the Legendre polynomials from [-1,1] to [0,1]
    w.input <- 2*w.norm - 1
    xx <- legendre.basis.fun(x=z.input,degree=nk)*sqrt(2) ## n X (nk+1) matrix; \sqrt{2} is used to guarantee the normalization to 1 of the shifted Legendre Polynomials
    ww <- legendre.basis.fun(x=w.input,degree=nk)*sqrt(2) ## n X (nk+1) matrix
    cmat <- t(xx)%*%ww%*%t(ww)%*%xx/(n^2)
    cmat2 <- cmat%*%t(cmat)
    eigenvalues <- eigen(cmat2)$values 
    rho2 <- 1/Re(eigenvalues[which.min(abs(eigenvalues))])
    crit <- rho2*( (nk+1)^3.5)/n - 1 ## The number of basis functions is nk+1 and not nk, as there is the constant basis function!
  }
  
  mhat.coef <- t(ww)%*%y/n
  qmat <- t(ww)%*%xx/n
  g.tilde.coef <-  solve( qr( qmat +tun*diag(nk+1)  , LAPACK=TRUE), mhat.coef )
  g.tilde.coef <- as.vector(g.tilde.coef)
  g.tilde <- xx%*%g.tilde.coef
  
  Ahat.minus1.star.psi.mat <-   solve( qr( qmat +tun*diag(nk+1)  , LAPACK=TRUE), t(ww) )
  Ahat.minus1.star.psi.mat2 <- Ahat.minus1.star.psi.mat^2
  Jhatfun <- rep(NA,nk)
  for(j in 2:(nk+1)){ ## Start from j=2, as j=1 is just the constant 
    ##   j <- 2
    ## xxj <- xx[,1:j]
    ## wwj <- ww[,1:j]
    g.tilde.coef.j <- g.tilde.coef[1:j]
    ## ghatj <- xxj%*%g.tilde.coef.j
    ## gamma.mat.Astarminus1.j <- gamma.mat.Astarminus1[,1:j]
    Ahat.minus1.star.psi.mat2_j <- Ahat.minus1.star.psi.mat2[1:j,]
    ## dim(Ahat.minus1.star.psi.mat2_j)
    Sumj.Ahat.minus1.star.psi.mat2_j <- colSums(Ahat.minus1.star.psi.mat2_j)
    ## length(Sumj.Ahat.minus1.star.psi.mat2_j)
    Jhatfun_j <- (2/3)*(log(n)/(n^2))*sum(  ((y-g.tilde)^2)*Sumj.Ahat.minus1.star.psi.mat2_j) - sum((g.tilde.coef.j)^2)  
    Jhatfun[j-1] <- as.numeric(Jhatfun_j) ## We put j-1 as an index, as we start from j=2
  }
  nk.star <- which.min(Jhatfun)+1 # truncation parameter selected with Horowitz method
  xx_nk <- xx[, 1:(nk.star)] 
  ## ww_nk <- ww[, 1:nk]
  ## xx <- hermite.fun(x=x,d=nk)
  ## ww <- hermite.fun(x=w,d=nk)
  ghat.coef <- g.tilde.coef[1:nk.star]
  ghat.horow.data <- as.vector(xx_nk%*%ghat.coef)
  df.mat <- xx_nk%*%(solve( qr( qmat +tun*diag(nk+1)  , LAPACK=TRUE), t(ww)/n )[1:nk.star, ])
  df <- sum(diag(df.mat))
  
  
  zeval.input <- 2*zeval.norm - 1
  zz.eval <- legendre.basis.fun(x=zeval.input,degree=(nk.star-1))*sqrt(2) ## length(zeval) X nk.star matrix (Recall that legendre basis contain constant which is the 0 degree )
  ghat.horow.eval <- as.vector(zz.eval%*%ghat.coef)
  
  return(list(ghat.horow.eval=ghat.horow.eval , ghat.horow.data=ghat.horow.data, 
              nk=nk.star, df=df))
}


