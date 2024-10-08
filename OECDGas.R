# set your working directory

rm(list=ls())

source("bits_comparison.R")
library(pracma)
library(compiler)
library(orthopolynom)
#library(MASS)
library(foreach)
#library(matrixcalc)
library(parallel)
library(doParallel)
library(doRNG)
library(scatterplot3d)
library(plotly)
library(mgcv)
library(dplyr)
library(kableExtra)
# load the data set and get an overview
library(AER)
library(dplyr)

data("OECDGas")
# summary(OECDGas)

# first-differencing the data to remove country fixed-effects (see README for the model specification)

OECDGas_diff <- OECDGas %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(
    d_gas = diff(c(NA, gas)),
    d_price = diff(c(NA, price)), ## for IV adn OLS
    d_income = diff(c(NA, income)), ## for IV and OLS
    d_cars = diff(c(NA, cars)),
    lag_price = lag(price),  # Compute difference with lagged price
    #lag_price_2 = c(NA, lag(d_price)), # Second lag
    price_level = exp(price), 
    lag_price_level = exp(lag_price)
 )%>%
  ungroup()

OECDGas_diff <- OECDGas_diff %>%
  arrange(country, year) %>%
  group_by(year) %>%
  mutate(
    avg_price_other_countries = ifelse(sum(country != first(country) & year == first(year)) > 0,
                                       mean(price[country != first(country) & year == first(year)], na.rm = TRUE),
                                       NA_real_),
    lag_avg_price_other_countries = ifelse(sum(country != first(country) & year == first(year)) > 0,
                                       mean(lag_price[country != first(country) & year == first(year)], na.rm = TRUE),
                                       NA_real_)
  )%>%
  ungroup()

OECDGas_diff <- OECDGas_diff %>%
  arrange(country, year) %>%
  group_by(year) %>%
  mutate(
    avg_level = ifelse(sum(country != first(country) & year == first(year)) > 0,
                                       mean(price_level[country != first(country) & year == first(year)], na.rm = TRUE),
                                       NA_real_),
    lag_avg_level = ifelse(sum(country != first(country) & year == first(year)) > 0,
                                           mean(lag_price_level[country != first(country) & year == first(year)], na.rm = TRUE),
                                           NA_real_), 
    avg_ln = log(avg_level), 
    lag_avg_ln = log(lag_avg_level)
  )%>%
  ungroup()

# Remove rows with NA in any column
OECDGas_diff <- na.omit(OECDGas_diff)

# specifying the distribution

kernel <- "laplace"

# number of folds used in CV
nfold <- 4

# variables used to assess lambda on
sdu <- 1
ll <- 30

pmin <- .00001
pmax <- .5

pp <- seq(pmin,pmax,length=ll)
lambda <- (pp/(1-pp))/(sdu^2)

# number of repetitions
nres <- 50

for (sim_n in 1:nres){
    OECDGas_diff <- OECDGas_diff[sample(nrow(OECDGas_diff)),]
    n <- nrow(OECDGas_diff)

    z1 <- OECDGas_diff$price
    z2 <- OECDGas_diff$lag_price
    
    # instruments
    w1 <- OECDGas_diff$avg_price_other_countries
    w2 <- OECDGas_diff$lag_avg_price_other_countries

    
    y <- OECDGas_diff$d_gas
    
    sfold <- as.integer(n/nfold)
    
    W <- wmat3(w1, w2, h=1, ker=kernel, knorm="sq", remove=FALSE)/n #this is Omega if only Z1 is endogenous
    
    obj <- rep(0,ll)
    
    nulmat <- matrix(0, (nfold-1)*sfold, (nfold-1)*sfold)
    
    # finding the optimal lambda
    # n-fold CV
    for (m in (1:ll)){
      
      l <- lambda[m]
      giv<- rep(0,n)
      g1iv <- rep(0, n)
      g2iv <- rep(0, n)
      
      for (k in (0:(nfold-1))){
        
        Wf <- W[-((k*sfold+1):((k+1)*sfold)), -((k*sfold+1):((k+1)*sfold))]
       
        
        matsf <- tpsmat2(z1[-((k*sfold+1):((k+1)*sfold))], z2[-((k*sfold+1):((k+1)*sfold))])
        tmatf1 <- matsf$tmat1
        tmatf2 <- matsf$tmat2
        ematf <- matsf$emat
        
        # e1 <- ematf[1:((nfold-1)*sfold),1:((nfold-1)*sfold)]
        # e2 <- ematf[1:((nfold-1)*sfold),(sfold+1):n]
        
        W1 <- cbind(Wf, nulmat)
        W2 <- cbind(nulmat, Wf)
        Wf <- rbind(W1, W2)
        
        # the matrix we need to invert
        bigmativf <- bmat2(l, ematf, tmatf1, tmatf2, Wf, n-sfold)
        bigwyf <- c(Wf%*%c(y[-((k*sfold+1):((k+1)*sfold))],y[-((k*sfold+1):((k+1)*sfold))])
                    ,0,0,0)
        parestivf <- solve(qr(bigmativf,LAPACK=TRUE),bigwyf)
        
        
        # recovering the function g
        for (j in ((k*sfold+1):((k+1)*sfold))){
          g1iv[j] <-  parestivf[(2*(n-sfold)+1)] + z1[j]*parestivf[(2*(n-sfold)+2)] + (1/12)*sum(parestivf[1:(n-sfold)]*(abs(z1[j]-z1[-((k*sfold+1):((k+1)*sfold))]))^3)
          g2iv[j] <-  z2[j]*parestivf[(2*(n-sfold)+3)] + (1/12)*sum(parestivf[1:(n-sfold)]*(abs(z2[j]-z2[-((k*sfold+1):((k+1)*sfold))]))^3)
          giv[j] <- g1iv[j] + g2iv[j]
        }
        
      }
      
      civ <- t(y-giv)%*%W%*%(y-giv)
      obj[m] <- civ

    }
    
    idl <- which(obj==min(obj)) ## the index for which lambda is minimized
    lcviv <- lambda[idl] # take that lambda
    
    nulmat <- matrix(0, n, n)
    
    # full sample estimates
    
    mats <- tpsmat2(z1, z2)
    tmat1 <- mats$tmat1
    tmat2 <- mats$tmat2
    emat <- mats$emat
    
    e1 <- emat[(n+1):(2*n),1:n]
    e2 <- emat[1:n,(n+1):(2*n)]
    
    W1 <- cbind(W, nulmat)
    W2 <- cbind(nulmat, W)
    W <- rbind(W1, W2)
    
    # the matrix we need to invert
    bigmativ <- bmat2(lcviv, emat, tmat1, tmat2, W, n)
    bigwy <- c(W%*%c(y,y),0,0,0)
    
    parestiv <- solve(qr(bigmativ,LAPACK=TRUE),bigwy)
    
    Z1 <- rbind(rep(1, n), z1)
    Z2 <- z2
    
    e1d1 <- cbind(e1, nulmat)%*%parestiv[1:(2*n)]
    e2d2 <- cbind(nulmat, e2)%*%parestiv[1:(2*n)]
    z1a1 <- t(Z1)%*%matrix(parestiv[((2*n)+1):((2*n)+2)],nrow=2,ncol=1)
    z2a2 <- matrix(Z2, nrow=n, ncol=1)%*%matrix(parestiv[(2*n)+3], nrow=1, ncol=1)
    
    # recovering the functions g1 and g2
    
    g1iv <- e1d1 + z1a1
    g2iv <- e2d2 + z2a2
    # 
    giv <- rep(0, n)
    giv <- g1iv + g2iv
    
    leval<-100
    
    zeval <- seq(-1,0, length=leval)
    
    ghativ1_aux <- unlist(sapply(zeval,tpseval21,knots=z1,a=parestiv[(2*n+1):(2*n+2)],d=parestiv[1:n]))
    ghativ2_aux <- unlist(sapply(zeval,tpseval22,knots=z2,a=parestiv[2*n+3],d=parestiv[(n+1):(2*n)]))
    
    # normalizing 
    
    ghativ1_aux <- ghativ1_aux - (ghativ1_aux[leval/2]+ghativ1_aux[leval/2 + 1])/2
    ghativ2_aux <- ghativ2_aux - (ghativ2_aux[leval/2]+ghativ2_aux[leval/2 + 1])/2
    
    ## Estimating GAM
    ## Imposing g1(0)=g2(0)=0
    
    gam <- mgcv::gam(y ~ s(z1) + s(z2))
    prediction1 <- predict(gam, data.frame(z1=zeval, z2=rep(0, leval)))
    prediction2 <- predict(gam, data.frame(z1=rep(0, leval), z2=zeval))
    

    if(sim_n==1){
      res <- rbind(cbind(ghativ1_aux, ghativ2_aux), c(lcviv, lcviv))
      res <- cbind(res, rbind(cbind(prediction1, prediction2), c(0, 0)))
    }else{
      res <- cbind(res, rbind(cbind(ghativ1_aux, ghativ2_aux), c(lcviv, lcviv)))
      res <- cbind(res, rbind(cbind(prediction1, prediction2), c(0, 0)))
    }
}

  filename1 <- paste ("OECDGas", nres, ".R")
  filenameplot1 <- paste ("OECDGas_plot_1_", nres, ".pdf")
  
  #filename2 <- paste ("sLuka2", "_rho_zw1"  , rho_zw1,"_rho_zw2"  , rho_zw2, "_rhouv", rho_uv, "_case",case,"n",n, ".R" ,sep = "")
  filenameplot2 <- paste ("OECDGas_plot_2_", nres)
  
  # the multiples of 2 are the knots for g2, and odd columns are the knots for g1
  save(res, file=filename1)
  
  
  #########################
  ##  PLOTTING & RESULTS  #
  #########################
  
  load(filename1)
  
  ## dividing the data in 4 parts for g1 and g2, gam1 and gam2
  
  s1 <- seq(1, 4*nres, by=4)
  s2 <- seq(2, 4*nres, by=4)
  s3 <- seq(3, 4*nres, by=4)
  s4 <- seq(4, 4*nres, by=4)
  leval <-100
  g1hativ <- res[1:leval, s1]
  g2hativ <- res[1:leval, s2]
  
  meang1 <- apply(g1hativ, 1, mean)
  meang2 <- apply(g2hativ, 1, mean)
  
  ## NORMALIZE by imposing g_i(0)=0
  g1_norm <- g1hativ - (meang1[leval/2]+meang1[leval/2+1])/2
  g2_norm <- g2hativ - (meang2[leval/2]+meang2[leval/2+1])/2
  
  # update meang1 & meang2
  meang1_norm <- apply(g1hativ, 1, mean) - (meang1[leval/2]+meang1[leval/2+1])/2
  meang2_norm <- apply(g2hativ, 1, mean) - (meang2[leval/2]+meang2[leval/2+1])/2
  
  ## Using GAM to compare our estimator (package "mgcv")
  
  gam1 <- res[1:leval,s3]
  gam2 <- res[1:leval,s4]
  meangam1 <- apply(gam1, 1, mean)
  meangam2 <- apply(gam2, 1, mean)

  ## saving
  #pdf(filenameplot1)
  #plot(zeval,f1, type = "l",col="black", lty="solid", lwd=2,xlab = "Z1", ylab="Y", ylim=c(lbound1,ubound1), cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
  plot(zeval,meang1_norm, type = "l", lty="solid",col="#708090", xlab = "ln(price)", ylab="∆ln(gas)", lwd=2)
  legend("topright", legend = c("SS(t)", "SS(t-1)","SSWE(t)" ,"SSWE(t-1)" ), lty=c("solid", "dashed", "solid","dashed"), col = c("red", "red","#708090", "#708090"), cex = 0.8)
  lines(zeval,meangam1, type = "l", lty="solid",col="red", lwd=2)
  lines(zeval,-meang2_norm, type = "l", lty="dotdash",col="#708090", lwd=1)
  lines(zeval,-meangam2, type = "l", lty="dashed",col="red", lwd=1)
  # legend("topleft", legend=c("True", "Smoothing splines", "Smoothing splines without endogeneity"),
  #       col=c("black", "#708090", "red"), lty=c("solid", "dotdash","dashed"), lwd=2, cex=1.5)
