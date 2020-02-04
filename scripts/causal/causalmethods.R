library(np)
library(MASS)
library(glmpath)
library(glmnet)

# adaptive LASSO
VTM <- function(vc, dm) {
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

alasso = function(data,wgt=NULL,rtn="EST",nopen.ind=NULL,regularize=TRUE,
                  BIC.factor=0.1,offset=NULL,fam="binomial",relax=FALSE) {
  fam0=fam;Wi=wgt;
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); pp = ncol(x)
  if(is.null(Wi)){Wi=rep(1,nn)}; if(is.null(offset)){offset=rep(0,nn)}
  if(regularize){
    ## adaptive lasso (aLASSO) with initial estimator obtained via ridge (bini); w.b creates adaptive weights for aLASSO ##
    lam.ridge = log(pp)/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef; 
    # bini = as.vector(coef(glmnet(x,y,weights=Wi,alpha=0,standardize=F,lambda=lam.ridge,family=fam0,offset=offset))) # To check: can remove offset
    bini = as.vector(coef(glmnet(x,y,weights=Wi,alpha=0,standardize=F,lambda=lam.ridge,family=fam0)))
    w.b = 1/abs(bini[-1]); x.t = x/VTM(w.b,nrow(x))
    
    ## glmpath provides solution path for a range of penalty parameters ##
    tmpfit = glmpath(x.t,y,nopenalty.subset=nopen.ind,family=fam0,weight=Wi,standardize=F,min.lambda=0,
                     lambda2=lam.ridge,offset=offset)
    #tmpfit = glmnet(x.t,y,family=fam0,weight=Wi,standardize=F,alpha=1,nlambda=500)
    lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
    b.all = predict(tmpfit, s=lam.all, type="coefficients",mode="lambda",offset=offset)
    b.all = b.all/VTM(c(1,w.b),nrow(b.all)); m0 = length(lam.all)
    ## ================================================================================ ##
    ## calculates degree of freedom for all beta's (corresponding to different lam.all) ##
    ## ================================================================================ ##
    df.all = apply(b.all[,-1,drop=F]!=0,1,sum); x = as.matrix(x)
    ## =============================================================================================================== ##
    ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
    ## =============================================================================================================== ##
    if(length(unique(y))==2){n.eff = sum(y)}else{n.eff = nn}
    BIC.lam = -2*apply(predict(tmpfit,newx=x.t,newy=y,s=lam.all,type="loglik",mode="lambda",offset=offset),2,sum)+min(n.eff^BIC.factor,log(n.eff))*df.all 
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
    if(relax){
      ind.out = bhat[-1] == 0; 
      tmpb = glm(y~cbind(1,x[,!ind.out])-1,family=fam0,weight=Wi)$coef
      bhat[c(T,!ind.out)] = tmpb; bhat[c(F,ind.out)]=0
    }
  } else{
    ## ========================================================================= ##
    ## if regularize = F, then use standard logistic regression w/o penalization ##
    ## ========================================================================= ##
    bhat=bini=glm(y~x,family=fam0,weight=Wi)$coef;lamhat = 0; lam.all=BIC.lam=b.all=NULL
  }
  out = c("b"=bhat, "bini"=bini,"lamhat"=lamhat,"lam.all"=lam.all,"BIC.lam"=BIC.lam,"b.all"=b.all)
  if(rtn=="EST"){return(bhat)}else{return(list(out,"b.all"=b.all,"lam.all"=lam.all,"fit"=tmpfit,"BIC.lam"=BIC.lam))}
}

# ridge
ridge <- function(dat, fam, wgt=NULL) {
  if(is.null(wgt)) wgt=rep(1,nrow(dat))
  
  y <- dat[,1]
  x <- data.matrix(dat[,-1])
  
  n <- length(y)
  p <- ncol(x)
  
  fit <- cv.glmnet(x, y, weights=wgt, family=fam, nfolds=5, alpha=0, nlambda=100) # length=100, nfold=10
  rv <- as.numeric(coef(fit, s=fit$lambda.min))
}

# methods
expit <- function(v) {
  1/(1+exp(-v))
}

ipw <- function(Yi,Ti,pihat,normalize=T) {
  n <- length(Yi)
  
  Wi1 <- (Ti/pihat)/n
  Wi0 <- ((1-Ti)/(1-pihat))/n
  
  if(normalize) {
    Wi1 <- Wi1/sum(Wi1)
    Wi0 <- Wi0/sum(Wi0)
  }
  
  delta <- sum((Wi1-Wi0)*Yi)
}

aipw <- function(Yi,Ti,pihat,muhat1,muhat0) {
  n <- length(Yi)
  
  Wi1 <- (Ti/pihat)/n
  Wi0 <- ((1-Ti)/(1-pihat))/n
    
  mu1 <- sum(Wi1*(Yi-muhat1))+mean(muhat1)
  mu0 <- sum(Wi0*(Yi-muhat0))+mean(muhat0)
  
  delta <- mu1-mu0
}

dips <- function(Yi,Ti,Xi,Gi,fam="gaussian") {
  Zi.ps <- model.matrix(Ti~Xi, data=data.frame(Xi))
  Zi.rp <- model.matrix(Yi~Ti+Xi, data=data.frame(Xi,Ti))
  
  b.ps.alas <- alasso(cbind(Ti, Zi.ps[,-1]), fam="binomial")
  b.rp.alas <- alasso(cbind(Yi, Zi.rp[,-1]), fam=fam,wgt=Gi)
  b.ps.alas[is.na(b.ps.alas)] <- 0
  b.rp.alas[is.na(b.rp.alas)] <- 0
  
  ps.alas <- expit(Zi.ps %*% b.ps.alas)
  xb.alas <- cbind(Zi.rp[,-c(2)]) %*% b.rp.alas[-c(2)]
  scr.xb <- as.numeric(pnorm(xb.alas, mean=mean(xb.alas), sd=sd(xb.alas)))
  scr.ps <- as.numeric(pnorm(ps.alas, mean=mean(ps.alas), sd=sd(ps.alas)))
  
  if(all(b.ps.alas[-1]==0) & !all(b.rp.alas[-c(1,2)]==0)) {
    bwd <- sd(scr.xb)/n^(1/6)
    pi.dips <- npreg(bws=bwd, txdat=cbind(scr.xb), tydat=Ti, exdat=cbind(scr.xb), ckerorder=4)$mean
  } else if(!all(b.ps.alas[-1]==0) & all(b.rp.alas[-c(1,2)]==0)) {
    bwd <- sd(scr.ps)/n^(1/6)
    pi.dips <- npreg(bws=bwd, txdat=cbind(scr.ps), tydat=Ti, exdat=cbind(scr.ps), ckerorder=4)$mean
  } else if(all(b.ps.alas[-1]==0) & all(b.rp.alas[-c(1,2)]==0)) {
    pi.dips <- mean(Ti)
  } else{
    bwd <- apply(cbind(scr.ps,scr.xb),2,sd)/n^(1/6)
    pi.dips <- npreg(bws=bwd, txdat=cbind(scr.ps, scr.xb), tydat=Ti, exdat=cbind(scr.ps, scr.xb), ckerorder=4)$mean
  }
  
  pi.dips
}

aipw.ss <- function(Yi,Ti,Ri,pihat,nuhat,muhat1,muhat0,xihat1,xihat0) {
  Wi.Ti1 <- (Ti/pihat)
  Wi.Ti0 <- ((1-Ti)/(1-pihat))
  Wi.Ri <- (Ri/nuhat)
  
  mu1 <- mean(Wi.Ri*Wi.Ti1*Yi - (Wi.Ri-1)*Wi.Ti1*xihat1 - (Wi.Ti1-1)*muhat1)
  mu0 <- mean(Wi.Ri*Wi.Ti0*Yi - (Wi.Ri-1)*Wi.Ti0*xihat0 - (Wi.Ti0-1)*muhat0)
  
  delta <- mu1-mu0
}

# function for computing truncated cubic 
tcub=function(x,xk){#xk is the knot location, x is the variable
  ((x>xk)*(x-xk))^3 
}
# function for computing natural spline basis
nsbasis=function(x,nknot){
  
  knots=quantile(x,seq(0,1,length=nknot))
  
  H=x
  if (length(unique(x))>4){# inactive if x is a discrete covariate with less than 4 levels;
    for (i in 1:(nknot-2)){
      d_i=tcub(x,knots[i])/(knots[nknot]-knots[i])
      d_k=tcub(x,knots[nknot-1])/(knots[nknot]-knots[nknot-1])
      
      Hnew=d_i-d_k
      H=cbind(H,Hnew)
    }
  }
  H
}

linear <- function(v) {
  rv <- v
}

g <- linear # link function

dips.ss <- function(Yi,Ti,Ri,Wi,Xi,Gi,offLog=F,fam="gaussian") {
  # n <- length(Yi)
  Zi.ps <- model.matrix(Ti~Xi, data=data.frame(Xi))
  Zi.rp <- stats::model.matrix(Yi~Ti+Xi, data=data.frame(Xi,Ti))
  p <- dim(Zi.ps)[2]-1
  
  if(ncol(Xi)==1) {
    b.ps.alas <- glm.fit(Zi.ps,Ti,family=binomial())$coef
    b.ps.alas[is.na(b.ps.alas)] <- 0
    ps.alas <- expit(Zi.ps %*% b.ps.alas)
  } else {
    b.ps.alas <- alasso(cbind(Ti, Zi.ps[,-1]), fam="binomial")
    b.ps.alas[is.na(b.ps.alas)] <- 0
    ps.alas <- expit(Zi.ps %*% b.ps.alas)
  }
  b.rp.alas <- alasso(cbind(Yi, Zi.rp[,-1])[Ri==1,], fam=fam, wgt=Gi[Ri==1])
  b.rp.alas[is.na(b.rp.alas)] <- 0
  
  xb.alas <- cbind(Zi.rp[,-c(2)]) %*% b.rp.alas[-c(2)]
  scr.xb <- as.numeric(pnorm(xb.alas, mean=mean(xb.alas), sd=sd(xb.alas)))
  scr.ps <- as.numeric(pnorm(ps.alas, mean=mean(ps.alas), sd=sd(ps.alas)))
  
  if(all(b.ps.alas[-1]==0) & !all(b.rp.alas[-c(1,2)]==0)) {
    bwd <- sd(scr.xb)/n^(1/6)
    pi.ric <- npreg(bws=bwd, txdat=cbind(scr.xb), tydat=Ti, exdat=cbind(scr.xb), ckerorder=4)$mean
  } else if(!all(b.ps.alas[-1]==0) & all(b.rp.alas[-c(1,2)]==0)) {
    bwd <- sd(scr.ps)/n^(1/6)
    pi.ric <- npreg(bws=bwd, txdat=cbind(scr.ps), tydat=Ti, exdat=cbind(scr.ps), ckerorder=4)$mean
  } else if(all(b.ps.alas[-1]==0) & all(b.rp.alas[-c(1,2)]==0)) {
    pi.ric <- mean(Ti)
  } else{
    bwd <- apply(cbind(scr.ps,scr.xb),2,sd)/n^(1/6)
    pi.ric <- npreg(bws=bwd, txdat=cbind(scr.ps, scr.xb), tydat=Ti, exdat=cbind(scr.ps, scr.xb), ckerorder=4)$mean
  }
  
  Ui.ric <- Ti/pi.ric - (1-Ti)/(1-pi.ric)
  
  if(offLog==T) {
    Zi <- cbind(
      bdBsMt(offLog(Xi[,which(lapply(apply(Xi,2,unique),length)>=3)])), Xi[,which(lapply(apply(Xi,2,unique),length)<3)],
      bdBsMt(offLog(Wi[,which(lapply(apply(Wi,2,unique),length)>=3)])), Wi[,which(lapply(apply(Wi,2,unique),length)<3)],
      Ti,
      Ui.ric
    )
  } else {
    if (is.null(dim(Wi))) {
      Zi <- cbind(nsbasis(Xi,6),Wi,Ti,Ui.ric)
    } else {
      Zi <- cbind(nsbasis(Xi,6),nsbasis(Wi,6),Ti,Ui.ric)
    }
  }
  bhat.rdg <- ridge(cbind(Yi[Ri==1], Z=Zi[Ri==1,]), fam, wgt=Gi[Ri==1])
  Yi.str <- g(cbind(1,Zi)%*%bhat.rdg)
  
  delta <- sum(Gi*(Ti*Yi.str/pi.ric))/sum(Gi*(Ti/pi.ric)) - sum(Gi*((1-Ti)*Yi.str/(1-pi.ric)))/sum(Gi*((1-Ti)/(1-pi.ric))) # normalized
}


# Simulated data (Kang and Schafer)
simulate = FALSE

if (simulate) {
  set.seed(100)
  
  n <- 5000
  
  p.Zi <- 4
  p.Wi <- 5
  mu.Zi <- rep(0,p.Zi)
  mu.Wi <- rep(0,p.Wi)
  Sigma.Zi <- diag(p.Zi)
  Sigma.Wi <- diag(p.Wi)
  
  alpha <- c(0,-1,.5,-.25,-.1)
  beta <- c(210,27.4,13.7,13.7,13.7)
  xi <- c(10,20,30,40,50)
  
  Zi <- mvrnorm(n,mu.Zi,Sigma.Zi)
  Wi <- mvrnorm(n,mu.Wi,Sigma.Wi)
  Yi <- cbind(1,Zi) %*% beta + Wi %*%xi + rnorm(n,0,1)
  
  pi <- expit(cbind(1,Zi) %*% alpha)
  
  Xi <- cbind(
    exp(Zi[,1]/2),
    Zi[,2]/(1+exp(Zi[,1])) + 10,
    (Zi[,1]*Zi[,3]/25 + .6)^3,
    (Zi[,2] + Zi[,4] + 20)^2
  )
  Xi <- Zi
  
  Ti <- rbinom(n,1,pi)
  Ri <- rbinom(n,1,.2)
  
  # Methods
  pihat.glm <- predict(glm(Ti~Xi,family=binomial()),type="response")
  muhat1.glm <- predict(glm(Yi~Xi+Ti),newdata=data.frame(Yi,Ti=1,Xi),type="response")
  muhat0.glm <- predict(glm(Yi~Xi+Ti),newdata=data.frame(Yi,Ti=0,Xi),type="response")
  xihat1.glm <- predict(glm(Yi~Wi+Xi+Ti),newdata=data.frame(Yi,Ti=1,Xi,Wi),type="response")
  xihat0.glm <- predict(glm(Yi~Wi+Xi+Ti),newdata=data.frame(Yi,Ti=0,Xi,Wi),type="response")
  
  
  pihat.dips <- dips(Yi,Ti,Xi,Gi=rep(1,n))
  nuhat <- mean(Ri)
  
  print(ipw(Yi,Ti,pihat.glm))
  print(aipw(Yi,Ti,pihat.glm,muhat1.glm,muhat0.glm))
  print(ipw(Yi,Ti,pihat.dips))
  
  print(aipw.ss(Yi,Ti,Ri,pihat.glm,nuhat,muhat1.glm,muhat0.glm,xihat1.glm,xihat0.glm))
  print(dips.ss(Yi,Ti,Ri,Wi,Xi,Gi=rep(1,n)))
  
  # Test
  Yi[Ri==0] <- rnorm(sum(Ri==0)) # alter the outcomes for unlabeled data
  print(dips.ss(Yi,Ti,Ri,Wi,Xi,Gi=rep(1,n)))
}

