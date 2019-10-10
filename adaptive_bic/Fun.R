###################Detialed regression functions##################################


g.logit = function(xx){exp(xx)/(exp(xx)+1)}
logit = function(xx){log(xx/(1-xx))}
dg.logit = function(xx){exp(xx)/(exp(xx)+1)^2}

logitlik.fun = function(bet.mat,dat){
  yi = dat[,1]; xi = dat[,-1]; pi.mat = g.logit(cbind(1,xi)%*%bet.mat) ## N x B
  apply(log(pi.mat)*yi + log(1-pi.mat)*(1-yi),2,sum)
}

A.fun = function(bet,dat){
  yy = dat[,1]; xx.vec = cbind(1,dat[,-1])
  -t(c(dg.logit(xx.vec%*%bet))*xx.vec)%*%xx.vec/length(yy)
}

###funtion used for ordinal approximate lasso
A.ridge.fun = function(par, dat,lambda, K){
  theta=par[1:(K-1)]
  eta=par[K:length(par)]
  eta=matrix(eta,nrow=length(eta),1)
  yy = dat[,1]; xx.vec = data.matrix(dat[,-1])
  n=length(yy)
  xeta=c(data.matrix(xx.vec)%*%eta)
  K=length(unique(yy))
  
  b.dd=matrix(0, nrow=length(c(theta,eta)), ncol=length(c(theta,eta)))
  b.dd[1,1]=sum(-(1*(yy==2)+1*(yy==1))*c(dg.logit(theta[1]+xeta))-(yy==2)*(exp(theta[1]+theta[2])/(exp(theta[2])-exp(theta[1]))^2))
  b.dd[1,2]=sum((yy==2)*(exp(theta[1]+theta[2])/(exp(theta[2])-exp(theta[1]))^2))
  b.dd[1,(K:dim(b.dd)[2])]=b.dd[(K:dim(b.dd)[2]),1]=colSums(-(1*(yy==2)+1*(yy==1))*diag(c(dg.logit(theta[1]+xeta)))%*%xx.vec)
  
  theta2.2=NULL
  for(k in 2:(K-2)){
    tmp=sum(-(1*(yy==(k+1))+1*(yy==(k)))*c(dg.logit(theta[k]+xeta))-
              (yy==k)*exp(theta[k]+theta[k-1])/(exp(theta[k])-exp(theta[k-1]))^2-
              (yy==(k+1))*exp(theta[k]+theta[k+1])/(exp(theta[k])-exp(theta[k+1]))^2)
    theta2.2=c(theta2.2, tmp)
  }
  diag(b.dd)[2:(K-2)]=theta2.2
  
  for(k in 2:(K-1)){
    tmp=sum((yy==(k))*exp(theta[k]+theta[k-1])/(exp(theta[k])-exp(theta[k-1]))^2)
    b.dd[k, k-1]=b.dd[k-1,k]=tmp
  }
  
  b.dd[K-1, K-1]=sum(-(yy==(K-1))*exp(theta[K-1]+theta[K-2])/(exp(theta[K-1])-exp(theta[K-2]))^2-
                       (1*(yy==K)+1*(yy==(K-1)))*c(dg.logit(theta[K-1]+xeta)))
  
  for(k in 2:(K-2)){
    tmp=colSums(-(1*(yy==(k+1))+1*(yy==k))*diag(c(dg.logit(theta[k]+xeta)))%*%xx.vec)
    b.dd[K:dim(b.dd)[1], k]=b.dd[k,K:dim(b.dd)[2]]=tmp
  } 
  
  b.dd[K:dim(b.dd)[1], K-1]=b.dd[K-1,K:dim(b.dd)[2]]=colSums((-(1*(yy==K)+1*(yy==(K-1)))*diag(c(dg.logit(theta[K-1]+xeta)))%*%xx.vec))
  
  for(i in 1:length(eta)){
    for(j in 1:i){
      tmp1=sum((yy==1)*(-xx.vec[,i]*xx.vec[,j]*c(dg.logit(theta[1]+xeta))))
      tmp2=0
      for(k in 2:(K-1)){
        tmp2=tmp2+sum((yy==k)*(-xx.vec[,i]*xx.vec[,j]*c(dg.logit(theta[k-1]+xeta))-xx.vec[,i]*xx.vec[,j]*c(dg.logit(theta[k]+xeta))))
      }
      tmp3=sum((yy==K)*(-xx.vec[,i]*xx.vec[,j]*c(dg.logit(theta[K-1]+xeta))))
      if(i==j){
        b.dd[K-1+i,K-1+j]=b.dd[K-1+j, K-1+i]=sum(tmp1, tmp2, tmp3)-n*lambda}
      if(i!=j){
        b.dd[K-1+i,K-1+j]=b.dd[K-1+j, K-1+i]=sum(tmp1, tmp2, tmp3)}
    }
  }
  b.dd
}

loglik.mat=function(par.mat, Y, X, K){
  logliki.mat=function(par.mat, Yi, Xi, K){
    theta.mat=par.mat[1:(K-1),]
    eta.mat=par.mat[K:dim(par.mat)[1],]
    L1=(Yi==1)*(theta.mat[1,]+Xi%*%eta.mat-log(1+exp(theta.mat[1,]+Xi%*%eta.mat)))
    L2=0
    for (k in 2:(K-1)){
      tmp=(Yi==k)*(Xi%*%eta.mat+log(exp(theta.mat[k,])-exp(theta.mat[k-1,]))-log(1+exp(theta.mat[k,]+Xi%*%eta.mat))-log(1+exp(theta.mat[k-1,]+Xi%*%eta.mat)))
      L2=L2+tmp
    }
    L3=-(Yi==K)*log(1+exp(theta.mat[K-1,]+Xi%*%eta.mat))
    loglik=as.vector((L1+L2+L3))
    loglik}
  colSums(matrix(unlist(lapply(1:length(Y), function(i) logliki.mat(par.mat, Y[i], X[i,], K))), ncol=dim(par.mat)[2], byrow=T))
}

Est.ALASSO.GLM.Approx = function(data,Wi=NULL,rtn="EST",adap=T,BIC.factor=0.1,fam0,offset=NULL,N_groups){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); pp = ncol(x); if(is.null(Wi)){Wi=rep(1,nn)};
  bini = glmnet(x,y,weights=Wi,family=fam0,alpha=0,offset=offset)
  lam.xx=svd(t(x)%*%x/nn)$d; tmpdf = apply(lam.xx/(lam.xx+VTM(bini$lambda,pp)),2,sum)
  tmpind = which.min(deviance(bini)+2*tmpdf); 
  lam.ridge = bini$lambda[tmpind]; c(bini$a0[tmpind], bini$beta[,tmpind]); 
  bini = c(bini$a0[tmpind], bini$beta[,tmpind]); 
  if(adap){
    w.b = 1/abs(bini[-1])
  }
  else {
    w.b = rep(1,pp)
  }
  
  Ahalf = svd(-A.fun(bini,data)+diag(c(0,rep(1,pp)))); Ahalf = Ahalf$u%*%diag(sqrt(Ahalf$d))%*%t(Ahalf$v)
  ynew = Ahalf%*%bini; xnew = Ahalf
  tmpfit = glmnet(x=xnew, y=ynew, family='gaussian',penalty.factor = c(0,w.b),
                  alpha=1, lambda = 10^seq(-4,3,0.01), intercept=F)
  BIC.lam = -2*logitlik.fun(tmpfit$beta,data)/nn+min(N_groups^BIC.factor,log(N_groups))*tmpfit$df/N_groups
  m.opt = which.min(BIC.lam); bhat.modBIC = tmpfit$beta[,m.opt]; lamhat = tmpfit$lambda[m.opt]
  bhat.modBIC
}

Est.ALASSO.POR.Approx = function(Y, X, K, BIC.factor=0.1){
  nn=length(Y);
  dat=cbind(Y, X)
  lam.ridge=dim(X)[2]/dim(X)[1]
  #lam.xx=svd(t(X)%*%X/nn)$d;  lam.99 = min(lam.xx[cumsum(lam.xx)/sum(lam.xx) < 0.99])
  fit.ini <- ordinalNet(X, Y, family="cumulative", link="logit",alpha=0,
                        parallelTerms=TRUE, nonparallelTerms=FALSE, lambdaVals=lam.ridge)
  bini=coef(fit.ini)
  w.eta = 1/abs(bini[K:length(bini)])
  Ahalf = svd(-A.ridge.fun(bini,dat, lam.ridge, K)+diag(c(rep(0,(K-1)), rep(1, length(w.eta)))))
  
  Ahalf = Ahalf$u%*%diag(sqrt(Ahalf$d))%*%t(Ahalf$v)
  ynew = Ahalf%*%bini; xnew = Ahalf
  tmpfit = glmnet(x=xnew, y=ynew, family='gaussian',penalty.factor = c(rep(0,(K-1)), w.eta),
                  alpha=1, lambda = 10^seq(-4,3,0.01), intercept=F)
  
  BIC.lam=-2*loglik.mat(tmpfit$beta,Y, X, K)+min(nn^BIC.factor,log(nn))*tmpfit$df
  m.opt = which.min(BIC.lam); bhat.modBIC = tmpfit$beta[,m.opt]; lamhat = tmpfit$lambda[m.opt]
  bhat.modBIC
}


# ## Simulate the data
# NCols= 130
# NRows= 2500
# 
# X = matrix(runif(NCols*NRows), ncol=NCols)
# bet = rep(0, NCols); bet[sample(1:NCols, 10)] = 1
# K = 5
# prbs = c(210, 346, 654, 550, 240)
# prbs = prbs/sum(prbs)
# Y = Y_true = factor(replicate(NRows, sample(0:(K-1), 1, prob=prbs)),
#                     ordered=TRUE, levels=0:(K-1))
# intercepts = polr(Y~X %*% bet)$zeta
# Y[sample(1:NRows, 10)] = NA
# 
