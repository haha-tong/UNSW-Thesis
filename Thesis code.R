library(TwoCop)
library(MixSAL)
library(expm)
library(mvtnorm)
library(Matrix)
library(fMultivar)
library(sn)
library(TwoCop)
library(copula)
library(gclus)
library(LaplacesDemon)
library(ald)
library(scatterplot3d)
library(ggplot2)
library(VineCopula)
## Let's plot an Asymmetric Laplace Distribution!
##Density
library(ald)
par(mfrow = c(2,2))
loop.vector <- 1:4

sseq = seq(-40,80,0.5)
p = c(0.2,0.4)
mu = c(0,10)
dens = dALD(y=sseq,mu=0,sigma=sqrt(2),p=0.2)
dens2 = dALD(y=sseq,mu=10,sigma=sqrt(2),p=0.2)
dens3 = dALD(y=sseq,mu=0,sigma=sqrt(2),p=0.6)
dens4 = dALD(y=sseq,mu=10,sigma=sqrt(2),p=0.6)
plot(sseq,dens,type = "l",lwd=2,col="red",xlab="x",ylab="f(x)", main="ALD Density function(mu=0,tau=0.2)")
plot(sseq,dens2,type = "l",lwd=2,col="pink",xlab="x",ylab="f(x)", main="ALD Density function(mu=10,tau=0.2)")
plot(sseq,dens3,type = "l",lwd=2,col="blue",xlab="x",ylab="f(x)", main="ALD Density function(mu=0,tau=0.6)")
plot(sseq,dens4,type = "l",lwd=2,col="black",xlab="x",ylab="f(x)", main="ALD Density function(mu=10,tau=0.6)")
#Look that is an special case of the skewed family in Galarza (2016)
#with sigma_new = 2*sigma
require(lqr)
dens2 = dSKD(y = sseq,mu = 0,sigma = sqrt(2),p = 1,dist = "laplace")
points(sseq,dens2,pch="+",cex=0.75)
## Distribution Function
df = pALD(q=sseq,mu=0,sigma=3,p=0.75)
plot(sseq,df,type="l",lwd=2,col="blue",xlab="x",ylab="F(x)", main="ALD Distribution function")
abline(h=1,lty=2)
##Inverse Distribution Function
prob = seq(0,1,length.out = 1000)
idf = qALD(prob=prob,mu=50,sigma=3,p=0.75)
plot(prob,idf,type="l",lwd=2,col="gray30",xlab="x",ylab=expression(F^{-1}~(x)))
title(main="ALD Inverse Distribution function")
abline(v=c(0,1),lty=2)
#Random Sample Histogram
sample = rALD(n=10000,mu=50,sigma=3,p=0.75)
hist(sample,breaks = 70,freq = FALSE,ylim=c(0,max(dens)),main="")
title(main="Histogram and True density")
lines(sseq,dens,col="red",lwd=2)

library(LaplacesDemon)
x <- daml(c(1,2,3), c(0,1,2), diag(3))
X <- raml(1000, c(0,1,2), diag(3))
joint.density.plot(X[,1], X[,2], color=FALSE)
plot(X)

set.seed(1234)
B = 100
n = 1000
p = 2
k = 2 

Y = matrix(rnorm(n*p), n, p)
X = matrix(rnorm(n*k), n, k)
beta = matrix(rnorm(p*k), p, k)
tau = c(0.5,0.5)
Psi = matrix(c(1,0.5,
               0.4,1),p,p) 

sigma = matrix(NA, p, 1)
xi = matrix(NA, p, 1)
for (j in 1:p)
{
  sigma[j]= sqrt(2/(tau[j]*(1-tau[j])))
  xi[j]= (1-2*tau[j])/(tau[j]*(1-tau[j]))
}
sigma = as.vector(sigma)
Lambda = diag(diag(sigma%*%matrix(1, 1,p)))

nu = (2-p)/2

Sigma = Lambda %*% Psi %*% Lambda

delta = c(0.13,0.23)
D = diag(delta, nrow = length(delta))

eps = matrix(NA, nrow=n, ncol=p)

Y = array(NA,c(n,p,B))

for (b in 1:B)
{
  set.seed=1234
  W = rexp(n, rate=1)
  Z = t(rmvnorm(1,c(rep(0,p)),diag(p)))
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma) %*% Z) 
    
    Y[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}

plot(density(Y))
dens = density(Y)
polygon(density(Y), col="brown", border="black")
plot(Y)
Owen <- function(x, thre, gradiente)
{
  grad <- t(gradiente)
  par <- t(as.matrix(x))
  eps <- 1/thre
  
  z <- 1+par%*%grad
  ans <- z
  lo <- ((z)<eps)
  
  ans[ lo  ] <- log(eps) - 1.5 + 2*(z[lo])/eps - 0.5*(z[lo]/eps)^2
  ans[ !lo ] <- log( z[!lo] )
  -sum(ans)
}


## The function for the computation of empirical likelihood

EL <- function(SC)
{
  n <- NROW(SC)
  p <- NCOL(SC)
  
  ##find Lagrange multiplier
  if(p!=1)
  {
    OBJ <- optim(par=rep(0,p), fn=Owen, thre=n, gradiente=SC, control=list(maxit=1e4))
    molt_lagr_0 <- OBJ$pa
    conv <- OBJ$conv
  }
  else
  {
    molt_lagr_0 <- optim(par=rep(0,p), fn=Owen, thre=n, gradiente=SC, method="Brent", lower=-1e1, upper=1e1)$pa
  }
  
  ##weights
  w_emp_0 <- as.numeric( 1/(n*(1+molt_lagr_0%*%t(SC))) )
  
  if(p!=1)
  {
    list(w0=w_emp_0, conv=conv, elr=-2*sum(log(w_emp_0))-2*n*log(n), el=log(w_emp_0))
  }
  else
  {
    list(w0=w_emp_0, elr=-2*sum(log(w_emp_0))-2*n*log(n), el=log(w_emp_0) )
  }
}

ELCOP<- function(sample){
  n<-dim(sample)[1] 
  p<-dim(sample)[2]
  rr<-rank(sample[,1])
  ss<-rank(sample[,2])
  vv<-rr*ss
  #print(vv)
  for(b in 1:B){
    cor(uu[,1,b], uu[,2,b],method="spearman")
  }
  S=1000
  rh<-runif(S, -1,1)
  omega<-rep(0,S)
  for (s in 1:S) {
    #print(rh[s])
    estim= 12 * vv /(n^2-1) - 3*(n+1)/(n-1) - rh[s]
    omega[s]<-exp(-EL(estim)$elr)
    #print(omega[s])
  }
  psam<-sample(rh, size=S, rep=T, prob=omega)
  sintesi<-c(quantile(psam,.05), quantile(psam,.5),  quantile(psam,.95))
  return(as.numeric(sintesi))  
}


uu = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu[,,b] = pobs(Y[,,b])
}

List1 = list()
for(b in 1:B){
  res = as.numeric(ELCOP(uu[,,b]))
  List1[[length(List1)+1]] = res
  
}
res1 = unlist(List1)
res1 = matrix(res1, nrow = 100, byrow = T)
head(res1)

ELCOP.vs.freq<- function(sample,S=1000){
  
  n<-dim(sample)[1] 
  p<-dim(sample)[2]
  
  # Ranghi
  rr<-rank(sample[,1])
  ss<-rank(sample[,2])
  vv<-rr*ss
  
  # Frequentist estimate
  rhon=12 * sum(vv) /(n*(n+1)*(n-1)) - 3*(n+1)/(n-1)
  
  # Frequentist estimate of the variance
  An=sum(rr/(n+1)*ss/(n+1))/n
  Bn=sum((rr/(n+1))^2*(ss/(n+1))^2)/n
  
  Cn=0
  for(i in 1:length(rr))
  {
    const=rr[i]/(n+1)*ss[i]/(n+1)
    for(j in 1:length(rr))
    {
      for(k in 1:length(33))
      {
        Cn=Cn+const*(rr[k]<=rr[i])*(ss[k]<=ss[j])				
      }
    }
  }	
  Cn=Cn/n^3+0.25-An
  
  Dn=0
  for(i in 1:length(rr))
  {
    for(j in 1:length(ss))
    {
      Dn=Dn+ss[i]*ss[j]/(n+1)^2*max(rr[i]/(n+1),rr[j]/(n+1))
    }
  }
  Dn=Dn/n^2
  
  En=0
  for(i in 1:length(rr))
  {
    for(j in 1:length(ss))
    {
      En=En+rr[i]*rr[j]/(n+1)^2*max(ss[i]/(n+1),ss[j]/(n+1))
    }
  }
  En=En/n^2
  
  sigman=144*(-9*An^2+Bn+2*Cn+2*Dn+2*En)
  
  ICn=c()
  if(sigman>=0)
  {
    ICn=c(rhon-1.96*sqrt(sigman/n),rhon+1.96*sqrt(sigman/n))
  } else {
    ICn=rep(NA,2)
  }
  
  ### ABC 
  
  rh<-runif(S, -1,1)
  omega<-rep(0,S)
  for (s in 1:S) {
    estim= 12 * vv /(n^2-1) - 3*(n+1)/(n-1) - rh[s]
    omega[s]<-exp(-EL(estim)$elr)
  }
  psam<-sample(rh, size=S, rep=T, prob=omega)
  sintesi<-c(quantile(psam,.025), quantile(psam,.5),  quantile(psam,.975))
  return(list(rhon=rhon,ICn=ICn,ABC=sintesi,post=psam))  
}
List2 = list()
for(b in 1:B)
{
  res2 = ELCOP.vs.freq(uu[,,b])
  List2[[length(List2)+1]] = res2
}
res2 = unlist(List2)
res2 = matrix(res2, nrow = 100, byrow = T)
#colnames(res2) = c("rhon","ICn1","ICn2","2.5%","50%","97.5%")
head(res2)
List2a = list()
for(b in 1:B)
{
  postSMAL = List2[[b]]$post
  List2a[[length(List2a)+1]] = postSMAL
}
postSMAL = unlist(List2a)
postSMAL = matrix(postSMAL, nrow = 100, byrow = T)
plot(density(postSMAL))

ELCOP.vs.freq.NPTDC<- function(X,varU.freq,varL.freq)
{
  
  n=dim(X)[1]
  t=sqrt(n)
  p=dim(X)[2]
  
  ### Frequentist estimate
  uU.mat=matrix(rep((1-t/n),p),ncol=p)
  uL.mat=matrix(rep(t/n,p),ncol=p)
  lambda_U=2-((1-C.n(u=uU.mat,X=X))/(t/n))
  lambda_L=C.n(u=uL.mat,X=X)/(t/n)
  
  # Frequentist intervals
  lambdaU_IC=lambda_U+c(-1.96,1.96)*sqrt(varU.freq/n)
  lambdaL_IC=lambda_L+c(-1.96,1.96)*sqrt(varL.freq/n)
  
  ### ABC
  
  S=1000
  lambdaUh=runif(S, -1,1)
  lambdaLh=runif(S, -1,1)
  omegaU=rep(0,S)
  omegaL=rep(0,S)
  for (s in 1:S) {
    estimU= lambda_U - lambdaUh[s]
    estimL= lambda_L - lambdaLh[s]
    omegaU[s]=exp(-EL(estimU)$elr)
    omegaL[s]=exp(-EL(estimL)$elr)
  }
  psamU=sample(lambdaUh, size=S, rep=T, prob=omegaU)
  psamL=sample(lambdaLh, size=S, rep=T, prob=omegaL)
  outU=c(quantile(psamU,.025), quantile(psamU,.5),  quantile(psamU,.975))
  outL=c(quantile(psamL,.025), quantile(psamL,.5),  quantile(psamL,.975))
  return(list(lambda.freq=c(lambda_U,lambda_L),postU=psamU,postL=psamL))
}

bootvar.TDC<- function(sam,B=1000){
  
  n<-dim(sam)[1]
  t=sqrt(n)
  d<-dim(sam)[2]
  
  U.mat=matrix(rep((1-t/n),d),ncol=d)
  L.mat=matrix(rep((t/n),d),ncol=d)
  lambda_U=c()
  lambda_L=c()
  
  for(j in 1:B)
  {
    # upperTDC estimator
    idx=sample(1:n,size=n,rep=T)
    sam.boot=sam[idx,]
    lambda_U[j]=2-((1-C.n(u=U.mat,X=sam.boot))/(t/n))
    lambda_L[j]=C.n(u=L.mat,X=sam.boot)/(t/n)
    
  }
  return(cbind(lambda_U,lambda_L))  
  
}

List3 = list()

for(b in 1:B)
{
  res3 = bootvar.TDC(uu[,,b])
  List3[[length(List3)+1]] = res3
}
res3 = array(as.numeric(unlist(List3)), dim=c(1000, 2, 100))

List4 = list()
for(b in 1:B)
{
  
  res4 = var.TDC=apply(res3[,,b],2,var)
  List4[[length(List4)+1]] = res4
}
res4 = unlist(List4)
res4 = matrix(res4, nrow = 100, byrow = T)

List5 = list()
for(b in 1:B)
{
  res5 = ELCOP.vs.freq.NPTDC(uu[,,b],varU.freq=res4[b,1], varL.freq=res4[b,2])
  List5[[length(List5)+1]] = res5
}

head(List5[[1]])
List5a = list()
for(b in 1:B)
{
  postUSMAL = List5[[b]]$postU
  List5a[[length(List5a)+1]] = postUSMAL
}
postUSMAL = unlist(List5a)
postUSMAL = matrix(postUSMAL, nrow = 100, byrow = T)
plot(density(postUSMAL))

List5b = list()
for(b in 1:B)
{
  postLSMAL = List5[[b]]$postL
  List5b[[length(List5b)+1]] = postLSMAL
}
postLSMAL = unlist(List5b)
postLSMAL = matrix(postLSMAL, nrow = 100, byrow = T)





## Simulate from Clayton
cc=claytonCopula(2,dim=2)
persp(claytonCopula(2,dim=2),dCopula)
ssize=1000
uu.clayton=rCopula(ssize,cc)
scatterplot3d(uu.clayton)
uu2 <- uu.clayton
res6 = as.numeric(ELCOP(uu2))
res7 = ELCOP.vs.freq(uu2)
postCC = res7$post
plot(density(postCC))
res8  = bootvar.TDC(uu2)
res9  = var.TDC=apply(res8,2,var)
res10 =  ELCOP.vs.freq.NPTDC(uu2,varU.freq=res9[1], varL.freq=res9[2])
postUCC = res10$postU
postLCC = res10$postL
plot(density(postUCC))
plot(density(postLCC))
## Simulate from Normal copula
norm.cop = normalCopula(0.8,dim=2)
uu3 = rCopula(ssize,norm.cop)
scatterplot3d(uu3)
res11 = as.numeric(ELCOP(uu3))
res12 = ELCOP.vs.freq(uu3)
postnorm.cop = res12$post
plot(density(postnorm.cop))
res13  = bootvar.TDC(uu3)
res14  = var.TDC=apply(res13,2,var)
res15 =  ELCOP.vs.freq.NPTDC(uu3,varU.freq=res14[1], varL.freq=res14[2])
postUnorm.cop = res15$postU
postLnorm.cop = res15$postL
plot(density(postUnorm.cop))
plot(density(postLnorm.cop))
## Simulate from the Gumbel copula
gc<-gumbelCopula(2,dim=2)
persp(gc,dCopula)
uu4 = rCopula(ssize,gc)
scatterplot3d(uu4)
res16 = as.numeric(ELCOP(uu4))
res17 = ELCOP.vs.freq(uu4)
postGC = res17$post
plot(density(postGC))
res18  = bootvar.TDC(uu4)
res19  = var.TDC=apply(res18,2,var)
res20 =  ELCOP.vs.freq.NPTDC(uu4,varU.freq=res19[1], varL.freq=res19[2])
postUGC = res20$postU
postLGC = res20$postL
plot(density(postUGC))
plot(density(postLGC))


aa = seq(-1,1,0.5)
bb = seq(-1,1,0.5)
mm = expand.grid(aa,bb)
Psi_g = matrix(NA,p,p)
list_mat <- list()
for(i in 1:25){
  Psi_g = matrix(c(1,mm$Var1[i],
                    mm$Var2[i],1),p,p)
  list_mat[[i]] <- Psi_g
}

## Equal correlation
Psi_g1 = list_mat[[1]]
Psi_g2 = list_mat[[7]]
Psi_g3 = list_mat[[13]]
Psi_g4 = list_mat[[19]]
Psi_g5 = list_mat[[25]]

## Opposite correlation
Psi_ng1 = list_mat[[5]]
Psi_ng2 = list_mat[[9]]
Psi_ng3 = list_mat[[17]]
Psi_ng4 = list_mat[[21]]

Sigma_g1 = Lambda %*% Psi_g1 %*% Lambda
Y_g1 = array(NA,c(n,p,B))

for (b in 1:B)
{
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma_g1) %*% Z) 
    
    Y_g1[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}

Sigma_g2 = Lambda %*% Psi_g2 %*% Lambda
Y_g2 = array(NA,c(n,p,B))

for (b in 1:B)
{
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma_g2) %*% Z) 
    
    Y_g2[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}

Sigma_g3 = Lambda %*% Psi_g3 %*% Lambda
Y_g3 = array(NA,c(n,p,B))

for (b in 1:B)
{
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma_g3) %*% Z) 
    
    Y_g3[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}

Sigma_g4 = Lambda %*% Psi_g4 %*% Lambda
Y_g4 = array(NA,c(n,p,B))

for (b in 1:B)
{
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma_g4) %*% Z) 
    
    Y_g4[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}

Sigma_g5 = Lambda %*% Psi_g5 %*% Lambda
Y_g5 = array(NA,c(n,p,B))

for (b in 1:B)
{
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma_g5) %*% Z) 
    
    Y_g5[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}

Sigma_ng1 = Lambda %*% Psi_ng1 %*% Lambda
Y_ng1 = array(NA,c(n,p,B))

for (b in 1:B)
{
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma_ng1) %*% Z) 
    
    Y_ng1[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}  
  
Sigma_ng2 = Lambda %*% Psi_ng2 %*% Lambda
Y_ng2 = array(NA,c(n,p,B))

for (b in 1:B)
{
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma_ng2) %*% Z) 
    
    Y_ng2[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}  

Sigma_ng3 = Lambda %*% Psi_ng3 %*% Lambda
Y_ng3 = array(NA,c(n,p,B))

for (b in 1:B)
{
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma_ng3) %*% Z) 
    
    Y_ng3[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}  

Sigma_ng4 = Lambda %*% Psi_ng4 %*% Lambda
Y_ng4 = array(NA,c(n,p,B))

for (b in 1:B)
{
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma_ng4) %*% Z) 
    
    Y_ng4[i,,b]= beta %*% matrix(X[i,],ncol=1)+ eps[i,]
  }
}  

uu_g1 = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu_g1[,,b] = pobs(Y_g1[,,b])
}

uu_g2 = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu_g2[,,b] = pobs(Y_g2[,,b])
}

uu_g3 = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu_g3[,,b] = pobs(Y_g3[,,b])
}

uu_g4 = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu_g4[,,b] = pobs(Y_g4[,,b])
}

uu_g5 = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu_g5[,,b] = pobs(Y_g5[,,b])
}

uu_ng1 = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu_ng1[,,b] = pobs(Y_ng1[,,b])
}

uu_ng2 = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu_ng2[,,b] = pobs(Y_ng2[,,b])
}

uu_ng3 = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu_ng3[,,b] = pobs(Y_ng3[,,b])
}

uu_ng4 = array(NA,c(n,p,B))
for(b in 1:B)
{
  uu_ng4[,,b] = pobs(Y_ng4[,,b])
}

List2_g1 = list()
for(b in 1:B)
{
  res2_g1 = ELCOP.vs.freq(uu_g1[,,b])
  List2_g1[[length(List2_g1)+1]] = res2_g1
}
res2_g1 = unlist(List2_g1)
res2_g1 = matrix(res2_g1, nrow = 100, byrow = T)
List2a_g1 = list()
for(b in 1:B)
{
  postSMAL_g1 = List2_g1[[b]]$post
  List2a_g1[[length(List2a_g1)+1]] = postSMAL_g1
}
postSMAL_g1 = unlist(List2a_g1)
postSMAL_g1 = matrix(postSMAL_g1, nrow = 100, byrow = T)
rhom_g1 = mean(postSMAL_g1)
List3_g1 = list()

for(b in 1:B)
{
  res3_g1 = bootvar.TDC(uu_g1[,,b])
  List3_g1[[length(List3_g1)+1]] = res3_g1
}
res3_g1 = array(as.numeric(unlist(List3_g1)), dim=c(1000, 2, 100))

List4_g1 = list()
for(b in 1:B)
{
  
  res4_g1 = var.TDC=apply(res3_g1[,,b],2,var)
  List4_g1[[length(List4_g1)+1]] = res4_g1
}
res4_g1 = unlist(List4_g1)
res4_g1 = matrix(res4_g1, nrow = 100, byrow = T)

List5_g1 = list()
for(b in 1:B)
{
  res5_g1 = ELCOP.vs.freq.NPTDC(uu_g1[,,b],varU.freq=res4_g1[b,1], varL.freq=res4_g1[b,2])
  List5_g1[[length(List5_g1)+1]] = res5_g1
}

List5a_g1 = list()
for(b in 1:B)
{
  postUSMAL_g1 = List5_g1[[b]]$postU
  List5a_g1[[length(List5a_g1)+1]] = postUSMAL_g1
}
postUSMAL_g1 = unlist(List5a_g1)
postUSMAL_g1 = matrix(postUSMAL_g1, nrow = 100, byrow = T)
postUSMALm_g1 = mean(postUSMAL_g1)

List5b_g1 = list()
for(b in 1:B)
{
  postLSMAL_g1 = List5_g1[[b]]$postL
  List5b_g1[[length(List5b_g1)+1]] = postLSMAL_g1
}
postLSMAL_g1 = unlist(List5b_g1)
postLSMAL_g1 = matrix(postLSMAL_g1, nrow = 100, byrow = T)
postLSMALm_g1 = mean(postLSMAL_g1)
plot(density(postLSMAL_g1))

List2_g2 = list()
for(b in 1:B)
{
  res2_g2 = ELCOP.vs.freq(uu_g2[,,b])
  List2_g2[[length(List2_g2)+1]] = res2_g2
}
res2_g2 = unlist(List2_g2)
res2_g2 = matrix(res2_g2, nrow = 100, byrow = T)
List2a_g2 = list()
for(b in 1:B)
{
  postSMAL_g2 = List2_g2[[b]]$post
  List2a_g2[[length(List2a_g2)+1]] = postSMAL_g2
}
postSMAL_g2 = unlist(List2a_g2)
postSMAL_g2 = matrix(postSMAL_g2, nrow = 100, byrow = T)
rhom_g2 = mean(postSMAL_g2)
List3_g2 = list()

for(b in 1:B)
{
  res3_g2 = bootvar.TDC(uu_g2[,,b])
  List3_g2[[length(List3_g2)+1]] = res3_g2
}
res3_g2 = array(as.numeric(unlist(List3_g2)), dim=c(1000, 2, 100))

List4_g2 = list()
for(b in 1:B)
{
  
  res4_g2 = var.TDC=apply(res3_g2[,,b],2,var)
  List4_g2[[length(List4_g2)+1]] = res4_g2
}
res4_g2 = unlist(List4_g2)
res4_g2 = matrix(res4_g2, nrow = 100, byrow = T)

List5_g2 = list()
for(b in 1:B)
{
  res5_g2 = ELCOP.vs.freq.NPTDC(uu_g2[,,b],varU.freq=res4_g2[b,1], varL.freq=res4_g2[b,2])
  List5_g2[[length(List5_g2)+1]] = res5_g2
}

List5a_g2 = list()
for(b in 1:B)
{
  postUSMAL_g2 = List5_g2[[b]]$postU
  List5a_g2[[length(List5a_g2)+1]] = postUSMAL_g2
}
postUSMAL_g2 = unlist(List5a_g2)
postUSMAL_g2 = matrix(postUSMAL_g2, nrow = 100, byrow = T)
postUSMALm_g2 = mean(postUSMAL_g2)

List5b_g2 = list()
for(b in 1:B)
{
  postLSMAL_g2 = List5_g2[[b]]$postL
  List5b_g2[[length(List5b_g2)+1]] = postLSMAL_g2
}
postLSMAL_g2 = unlist(List5b_g2)
postLSMAL_g2 = matrix(postLSMAL_g2, nrow = 100, byrow = T)
postLSMALm_g2 = mean(postLSMAL_g2)
plot(density(postLSMAL_g2))

List2_g3 = list()
for(b in 1:B)
{
  res2_g3 = ELCOP.vs.freq(uu_g3[,,b])
  List2_g3[[length(List2_g3)+1]] = res2_g3
}
res2_g3 = unlist(List2_g3)
res2_g3 = matrix(res2_g3, nrow = 100, byrow = T)
List2a_g3 = list()
for(b in 1:B)
{
  postSMAL_g3 = List2_g3[[b]]$post
  List2a_g3[[length(List2a_g3)+1]] = postSMAL_g3
}
postSMAL_g3 = unlist(List2a_g3)
postSMAL_g3 = matrix(postSMAL_g3, nrow = 100, byrow = T)
rhom_g3 = mean(postSMAL_g3)
List3_g3 = list()

for(b in 1:B)
{
  res3_g3 = bootvar.TDC(uu_g3[,,b])
  List3_g3[[length(List3_g3)+1]] = res3_g3
}
res3_g3 = array(as.numeric(unlist(List3_g3)), dim=c(1000, 2, 100))

List4_g3 = list()
for(b in 1:B)
{
  
  res4_g3 = var.TDC=apply(res3_g3[,,b],2,var)
  List4_g3[[length(List4_g3)+1]] = res4_g3
}
res4_g3 = unlist(List4_g3)
res4_g3 = matrix(res4_g3, nrow = 100, byrow = T)

List5_g3 = list()
for(b in 1:B)
{
  res5_g3 = ELCOP.vs.freq.NPTDC(uu_g3[,,b],varU.freq=res4_g3[b,1], varL.freq=res4_g3[b,2])
  List5_g3[[length(List5_g3)+1]] = res5_g3
}

List5a_g3 = list()
for(b in 1:B)
{
  postUSMAL_g3 = List5_g3[[b]]$postU
  List5a_g3[[length(List5a_g3)+1]] = postUSMAL_g3
}
postUSMAL_g3 = unlist(List5a_g3)
postUSMAL_g3 = matrix(postUSMAL_g3, nrow = 100, byrow = T)
postUSMALm_g3 = mean(postUSMAL_g3)

List5b_g3 = list()
for(b in 1:B)
{
  postLSMAL_g3 = List5_g3[[b]]$postL
  List5b_g3[[length(List5b_g3)+1]] = postLSMAL_g3
}
postLSMAL_g3 = unlist(List5b_g3)
postLSMAL_g3 = matrix(postLSMAL_g3, nrow = 100, byrow = T)
postLSMALm_g3 = mean(postLSMAL_g3)
plot(density(postLSMAL_g3))

List2_g4 = list()
for(b in 1:B)
{
  res2_g4 = ELCOP.vs.freq(uu_g4[,,b])
  List2_g4[[length(List2_g4)+1]] = res2_g4
}
res2_g4 = unlist(List2_g4)
res2_g4 = matrix(res2_g4, nrow = 100, byrow = T)
List2a_g4 = list()
for(b in 1:B)
{
  postSMAL_g4 = List2_g4[[b]]$post
  List2a_g4[[length(List2a_g4)+1]] = postSMAL_g4
}
postSMAL_g4 = unlist(List2a_g4)
postSMAL_g4 = matrix(postSMAL_g4, nrow = 100, byrow = T)
rhom_g4 = mean(postSMAL_g4)
List3_g4 = list()

for(b in 1:B)
{
  res3_g4 = bootvar.TDC(uu_g4[,,b])
  List3_g4[[length(List3_g4)+1]] = res3_g4
}
res3_g4 = array(as.numeric(unlist(List3_g4)), dim=c(1000, 2, 100))

List4_g4 = list()
for(b in 1:B)
{
  
  res4_g4 = var.TDC=apply(res3_g4[,,b],2,var)
  List4_g4[[length(List4_g4)+1]] = res4_g4
}
res4_g4 = unlist(List4_g4)
res4_g4 = matrix(res4_g4, nrow = 100, byrow = T)

List5_g4 = list()
for(b in 1:B)
{
  res5_g4 = ELCOP.vs.freq.NPTDC(uu_g4[,,b],varU.freq=res4_g4[b,1], varL.freq=res4_g4[b,2])
  List5_g4[[length(List5_g4)+1]] = res5_g4
}

List5a_g4 = list()
for(b in 1:B)
{
  postUSMAL_g4 = List5_g4[[b]]$postU
  List5a_g4[[length(List5a_g4)+1]] = postUSMAL_g4
}
postUSMAL_g4 = unlist(List5a_g4)
postUSMAL_g4 = matrix(postUSMAL_g4, nrow = 100, byrow = T)
postUSMALm_g4 = mean(postUSMAL_g4)

List5b_g4 = list()
for(b in 1:B)
{
  postLSMAL_g4 = List5_g4[[b]]$postL
  List5b_g4[[length(List5b_g4)+1]] = postLSMAL_g4
}
postLSMAL_g4 = unlist(List5b_g4)
postLSMAL_g4 = matrix(postLSMAL_g4, nrow = 100, byrow = T)
postLSMALm_g4 = mean(postLSMAL_g4)
plot(density(postLSMAL_g4))

List2_g5 = list()
for(b in 1:B)
{
  res2_g5 = ELCOP.vs.freq(uu_g5[,,b])
  List2_g5[[length(List2_g5)+1]] = res2_g5
}
res2_g5 = unlist(List2_g5)
res2_g5 = matrix(res2_g5, nrow = 100, byrow = T)
List2a_g5 = list()
for(b in 1:B)
{
  postSMAL_g5 = List2_g5[[b]]$post
  List2a_g5[[length(List2a_g5)+1]] = postSMAL_g5
}
postSMAL_g5 = unlist(List2a_g5)
postSMAL_g5 = matrix(postSMAL_g5, nrow = 100, byrow = T)
rhom_g5 = mean(postSMAL_g5)
List3_g5 = list()

for(b in 1:B)
{
  res3_g5 = bootvar.TDC(uu_g5[,,b])
  List3_g5[[length(List3_g5)+1]] = res3_g5
}
res3_g5 = array(as.numeric(unlist(List3_g5)), dim=c(1000, 2, 100))

List4_g5 = list()
for(b in 1:B)
{
  
  res4_g5 = var.TDC=apply(res3_g5[,,b],2,var)
  List4_g5[[length(List4_g5)+1]] = res4_g5
}
res4_g5 = unlist(List4_g5)
res4_g5 = matrix(res4_g5, nrow = 100, byrow = T)

List5_g5 = list()
for(b in 1:B)
{
  res5_g5 = ELCOP.vs.freq.NPTDC(uu_g5[,,b],varU.freq=res4_g5[b,1], varL.freq=res4_g5[b,2])
  List5_g5[[length(List5_g5)+1]] = res5_g5
}

List5a_g5 = list()
for(b in 1:B)
{
  postUSMAL_g5 = List5_g5[[b]]$postU
  List5a_g5[[length(List5a_g5)+1]] = postUSMAL_g5
}
postUSMAL_g5 = unlist(List5a_g5)
postUSMAL_g5 = matrix(postUSMAL_g5, nrow = 100, byrow = T)
postUSMALm_g5 = mean(postUSMAL_g5)

List5b_g5 = list()
for(b in 1:B)
{
  postLSMAL_g5 = List5_g5[[b]]$postL
  List5b_g5[[length(List5b_g5)+1]] = postLSMAL_g5
}
postLSMAL_g5 = unlist(List5b_g5)
postLSMAL_g5 = matrix(postLSMAL_g5, nrow = 100, byrow = T)
postLSMALm_g5 = mean(postLSMAL_g5)

rhom_g = cbind(rhom_g1,rhom_g2,rhom_g3,rhom_g4,rhom_g5)
postUSMALm_g = cbind(postUSMALm_g1,postUSMALm_g2,postUSMALm_g3,postUSMALm_g4,postUSMALm_g5)
postLSMALm_g = cbind(postLSMALm_g1,postLSMALm_g2,postLSMALm_g3,postLSMALm_g4,postLSMALm_g5)

List2_ng1 = list()
for(b in 1:B)
{
  res2_ng1 = ELCOP.vs.freq(uu_ng1[,,b])
  List2_ng1[[length(List2_ng1)+1]] = res2_ng1
}
res2_ng1 = unlist(List2_ng1)
res2_ng1 = matrix(res2_ng1, nrow = 100, byrow = T)
List2a_ng1 = list()
for(b in 1:B)
{
  postSMAL_ng1 = List2_ng1[[b]]$post
  List2a_ng1[[length(List2a_ng1)+1]] = postSMAL_ng1
}
postSMAL_ng1 = unlist(List2a_ng1)
postSMAL_ng1 = matrix(postSMAL_ng1, nrow = 100, byrow = T)
rhom_ng1 = mean(postSMAL_ng1)
List3_ng1 = list()

for(b in 1:B)
{
  res3_ng1 = bootvar.TDC(uu_ng1[,,b])
  List3_ng1[[length(List3_ng1)+1]] = res3_ng1
}
res3_ng1 = array(as.numeric(unlist(List3_ng1)), dim=c(1000, 2, 100))

List4_ng1 = list()
for(b in 1:B)
{
  
  res4_ng1 = var.TDC=apply(res3_ng1[,,b],2,var)
  List4_ng1[[length(List4_ng1)+1]] = res4_ng1
}
res4_ng1 = unlist(List4_ng1)
res4_ng1 = matrix(res4_ng1, nrow = 100, byrow = T)

List5_ng1 = list()
for(b in 1:B)
{
  res5_ng1 = ELCOP.vs.freq.NPTDC(uu_ng1[,,b],varU.freq=res4_ng1[b,1], varL.freq=res4_ng1[b,2])
  List5_ng1[[length(List5_ng1)+1]] = res5_ng1
}

List5a_ng1 = list()
for(b in 1:B)
{
  postUSMAL_ng1 = List5_ng1[[b]]$postU
  List5a_ng1[[length(List5a_ng1)+1]] = postUSMAL_ng1
}
postUSMAL_ng1 = unlist(List5a_ng1)
postUSMAL_ng1 = matrix(postUSMAL_ng1, nrow = 100, byrow = T)
postUSMALm_ng1 = mean(postUSMAL_ng1)

List5b_ng1 = list()
for(b in 1:B)
{
  postLSMAL_ng1 = List5_ng1[[b]]$postL
  List5b_ng1[[length(List5b_ng1)+1]] = postLSMAL_ng1
}
postLSMAL_ng1 = unlist(List5b_ng1)
postLSMAL_ng1 = matrix(postLSMAL_ng1, nrow = 100, byrow = T)
postLSMALm_ng1 = mean(postLSMAL_ng1)
plot(density(postLSMAL_ng1))

List2_ng2 = list()
for(b in 1:B)
{
  res2_ng2 = ELCOP.vs.freq(uu_ng2[,,b])
  List2_ng2[[length(List2_ng2)+1]] = res2_ng2
}
res2_ng2 = unlist(List2_ng2)
res2_ng2 = matrix(res2_ng2, nrow = 100, byrow = T)
List2a_ng2 = list()
for(b in 1:B)
{
  postSMAL_ng2 = List2_ng2[[b]]$post
  List2a_ng2[[length(List2a_ng2)+1]] = postSMAL_ng2
}
postSMAL_ng2 = unlist(List2a_ng2)
postSMAL_ng2 = matrix(postSMAL_ng2, nrow = 100, byrow = T)
rhom_ng2 = mean(postSMAL_ng2)
List3_ng2 = list()

for(b in 1:B)
{
  res3_ng2 = bootvar.TDC(uu_ng2[,,b])
  List3_ng2[[length(List3_ng2)+1]] = res3_ng2
}
res3_ng2 = array(as.numeric(unlist(List3_ng2)), dim=c(1000, 2, 100))

List4_ng2 = list()
for(b in 1:B)
{
  
  res4_ng2 = var.TDC=apply(res3_ng2[,,b],2,var)
  List4_ng2[[length(List4_ng2)+1]] = res4_ng2
}
res4_ng2 = unlist(List4_ng2)
res4_ng2 = matrix(res4_ng2, nrow = 100, byrow = T)

List5_ng2 = list()
for(b in 1:B)
{
  res5_ng2 = ELCOP.vs.freq.NPTDC(uu_ng2[,,b],varU.freq=res4_ng2[b,1], varL.freq=res4_ng2[b,2])
  List5_ng2[[length(List5_ng2)+1]] = res5_ng2
}

List5a_ng2 = list()
for(b in 1:B)
{
  postUSMAL_ng2 = List5_ng2[[b]]$postU
  List5a_ng2[[length(List5a_ng2)+1]] = postUSMAL_ng2
}
postUSMAL_ng2 = unlist(List5a_ng2)
postUSMAL_ng2 = matrix(postUSMAL_ng2, nrow = 100, byrow = T)
postUSMALm_ng2 = mean(postUSMAL_ng2)

List5b_ng2 = list()
for(b in 1:B)
{
  postLSMAL_ng2 = List5_ng2[[b]]$postL
  List5b_ng2[[length(List5b_ng2)+1]] = postLSMAL_ng2
}
postLSMAL_ng2 = unlist(List5b_ng2)
postLSMAL_ng2 = matrix(postLSMAL_ng2, nrow = 100, byrow = T)
postLSMALm_ng2 = mean(postLSMAL_ng2)
plot(density(postLSMAL_ng2))

List2_ng3 = list()
for(b in 1:B)
{
  res2_ng3 = ELCOP.vs.freq(uu_ng3[,,b])
  List2_ng3[[length(List2_ng3)+1]] = res2_ng3
}
res2_ng3 = unlist(List2_ng3)
res2_ng3 = matrix(res2_ng3, nrow = 100, byrow = T)
List2a_ng3 = list()
for(b in 1:B)
{
  postSMAL_ng3 = List2_ng3[[b]]$post
  List2a_ng3[[length(List2a_ng3)+1]] = postSMAL_ng3
}
postSMAL_ng3 = unlist(List2a_ng3)
postSMAL_ng3 = matrix(postSMAL_ng3, nrow = 100, byrow = T)
rhom_ng3 = mean(postSMAL_ng3)
List3_ng3 = list()

for(b in 1:B)
{
  res3_ng3 = bootvar.TDC(uu_ng3[,,b])
  List3_ng3[[length(List3_ng3)+1]] = res3_ng3
}
res3_ng3 = array(as.numeric(unlist(List3_ng3)), dim=c(1000, 2, 100))

List4_ng3 = list()
for(b in 1:B)
{
  
  res4_ng3 = var.TDC=apply(res3_ng3[,,b],2,var)
  List4_ng3[[length(List4_ng3)+1]] = res4_ng3
}
res4_ng3 = unlist(List4_ng3)
res4_ng3 = matrix(res4_ng3, nrow = 100, byrow = T)

List5_ng3 = list()
for(b in 1:B)
{
  res5_ng3 = ELCOP.vs.freq.NPTDC(uu_ng3[,,b],varU.freq=res4_ng3[b,1], varL.freq=res4_ng3[b,2])
  List5_ng3[[length(List5_ng3)+1]] = res5_ng3
}

List5a_ng3 = list()
for(b in 1:B)
{
  postUSMAL_ng3 = List5_ng3[[b]]$postU
  List5a_ng3[[length(List5a_ng3)+1]] = postUSMAL_ng3
}
postUSMAL_ng3 = unlist(List5a_ng3)
postUSMAL_ng3 = matrix(postUSMAL_ng3, nrow = 100, byrow = T)
postUSMALm_ng3 = mean(postUSMAL_ng3)

List5b_ng3 = list()
for(b in 1:B)
{
  postLSMAL_ng3 = List5_ng3[[b]]$postL
  List5b_ng3[[length(List5b_ng3)+1]] = postLSMAL_ng3
}
postLSMAL_ng3 = unlist(List5b_ng3)
postLSMAL_ng3 = matrix(postLSMAL_ng3, nrow = 100, byrow = T)
postLSMALm_ng3 = mean(postLSMAL_ng3)
plot(density(postLSMAL_ng3))

List2_ng4 = list()
for(b in 1:B)
{
  res2_ng4 = ELCOP.vs.freq(uu_ng4[,,b])
  List2_ng4[[length(List2_ng4)+1]] = res2_ng4
}
res2_ng4 = unlist(List2_ng4)
res2_ng4 = matrix(res2_ng4, nrow = 100, byrow = T)
List2a_ng4 = list()
for(b in 1:B)
{
  postSMAL_ng4 = List2_ng4[[b]]$post
  List2a_ng4[[length(List2a_ng4)+1]] = postSMAL_ng4
}
postSMAL_ng4 = unlist(List2a_ng4)
postSMAL_ng4 = matrix(postSMAL_ng4, nrow = 100, byrow = T)
rhom_ng4 = mean(postSMAL_ng4)
List3_ng4 = list()

for(b in 1:B)
{
  res3_ng4 = bootvar.TDC(uu_ng4[,,b])
  List3_ng4[[length(List3_ng4)+1]] = res3_ng4
}
res3_ng4 = array(as.numeric(unlist(List3_ng4)), dim=c(1000, 2, 100))

List4_ng4 = list()
for(b in 1:B)
{
  
  res4_ng4 = var.TDC=apply(res3_ng4[,,b],2,var)
  List4_ng4[[length(List4_ng4)+1]] = res4_ng4
}
res4_ng4 = unlist(List4_ng4)
res4_ng4 = matrix(res4_ng4, nrow = 100, byrow = T)

List5_ng4 = list()
for(b in 1:B)
{
  res5_ng4 = ELCOP.vs.freq.NPTDC(uu_ng4[,,b],varU.freq=res4_ng4[b,1], varL.freq=res4_ng4[b,2])
  List5_ng4[[length(List5_ng4)+1]] = res5_ng4
}

List5a_ng4 = list()
for(b in 1:B)
{
  postUSMAL_ng4 = List5_ng4[[b]]$postU
  List5a_ng4[[length(List5a_ng4)+1]] = postUSMAL_ng4
}
postUSMAL_ng4 = unlist(List5a_ng4)
postUSMAL_ng4 = matrix(postUSMAL_ng4, nrow = 100, byrow = T)
postUSMALm_ng4 = mean(postUSMAL_ng4)

List5b_ng4 = list()
for(b in 1:B)
{
  postLSMAL_ng4 = List5_ng4[[b]]$postL
  List5b_ng4[[length(List5b_ng4)+1]] = postLSMAL_ng4
}
postLSMAL_ng4 = unlist(List5b_ng4)
postLSMAL_ng4 = matrix(postLSMAL_ng4, nrow = 100, byrow = T)
postLSMALm_ng4 = mean(postLSMAL_ng4)
plot(density(postLSMAL_ng4))
rhom_ng = cbind(rhom_ng1,rhom_ng2,rhom_ng3,rhom_ng4)
postUSMALm_ng = cbind(postUSMALm_ng1,postUSMALm_ng2,postUSMALm_ng3,postUSMALm_ng4)
postLSMALm_ng = cbind(postLSMALm_ng1,postLSMALm_ng2,postLSMALm_ng3,postLSMALm_ng4)
rhom_ng
postUSMALm_ng
postLSMALm_ng
##Simulate from Smal with k=4
set.seed(1234)


Ytest = matrix(rnorm(n*p), n, p)
Xtest = matrix(rnorm(n*k), n, 4)
betatest = matrix(rnorm(p*k), p, 4)

Ytest = array(NA,c(n,p,B))

for (b in 1:B)
{
  set.seed=1234
  W = rexp(n, rate=1)
  Z = t(rmvnorm(1,c(rep(0,p)),diag(p)))
  for (i in 1:n)
  {
    
    eps[i,]=as.matrix(D %*% xi * W[i]) + sqrt(W[i]) * (D %*% sqrtm(Sigma) %*% Z) 
    
    Ytest[i,,b]= betatest %*% matrix(Xtest[i,],ncol=1)+ eps[i,]
  }
}
uutest = array(NA,c(n,p,B))
for(b in 1:B)
{
  uutest[,,b] = pobs(Ytest[,,b])
}
List2test = list()
for(b in 1:B)
{
  res2test = ELCOP.vs.freq(uutest[,,b])
  List2test[[length(List2test)+1]] = res2test
}
res2test = unlist(List2test)
res2test = matrix(res2test, nrow = 100, byrow = T)
List2atest = list()
for(b in 1:B)
{
  postSMALtest = List2test[[b]]$post
  List2atest[[length(List2atest)+1]] = postSMALtest
}
postSMALtest = unlist(List2atest)
postSMALtest = matrix(postSMALtest, nrow = 100, byrow = T)
plot(density(postSMALtest))
List3test = list()

for(b in 1:B)
{
  res3test = bootvar.TDC(uutest[,,b])
  List3test[[length(List3test)+1]] = res3test
}
res3test = array(as.numeric(unlist(List3test)), dim=c(1000, 2, 100))

List4test = list()
for(b in 1:B)
{
  
  res4test = var.TDC=apply(res3test[,,b],2,var)
  List4test[[length(List4test)+1]] = res4test
}
res4test = unlist(List4test)
res4test = matrix(res4test, nrow = 100, byrow = T)

List5test = list()
for(b in 1:B)
{
  res5test = ELCOP.vs.freq.NPTDC(uutest[,,b],varU.freq=res4test[b,1], varL.freq=res4test[b,2])
  List5test[[length(List5test)+1]] = res5test
}

List5atest = list()
for(b in 1:B)
{
  postUSMALtest = List5test[[b]]$postU
  List5atest[[length(List5atest)+1]] = postUSMALtest
}
postUSMALtest = unlist(List5atest)
postUSMALtest = matrix(postUSMALtest, nrow = 100, byrow = T)
plot(density(postUSMALtest))

List5btest = list()
for(b in 1:B)
{
  postLSMALtest = List5test[[b]]$postL
  List5btest[[length(List5btest)+1]] = postLSMALtest
}
postLSMALtest = unlist(List5btest)
postLSMALtest = matrix(postLSMALtest, nrow = 100, byrow = T)
postLSMALmtest = mean(postLSMALtest)
plot(density(postLSMALtest))
## Real data
bd = data("body")
Y1 = body$CalfG
Y2 = body$ThighG
Ybody = cbind(Y1,Y2)
BMI = body$Weight/(body$Height)^2
Xbody = cbind(body$Weight,body$Age,BMI,body$Height)
betabody = matrix(rnorm(p*k), 2, 4)
uubody = pobs(Ybody)
res2body = ELCOP.vs.freq(uubody)
postbody = res2body$post
plot(density(postbody))
res3body = bootvar.TDC(uubody)
res4body = var.TDC=apply(res3body,2,var)
res5body = ELCOP.vs.freq.NPTDC(uubody,varU.freq=res4body[1], varL.freq=res4body[2])
postUbody = res5body$postU
postLbody = res5body$postL
plot(density(postUbody))
plot(density(postLbody))


##########################plots
plot(density(postSMAL_g1))
plot(density(postSMAL_g2))
plot(density(postSMAL_g3))
plot(density(postSMAL_g4))
plot(density(postSMAL_g5))
plot(density(postSMAL_ng1))
plot(density(postSMAL_ng2))
plot(density(postSMAL_ng3))
plot(density(postSMAL_ng4))
plot(density(postUSMAL_g1))
plot(density(postUSMAL_g2))
plot(density(postUSMAL_g3))
plot(density(postUSMAL_g4))
plot(density(postUSMAL_g5))
plot(density(postUSMAL_ng1))
plot(density(postUSMAL_ng2))
plot(density(postUSMAL_ng3))
plot(density(postUSMAL_ng4))
plot(density(postLSMAL_g1))
plot(density(postLSMAL_g2))
plot(density(postLSMAL_g3))
plot(density(postLSMAL_g4))
plot(density(postLSMAL_g5))
plot(density(postLSMAL_ng1))
plot(density(postLSMAL_ng2))
plot(density(postLSMAL_ng3))
plot(density(postLSMAL_ng4))

################################



ssss = BiCopPar2TailDep(1, 0.7)
BiCop(1, 0.7)$taildep
plot(ssss)
library(evd)
op <- par(mfrow = c(1, 3))
BiCopChiPlot(uu.clayton[,1], uu.clayton[,2],main="General chi-plot")
plot(ssss)
BiCopChiPlot(uu.clayton[,1], uu.clayton[,2], mode = "lower", xlim = c(-1,1),
             ylim = c(-1,1), main = "Lower chi-plot")
BiCopChiPlot(uu.clayton[,1], uu.clayton[,2], mode = "upper", xlim = c(-1,1),
             ylim = c(-1,1), main = "Upper chi-plot")

