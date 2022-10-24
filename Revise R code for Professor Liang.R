library(np)

###use multicore to booster our results
library(foreach)
library(doParallel)
library(iterators)
####




set.seed(123)

n <- 200

x1 <- cbind(rnorm(n,mean=1,sd=0.5),rnorm(n,mean=3,sd=1.5))
x2 <- cbind(rnorm(n,mean=1.5,sd=0.5),rnorm(n,mean=2,sd=1.5))

nx<- col(x1)

beta0 <- c(1.8,-3.2)
beta1 <- c(2,1.9)

e <- cbind(rnorm(n,mean=0,sd=1),rnorm(n,mean=0,sd=1.5))


#####need to add x1[,2] terms in building potential y0 and y1
y0 <- 12+5*x1[,1]^2*x1[,2]^3+x2[,1]*beta0[1]+x2[,2]*beta0[2]+e[,1]
y1 <- 2-3*x1[,1]^2+5*x1[,2]^3+x2[,1]*beta1[1]+x2[,2]*beta1[2]+e[,2]

z <- kronecker(rnorm(n,mean=0.5,sd=2),array(1,dim=c(1,1)))

u <- kronecker(rnorm(n,mean=0,sd=1),array(1,dim=c(1,1)))

gamma <- c(0.8,-2.1,2,1.5,3)

w <- (-1.5+x1[,1]*gamma[1]+x1[,2]*gamma[2]+x2[,1]*gamma[3]+x2[,2]*gamma[4]+z*gamma[5]>=u)

y <- y1*w+y0*(1-w)

kf<-function(v){
  apply(exp(-v^2/2)/sqrt(2*pi),1,prod)
}
##test other function
kf23<-function(v){
  apply(v^2,1,prod)
}
kf24<-function(v){
  apply(v,1,prod)
}
kf25<-function(v){
  apply(1/v,1,prod)
}

z24<-cbind(rnorm(5,mean=0,sd=1),rnorm(5,mean=1,sd=1))
v24<-kronecker(z24[1,1],array(1,dim=c(5,2)))

kf23(z24-v24)
kf25(z24-v24)


#####

ccdf1<-function(x1,y1,x2,y2,h1,h2){
  yr1<-NROW(y1)##nrow should be NROW()
  yr2<-NROW(y2)
  v1<-kronecker(x1,array(1,dim=c(yr1,1)))
  v2<-kronecker(x2,array(1,dim=c(yr2,1)))
  h1<-kronecker(h1,array(1,dim=dim(v1)))
  h2<-kronecker(h2,array(1,dim=dim(v2)))
  out<-sum(kf((y1-v1)/h1)*pnorm((v2-y2)/h2))/sum(kf((y1-v1)/h1))
  return(out)
}

nw<-length(w)

v<-rep(0,nw)
hw<-1.06*sd(w)/n^(.2)
hz<-1.06*sd(z)/n^(.2)
for (i in 1:nw){
  v[i]<-ccdf1(z[i,1],z,w[i,1],w,hz,hw)
}

##check v in np package


zfit23<- data.frame(w,z)
zfit23$w<-as.integer(zfit23$w)
w<-zfit23$w
z<-zfit23$z


bw23 <- npcdistbw(formula=w~z,bwmethod="normal-reference")
summary(bw23)
fitv23 <- fitted(npcdist(bws=bw23,newdata=zfit23))
zfit23$fitv23<-fitv23
zfit23$v<-v

##



G<-function(y0,x10,x20,w0,x1g,x2g,hx1g,hx2g,hyg,hvg){
  out<-rep(0,nw)
  for (i in 1:nw){
    xr1<-nrow(x1g)
    xr2<-nrow(x2g)
    v1<-kronecker(t(x10),array(1,dim=c(xr1,1)))
    v2<-kronecker(t(x20),array(1,dim=c(xr2,1)))
    vi<-kronecker(v[i],array(1,dim=c(length(v),1)))
    hx1n<-kronecker(t(hx1g),array(1,dim=c(xr1,1)))
    hx2n<-kronecker(t(hx2g),array(1,dim=c(xr2,1)))
    hvn<-kronecker(hvg,array(1,dim=dim(vi)))
    hyn<-kronecker(hyg,array(1,dim=dim(y)))
    out[i]<-sum(kf((v1-x1g)/hx1n)*kf((v2-x2g)/hx2n)*kf((vi-v)/hvn)*pnorm((y0-y)/hyn)*(w==w0))/sum(kf((v1-x1g)/hx1n)*kf((v2-x2g)/hx2n)*kf((vi-v)/hvn)*(w==w0))
  }
  return(mean(out))
}

x10m <- apply(x1,2,mean)
x20m <- apply(x2,2,mean)


hx1<-1.06*apply(x1,2,sd)/n^(.2)
hx2<-1.06*apply(x2,2,sd)/n^(.2)
hv<-1.06*sd(v)/n^(.2)
hy<-1.06*sd(y)/n^(.2)
tau <- 0.5

###
J <- seq(-1000,1000,by=0.01)
nJ <- length(J)
G0<-rep(0,nJ)
G1<-rep(0,nJ)



####### use uniroot function find y0
invgtreat<-function(q){
  uniroot(function(y0){G(y0,x10=x10m,x20=x20m,w0=1,x1g=x1,x2g=x2,hx1g=hx1,hx2g=hx2,hyg=hy,hvg=hv)-q},range(y))$root
}

invgtreat(tau)

invgcontrol<-function(q){
  uniroot(function(y0){G(y0,x10=x10m,x20=x20m,w0=0,x1g=x1,x2g=x2,hx1g=hx1,hx2g=hx2,hyg=hy,hvg=hv)-q},range(y))$root
}

invgcontrol(tau)


####rename J in min and max of y
J1 <- seq(min(y),max(y),by=0.01)
nJ1 <- length(J1)
G0<-rep(0,nJ1)
G1<-rep(0,nJ1)



for (j in 1:nJ1){
  G1[j] <- G(J1[j],x10=x10m,x20=x20m,w0=1,x1g=x1,x2g=x2,hx1g=hx1,hx2g=hx2,hyg=hy,hvg=hv)-tau
  if(abs(G1[j])<0.001) q1<-J1[j]
}


###use mutiple core to booster results.


registerDoParallel(8)  # use multicore, set to the number of our cores

foreach (j=1:nJ1, .combine=rbind) %dopar% {
  G1[j] <- G(J1[j],x10=x10m,x20=x20m,w0=1,x1g=x1,x2g=x2,hx1g=hx1,hx2g=hx2,hyg=hy,hvg=hv)-tau
  if(abs(G1[j])<0.001) q1<-J1[j]
  }



####true QTE
a<-12+5*0.9957148^2*3.0631817^3+1.8*1.515888-3.2*1.967172
b<-2-3*0.9957148^2+5*3.0631817^3+2*1.515888+1.9**1.967172

ctrue<-b-a

c<-139.82-142.1557

ctrue
c





#######


for (j in 1:nJ){
  G0[j] <- G(J[j],x10=x10m,x20=x20m,w0=0,x1g=x1,x2g=x2,hx1g=hx1,hx2g=hx2,hyg=hy,hvg=hv)-tau
  G1[j] <- G(J[j],x10=x10m,x20=x20m,w0=1,x1g=x1,x2g=x2,hx1g=hx1,hx2g=hx2,hyg=hy,hvg=hv)-tau
  if(G0[j]>0) q0<-J[j]
  if(G1[j]>0) q1<-J[j]
  if(G0[j]>0 & G1[j]>0) break
}


Treat <- q1-q0


#q0 <- J[min(which(G0>0))]
#q1 <- J[min(which(G1>0))]

#q0 <- uniroot(G,c(-1000,1000),tol=1e-6,x10=x10m,x20=x20m,w0=0,x1g=x1,x2g=x2,hx1g=hx1,hx2g=hx2,hyg=hy,hvg=hv)
#q1 <- uniroot(G,c(-1000,1000),tol=1e-6,x10=x10m,x20=x20m,w0=1,x1g=x1,x2g=x2,hx1g=hx1,hx2g=hx2,hyg=hy,hvg=hv)











































