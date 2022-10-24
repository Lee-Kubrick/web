


library(np)

set.seed(123)

n <- 200


z<- data.frame(y=c(rnorm(n,3,5)), x1=c(rnorm(n,mean=-2,sd=0.5)), x2=c(rbinom(n,1,0.75)))

               
ckf<-function(v){
  exp(-v^2/2)/sqrt(2*pi)
}

#yf<-function(u){
#  integrate(ckf, -Inf, u)$value
#}

dkf <-function(xi,x,lamd){
  res<-ifelse(xi!=x,lamd,1)
  return(res)
}

condist<-function(ys,y,hy,xs1,x1,hx1,xs2,x2,lamdx2){
res<-
  sum(pnorm((ys-y)/hy)*ckf((x1-xs1)/hx1)/hx1*
  dkf(xs2,x2,lamdx2))/
  (sum(ckf((x1-xs1)/hx1)/hx1*dkf(xs2,x2,lamdx2)))
return(res)
}


##choose bandwith from np package
y<-z$y
x1<- z$x1
x2<-z$x2

bw <- npcdistbw(formula=y~x1+ordered(x2))
summary(bw)
fitv <- fitted(npcdist(bws=bw,newdata=z))


bw1 <- npcdistbw(formula=y~x1+ordered(x2),bwmethod="normal-reference")
summary(bw1)
fitv1 <- fitted(npcdist(bws=bw1,newdata=z))


##step 1 calculate v

nx<-n
dfedu<-rep(0,nx)
ys<-z$y
xs1<-z$x1
xs2<-z$x2

v<-rep(0,nx)
for (i in 1:nx){
  v[i]<-condist(ys[i],y<-z$y,hy<-bw$ybw,xs1[i],x1<-z$x1,hx1<-bw$xbw[1],xs2[i],x2<-z$x2,lamdx2<-bw$xbw[2])
}



v1<-rep(0,nx)
for (i in 1:nx){
  v1[i]<-condist(ys[i],y<-z$y,hy<-bw1$ybw,xs1[i],x1<-z$x1,hx1<-bw1$xbw[1],xs2[i],x2<-z$x2,lamdx2<-bw1$xbw[2])
}


z1<-data.frame(z,v,fitv,v1,fitv1)

###step 2

condistv<-function(ys,y,hy,xs1,x1,hx1,xs2,x2,lamdx2,vs,v,hv){
  res<-function(ys,y,hy,xs1,x1,hx1,xs2,x2,lamdx2){
    sum(pnorm((ys-y)/hy)*ckf((x1-xs1)/hx1)/hx1*
          dkf(xs2,x2,lamdx2))*ckf((v-vs)/hv)/hv/
    (sum(ckf((x1-xs1)/hx1)/hx1*dkf(xs2,x2,lamdx2)*ckf((v-vs)/hv)/hv))}
  return(res)
}


g<-rep(0,nx)
for (i in 1:nx){
  g[i]<-condist(ys[i],y<-z$y,hy<-bw1$ybw,xs1[i],x1<-z$x1,hx1<-bw1$xbw[1],xs2[i],x2<-z$x2,lamdx2<-bw1$xbw[2])
}

###step 2 evaluative v kernel

kernelv<-function(vs,v,hv){
  res<-ckf((v-vs)/hv)/hv
  return(res)
  }


vs<-z1$v
weightv<-rep(0,nx)
for (i in 1:nx){
  weightv[i]<-kernelv(vs[i],v=z1$v,hv=0.02)
}
weightv


###step 2 
##find bandwith with new v


bw2 <- npcdistbw(formula=y~x1+ordered(x2)+v)
summary(bw2)


##step 3 sum for set v
condistaddv<-function(ys,y,hy,xs1,x1,hx1,xs2,x2,lamdx2,vs,v,hv){
  res<-
    sum(pnorm((ys-y)/hy)*ckf((x1-xs1)/hx1)/hx1*
          dkf(xs2,x2,lamdx2))*ckf((v-vs)/hv)/hv/
    (sum(ckf((x1-xs1)/hx1)/hx1*dkf(xs2,x2,lamdx2)*ckf((v-vs)/hv)/hv))
  return(res)
}



g<-function(ys,y,xs1,x1,xs2,x2,v){
res<-sum(condistaddv(ys,y,hy<-bw2$ybw,xs1,x1,hx1<-bw2$xbw[1],xs2,x2,lamdx2<-bw2$xbw[2],vs<-z1$v,v,hv<-bw2$xbw[3]))/n
return(res)
}

##step 4 find root 

#invcdf <- function(q){
#  uniroot(function(x){cdf(x) - q}, range(x))$root
#}

###
invgtreat<-function(q){
  uniroot(function(ys){g(ys,y<-z1$y,xs1<-(-1),x1<-z1$x1,xs2<-1,x2<-z1$x2,v<-z1$v)-q},range(z1$y))$root
}


invgcontrol<-function(q){
  uniroot(function(ys){g(ys,y<-z1$y,xs1<-(-1),x1<-z1$x1,xs2<-0,x2<-z1$x2,v<-z1$v)-q},range(z1$y))$root
}


qte<-invgtreat(0.5)-invgcontrol(0.5)

qte















