#######QTE estimation (Li Lu; 10/16/2022)############

####We use np() package only for bandwidth estimation (cross validation method, one can use rule-of-thumb as well)

##1.call np() package and set seed at 123
library(np)
set.seed(123)

##2.Generate n and pseudo-sample data for analysis
#z table has three variable, y, x1 and x2 w.r.t continuous, continuous, discrete type.
n <- 200
z<- data.frame(y=c(rnorm(n,3,5)), x1=c(rnorm(n,mean=-2,sd=0.5)), x2=c(rbinom(n,1,0.75)))

##3.Define the kernel function for continuous/discrete 
# Continous use normal Guassian/ discrete use Li and Racine kernel, the smooth CDF kernel didn't define here, since
# one could use pnorm() function;
ckf<-function(v){
  exp(-v^2/2)/sqrt(2*pi)
}
dkf <-function(xi,x,lamd){
  res<-ifelse(xi!=x,lamd,1)
  return(res)
}

##4.Define Kernel Conditional Distribution Estimation with Mixed Data Types, we leave all bandwidth to be calculated 
##from np() package
#ys is the unknown value to be set, y is the value from original sample data z$y, hy is the bandwidth for y;
#xs1 is the unknown value to be set, x1 is the value from original sample data z$x1, hx1 is the bandwidth for x1;
#xs2 is the unknown value to be set, x2 is the value from original sample data z$x2, lamdx2 is the smooth parameter;
condist<-function(ys,y,hy,xs1,x1,hx1,xs2,x2,lamdx2){
  res<-
    sum(pnorm((ys-y)/hy)*ckf((x1-xs1)/hx1)/hx1*
          dkf(xs2,x2,lamdx2))/
    (sum(ckf((x1-xs1)/hx1)/hx1*dkf(xs2,x2,lamdx2)))
  return(res)
}


##5.Choose bandwidth and smooth parameter from np() package,such as function of npcdistbw()
#We do run fit value for each vi to double check whether our result by using our method of kernel CDF function 
#from step 4 is same with np() package.
y<-z$y
x1<-z$x1
x2<-z$x2

bw <- npcdistbw(formula=y~x1+ordered(x2))
summary(bw)
fitv <- fitted(npcdist(bws=bw,newdata=z))

##6.Calculate vi by using our kernel CDF from step 4, also using bandwidth from np() package
#nx is for loop n
#ys is the value to be set, also need to evaluate from this point, so for each loop i set ys<-z$y[i]
#xs1 and xs2 is same logic with ys
#y, x1, x2 are all the know sample data, bandwidth could be called with bw$xbw[j] command
#Finally we could get vi for each yi,x1i,x2i.
nx<-n
ys<-z$y
xs1<-z$x1
xs2<-z$x2

v<-rep(0,nx)
for (i in 1:nx){
  v[i]<-condist(ys[i],y<-z$y,hy<-bw$ybw,xs1[i],x1<-z$x1,hx1<-bw$xbw[1],xs2[i],x2<-z$x2,lamdx2<-bw$xbw[2])
}

####Rename y x1 x2 to the same dimension of ys, xs1, xs2 to double check


#without kronecker product, which means ys[i] (just 1) above and y (N) are not in the same dimension

v1<-condist(ys[1],y<-z$y,hy<-bw$ybw,xs1[1],x1<-z$x1,hx1<-bw$xbw[1],xs2[1],x2<-z$x2,lamdx2<-bw$xbw[2])



#with kronecker product, which means ys[i] (N) above and y (N) are in the same dimension
ysdb<-kronecker(y[1],array(1,dim=c(n,1)))

xs1db<-kronecker(x1[1],array(1,dim=c(n,1)))

xs2db<-kronecker(x2[1],array(1,dim=c(n,1)))


vdb_1<-condist(ysdb,y<-z$y,hy<-bw$ybw,xs1db,x1<-z$x1,hx1<-bw$xbw[1],xs2db,x2<-z$x2,lamdx2<-bw$xbw[2])

####


##7.Let's double check our result for vi and fitvi in step 5, which using np() kernels function
#And they are same, which means our kernel funntion defined in step 4 is safe.
z1<-data.frame(z,v,fitv)


##8.Let's do the model F_{y|x1,x2,v} kernel estimation
#8.1 First, find optimal bandwidth through np() package;
bw2 <- npcdistbw(formula=y~x1+ordered(x2)+v)
summary(bw2)

#8.2 Second, define a condistaddv() function for add variable v, similar function defined in step 4
condistaddv<-function(ys,y,hy,xs1,x1,hx1,xs2,x2,lamdx2,vs,v,hv){
  res<-
    sum(pnorm((ys-y)/hy)*ckf((x1-xs1)/hx1)/hx1*
          dkf(xs2,x2,lamdx2))*ckf((v-vs)/hv)/hv/
    (sum(ckf((x1-xs1)/hx1)/hx1*dkf(xs2,x2,lamdx2)*ckf((v-vs)/hv)/hv))
  return(res)
}
#8.3 Third, define a function g() which appled all bandwidth, and also set value vs to the data point from z1$v
g<-function(ys,y,xs1,x1,xs2,x2,v){
  res<-sum(condistaddv(ys,y,hy<-bw2$ybw,xs1,x1,hx1<-bw2$xbw[1],
                       xs2,x2,lamdx2<-bw2$xbw[2],vs<-z1$v,v,hv<-bw2$xbw[3]))/n
  return(res)
}

##9.Find the inverse root of function g() at specific point
#9.1 If we see x2 is binary treatment variable, then define set value xs2 at 1 for treated, as for xs1 I just
#set at -1, and other realized value y, x1, x2, v just by plugging in z1$y, z1$x1, z1$x2, z1$v
#As for find inverse root, we use uniroot function for ys, and with range of original interval to search root 
invgtreat<-function(q){
  uniroot(function(ys){g(ys,y<-z1$y,xs1<-(-1),x1<-z1$x1,xs2<-1,x2<-z1$x2,v<-z1$v)-q},range(z1$y))$root
}

#9.2 find control for xs2<-0, other xs1 is same with step 9.1
#q is the quantile we need to set as well, such as 0.5 is the median. 
invgcontrol<-function(q){
  uniroot(function(ys){g(ys,y<-z1$y,xs1<-(-1),x1<-z1$x1,xs2<-0,x2<-z1$x2,v<-z1$v)-q},range(z1$y))$root
}

##10 Find the QTE at \tau <- 0.2, which is q <-.2, at point xs1 <- (-1) point, since vs is already sum out.
qte<-function(q){
  res<-invgtreat(q)- invgcontrol(q)
  return(res)
}
qte(0.2)












































































































