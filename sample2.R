
# Get the cdf by numeric integration
#cdf <- function(x){
#  integrate(f, -Inf, x)$value
#}
# Use a root finding function to invert the cdf
#invcdf <- function(q){
#  uniroot(function(x){cdf(x) - q}, range(x))$root
#}




library(np)

set.seed(123)


n <- 200


x <- sort(c(rnorm(n,mean=-2,sd=0.5),rnorm(n,mean=3,sd=1.5)))

x.seq <- seq(-5,9,length=1000)


plot(x.seq,0.5*dnorm(x.seq,mean=-2,sd=0.5)+0.5*dnorm(x.seq,mean=3,sd=1.5),xlab="X",
     ylab="Mixture of Normal Densities",type="l",main="",col="blue",lty=1)



kf<-function(v){
  exp(-v^2/2)/sqrt(2*pi)
}

kfcon<-function(v){
  exp(-v^2/4)/sqrt(4*pi)
}

df1<-function(h,x1,y1,n1){
  u<-kronecker(y1,array(1,dim=c(1,n1)))
  v<-kronecker(t(x1),array(1,dim=c(n1,1)))
  A<-kfcon((u-v)/h)
  B<-kf((u-v)/h)
  out<-sum(A)/n1^2/h-2*(sum(B)-sum(diag(B)))/n1/(n1-1)/h
  return(out)
}

str(hopt <- optimize(df1, c(0.2*sd(x)/n^(.2), 3*sd(x)/n^(.2)), tol = 0.0001, x1 = array(x), y1 = array(x), n1 = length(x)))

x1<-seq(min(x),max(x),by=0.01)
nx<-length(x1)
n<-length(x)

df1<-function(x1,y1,n1,h){
  out<-sum(kf((y1-x1)/h))/n1/h
  return(out)
}

dfedu<-rep(0,nx)
h<-1.06*sd(x)/n^(.2)
for (i in 1:nx){
  dfedu[i]<-df1(x1[i],x,n,hopt$minimum)
}

dev.new()

plot(x1,dfedu, type='l', main="Kernel Density",
     xlab="x", ylab="Probability Density")


dfedu<-rep(0,nx)
h<-1.06*sd(x)/n^(.2)
for (i in 1:nx){
  dfedu[i]<-df1(x1[i],x,n,h)
}

dev.new()

plot(x1,dfedu, type='l', main="Kernel Density",
     xlab="x", ylab="Probability Density")












hist(x,prob=TRUE,main="")

plot(density(x),main="")

data("faithful",package="datasets")

fhat <- npudens(~waiting+eruptions,data=faithful)

plot(fhat,view="fixed",xtrim=-0.1,theta=310,phi=30,main="")






































