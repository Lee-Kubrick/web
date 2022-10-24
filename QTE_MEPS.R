###########10/22/2022 Li Lu########################
library(np)
library(haven)
library(dplyr)
set.seed(123)


library(foreach)
library(doParallel)
library(iterators)

#####Input the data contains Y for the continous variable (expenditure)
#####X1 is the some endogenous (race, age, education, insurance) 
#####       and exogenous variable in the m(.) function,  y <- m(x) + error
#####Z is the instrumental variable (pain related condition) 
#####       enter in the selection equation and all of X1 will and x2(other controls) enter as well.
#####w is the treatment choice of opioid use 

mepsdata <- read.csv("E:/Paper/Professor Liang/Health application/Job market paper writting/Stata code for data into R_2/meps1819expandage.csv")

##Recode sex 1 female, 0 male
mepsdata$sex<-replace(mepsdata$sex,mepsdata$sex==1,0)
mepsdata$sex<-replace(mepsdata$sex,mepsdata$sex==2,1)


###Define variable for estimate v 



w<-mepsdata$pd

x1continous<-cbind(mepsdata$logincome)

x2discrete<-cbind(mepsdata$age,mepsdata$race,mepsdata$sex,mepsdata$region,mepsdata$insure,mepsdata$pain2,mepsdata$pain3,
              mepsdata$marriage,mepsdata$occupation,mepsdata$employ)

x3ordered<-cbind(mepsdata$educ,mepsdata$pain1,mepsdata$familysize)

##### Define the kernel function
kf<-function(v){
  apply(exp(-v^2/2)/sqrt(2*pi),1,prod)
}

###For discrete type kernel, such as Li and Racine l(xi; x; \lamda) = 1 if xi = x, and \lamda if xi != x.
dkf <-function(xi,x,lamd){
  apply(ifelse(xi!=x,lamd,1),1,prod)
}

###For discrete order type kernel, such as Li and Racine
###l(xi; x; \lamda) = 1 if xi = x, and \lamda^abs(xi-x) if abs(xi-x)>=1

dorderdf <-function(xi,x,lamd){
  apply(ifelse(abs(xi-x)>=1,lamd^(abs(xi-x)),1),1,prod)
}



#####

###use np() package calculate v 
##cross validation bw


vsample = subset(mepsdata, select = -c(DUPERSID,year,logexp,exp,employ))
attach(vsample)

bw <- npcdistbw(formula = factor(pd) ~ ordered(age) + logincome +
                 factor(race) + factor(sex) + factor(region) +  factor(insure) + 
                             factor(pain2)  + factor(pain3) + factor(marriage) + factor(occupation)  +
                  ordered(educ) + ordered(pain1) + ordered(familysize))
summary(bw)

fitv <- fitted(npcdist(bws=bw,newdata=vsample))

detach(vsample)
mepsdata$v <- fitv

## need add own function creat v


##





##Create G() function to find QTE
##Find the bandwith for estimate G() through np()


# gsample = subset(mepsdata, select = c(exp, race, age, educ, sex, region, insure, pd, v))
# attach(gsample)
# bw2 <- npcdistbw(formula = exp ~ factor(race) + age + ordered(educ) + factor(sex) + factor(region) + factor(insure)+
#                    factor(pd) + v)
# summary(bw2)
# detach(gsample)

#####Get variable for G()
n <- nrow(mepsdata)
yg<-mepsdata$exp

#x1cg<-cbind(mepsdata$age)

x2dg<-cbind(mepsdata$race,mepsdata$sex,mepsdata$region,mepsdata$insure)

x3og<-cbind(mepsdata$age,mepsdata$educ)

wg<-mepsdata$pd

vg<-mepsdata$v

##bandwidth
hy <- c(1.06*sd(yg)/n^(.2))
#hx1 <-c(1.06*sd(x1cg)/n^(.2))


lamdx2<-c(1.06*apply(x2dg,2,sd)/n^(.2))
lamdx3<-c(1.06*apply(x3og,2,sd)/n^(.2))


hv<-c(1.06*sd(vg)/n^(.2))
###Define the G() function to average out v
nw<-length(wg)

G<-function(ys,hy,x2dg,x2dgs,lamdx2,x3og,x3ogs,lamdx3,wgs,hv){
  out<-rep(0,nw)
  for (i in 1:nw){

    xr2<-nrow(x2dg)
    xr3<-nrow(x3og)

    u2<-kronecker(t(x2dgs),array(1,dim=c(xr2,1)))
    u3<-kronecker(t(x3ogs),array(1,dim=c(xr3,1)))
    u4<-kronecker(t(wgs),array(1,dim=c(length(wg),1)))
    vgi<-kronecker(vg[i],array(1,dim=c(length(vg),1)))
    hyn<-kronecker(hy,array(1,dim=c(length(yg),1)))

    lamdx2n<-kronecker(t(lamdx2),array(1,dim=c(xr2,1)))
    lamdx3n<-kronecker(t(lamdx3),array(1,dim=c(xr3,1)))
    hvn<-kronecker(t(hv),array(1,dim=c(length(vgi),1)))
    out<-sum(pnorm((ys-yg)/hyn)* dkf(x2dg,u2,lamdx2n)*
               dorderdf(x3og,u3,lamdx3n)*
               (wg==wgs)*
               kf((vg-vgi)/hvn))/sum(dkf(x2dg,u2,lamdx2n)*
                                       dorderdf(x3og,u3,lamdx3n)*
                                       (wg==wgs)*
                                       kf((vg-vgi)/hvn))
  }
  return(mean(out))
}

###set value
###treatment 
invgtreat<-function(ra,se,re,ins,ag,ed,q){
  uniroot(function(y1){G(y1,hy<-hy,
                         x2dg<-x2dg,x2dgs<-c(ra,se,re,ins),lamdx2<-lamdx2,
                         x3og<-x3og,x3ogs<-c(ag,ed),lamdx3<-lamdx3,
                         wgs<-1,hv<-hv)-q},c(-1000000,1000000))$root
}


##control
invgcontrol<-function(ra,se,re,ins,ag,ed,q){
  uniroot(function(y0){G(y0,hy<-hy,
                         x2dg<-x2dg,x2dgs<-c(ra,se,re,ins),lamdx2<-lamdx2,
                         x3og<-x3og,x3ogs<-c(ag,ed),lamdx3<-lamdx3,
                         wgs<-0,hv<-hv)-q},c(-1000000,1000000))$root
}


qte<-function(ra,se,re,ins,ag,ed,q){
  res<-
  invgtreat(ra,se,re,ins,ag,ed,q)-invgcontrol(ra,se,re,ins,ag,ed,q)
  return(res)
}




J <- seq(0.05,0.95,by=0.1)
nJ <- length(J)

#race and sex at (insure and region =1) at [age=1 25 44] at [educ=3]
QTE10<-rep(0,nJ)
QTE20<-rep(0,nJ)
QTE30<-rep(0,nJ)
QTE40<-rep(0,nJ)
QTE50<-rep(0,nJ)
QTE60<-rep(0,nJ)
QTE70<-rep(0,nJ)
QTE80<-rep(0,nJ)

#race and sex at (insure and region =2) at [age=1 25 44] at [educ=3]
QTE11<-rep(0,nJ)
QTE21<-rep(0,nJ)
QTE31<-rep(0,nJ)
QTE41<-rep(0,nJ)
QTE51<-rep(0,nJ)
QTE61<-rep(0,nJ)
QTE71<-rep(0,nJ)
QTE81<-rep(0,nJ)

#race and sex at (insure and region =3) at [age=1 25 44] at [educ=3]
QTE12<-rep(0,nJ)
QTE22<-rep(0,nJ)
QTE32<-rep(0,nJ)
QTE42<-rep(0,nJ)
QTE52<-rep(0,nJ)
QTE62<-rep(0,nJ)
QTE72<-rep(0,nJ)
QTE82<-rep(0,nJ)


#race and sex at (insure and region =4) at [age=1 25 44] at [educ=3]
QTE13<-rep(0,nJ)
QTE23<-rep(0,nJ)
QTE33<-rep(0,nJ)
QTE43<-rep(0,nJ)
QTE53<-rep(0,nJ)
QTE63<-rep(0,nJ)
QTE73<-rep(0,nJ)
QTE83<-rep(0,nJ)

#race and sex at (noinsure and region =1) at [age=1 25 44] at [educ=3]
QTE14<-rep(0,nJ)
QTE24<-rep(0,nJ)
QTE34<-rep(0,nJ)
QTE44<-rep(0,nJ)
QTE54<-rep(0,nJ)
QTE64<-rep(0,nJ)
QTE74<-rep(0,nJ)
QTE84<-rep(0,nJ)


#race and sex at (noinsure and region =2) at [age=1 25 44] at [educ=3]
QTE15<-rep(0,nJ)
QTE25<-rep(0,nJ)
QTE35<-rep(0,nJ)
QTE45<-rep(0,nJ)
QTE55<-rep(0,nJ)
QTE65<-rep(0,nJ)
QTE75<-rep(0,nJ)
QTE85<-rep(0,nJ)



#race and sex at (noinsure and region =3) at [age=1 25 44] at [educ=3]
QTE16<-rep(0,nJ)
QTE26<-rep(0,nJ)
QTE36<-rep(0,nJ)
QTE46<-rep(0,nJ)
QTE56<-rep(0,nJ)
QTE66<-rep(0,nJ)
QTE76<-rep(0,nJ)
QTE86<-rep(0,nJ)

#race and sex at (noinsure and region =4) at [age=1 25 44] at [educ=3]
QTE17<-rep(0,nJ)
QTE27<-rep(0,nJ)
QTE37<-rep(0,nJ)
QTE47<-rep(0,nJ)
QTE57<-rep(0,nJ)
QTE67<-rep(0,nJ)
QTE77<-rep(0,nJ)
QTE87<-rep(0,nJ)







# use multicore, set to the number of our cores, I set 8 here.

for (i in 1:nJ){
  QTE10[i]<-qte(1,0,1,1,1,3,J[i])
  QTE20[i]<-qte(2,0,1,1,1,3,J[i])
  QTE30[i]<-qte(3,0,1,1,1,3,J[i])
  QTE40[i]<-qte(4,0,1,1,1,3,J[i])
  QTE50[i]<-qte(1,1,1,1,1,3,J[i])
  QTE60[i]<-qte(2,1,1,1,1,3,J[i])
  QTE70[i]<-qte(3,1,1,1,1,3,J[i])
  QTE80[i]<-qte(4,1,1,1,1,3,J[i])
}  
  
plot(J,QTE10,xlab = "tau", ylab = "QTE at White and Male")
plot(J,QTE20,xlab = "tau", ylab = "QTE at Black and Male")
plot(J,QTE30,xlab = "tau", ylab = "QTE at Hispanic and Male")
plot(J,QTE40,xlab = "tau", ylab = "QTE at Other and Male")
plot(J,QTE50,xlab = "tau", ylab = "QTE at White and Female")
plot(J,QTE60,xlab = "tau", ylab = "QTE at Black and Female")
plot(J,QTE70,xlab = "tau", ylab = "QTE at Hispanic and Female")
plot(J,QTE80,xlab = "tau", ylab = "QTE at Other and Female")



for (i in 1:nJ){
  QTE11[i]<-qte(1,0,2,1,1,3,J[i])
  QTE21[i]<-qte(2,0,2,1,1,3,J[i])
  QTE31[i]<-qte(3,0,2,1,1,3,J[i])
  QTE41[i]<-qte(4,0,2,1,1,3,J[i])
  QTE51[i]<-qte(1,1,2,1,1,3,J[i])
  QTE61[i]<-qte(2,1,2,1,1,3,J[i])
  QTE71[i]<-qte(3,1,2,1,1,3,J[i])
  QTE81[i]<-qte(4,1,2,1,1,3,J[i])
}   
  



for (i in 1:nJ){
  QTE12[i]<-qte(1,0,3,1,1,3,J[i])
  QTE22[i]<-qte(2,0,3,1,1,3,J[i])
  QTE32[i]<-qte(3,0,3,1,1,3,J[i])
  QTE42[i]<-qte(4,0,3,1,1,3,J[i])
  QTE52[i]<-qte(1,1,3,1,1,3,J[i])
  QTE62[i]<-qte(2,1,3,1,1,3,J[i])
  QTE72[i]<-qte(3,1,3,1,1,3,J[i])
  QTE82[i]<-qte(4,1,3,1,1,3,J[i])
}
for (i in 1:nJ){
  QTE13[i]<-qte(1,0,4,1,1,3,J[i])
  QTE23[i]<-qte(2,0,4,1,1,3,J[i])
  QTE33[i]<-qte(3,0,4,1,1,3,J[i])
  QTE43[i]<-qte(4,0,4,1,1,3,J[i])
  QTE53[i]<-qte(1,1,4,1,1,3,J[i])
  QTE63[i]<-qte(2,1,4,1,1,3,J[i])
  QTE73[i]<-qte(3,1,4,1,1,3,J[i])
  QTE83[i]<-qte(4,1,4,1,1,3,J[i])
}
for (i in 1:nJ){
  
  
  QTE14[i]<-qte(1,0,1,0,1,3,J[i])
  QTE24[i]<-qte(2,0,1,0,1,3,J[i])
  QTE34[i]<-qte(3,0,1,0,1,3,J[i])
  QTE44[i]<-qte(4,0,1,0,1,3,J[i])
  QTE54[i]<-qte(1,1,1,0,1,3,J[i])
  QTE64[i]<-qte(2,1,1,0,1,3,J[i])
  QTE74[i]<-qte(3,1,1,0,1,3,J[i])
  QTE84[i]<-qte(4,1,1,0,1,3,J[i])
}
for (i in 1:nJ){
  
  QTE15[i]<-qte(1,0,2,0,1,3,J[i])
  QTE25[i]<-qte(2,0,2,0,1,3,J[i])
  QTE35[i]<-qte(3,0,2,0,1,3,J[i])
  QTE45[i]<-qte(4,0,2,0,1,3,J[i])
  QTE55[i]<-qte(1,1,2,0,1,3,J[i])
  QTE65[i]<-qte(2,1,2,0,1,3,J[i])
  QTE75[i]<-qte(3,1,2,0,1,3,J[i])
  QTE85[i]<-qte(4,1,2,0,1,3,J[i])
}
for (i in 1:nJ){
  
  QTE16[i]<-qte(1,0,3,0,1,3,J[i])
  QTE26[i]<-qte(2,0,3,0,1,3,J[i])
  QTE36[i]<-qte(3,0,3,0,1,3,J[i])
  QTE46[i]<-qte(4,0,3,0,1,3,J[i])
  QTE56[i]<-qte(1,1,3,0,1,3,J[i])
  QTE66[i]<-qte(2,1,3,0,1,3,J[i])
  QTE76[i]<-qte(3,1,3,0,1,3,J[i])
  QTE86[i]<-qte(4,1,3,0,1,3,J[i])
}
for (i in 1:nJ){
  QTE17[i]<-qte(1,0,4,0,1,3,J[i])
  QTE27[i]<-qte(2,0,4,0,1,3,J[i])
  QTE37[i]<-qte(3,0,4,0,1,3,J[i])
  QTE47[i]<-qte(4,0,4,0,1,3,J[i])
  QTE57[i]<-qte(1,1,4,0,1,3,J[i])
  QTE67[i]<-qte(2,1,4,0,1,3,J[i])
  QTE77[i]<-qte(3,1,4,0,1,3,J[i])
  QTE87[i]<-qte(4,1,4,0,1,3,J[i])
  
}






















































































































