install.packages('pracma')
library(pracma)
install.packages('zoo')
library(zoo)
library(ggplot2)
##################################### Q1 - Q3 #######################################
alpha = .35
beta = .99
sigma = 2
delta = .025
abar = 1
rbar = 1/beta + delta - 1
kbar = (alpha*abar/(rbar))**(1/(1-alpha))

klow = 0.01*kbar
khigh = 1.99*kbar
ks = 101
ktod = data.frame(linspace(klow, khigh, ks))
#ktom = data.frame(t(ktom))

k_mat = data.frame(rep(ktod,ks))
#ktod_m = data.frame(t(ktom_m))
#ktom_mat = as.matrix(ktom_m) 
#ktod_mat = as.matrix(ktod_m)

a = c(1.1, 0.678)
pmat = matrix(c(0.977,1-0.926,1-0.977,0.926),nrow=2,ncol=2)

dis = 1
tol = 1e-06
value_mat = matrix(0,nrow=2,ncol=ks)
VF_guess = matrix(0,nrow=2,ncol=ks)
VF = matrix(0,nrow=2,ncol=ks)
PF = matrix(0,nrow=2,ncol=ks)

Sys.time()
while (dis > tol){
  VF = matrix(0,nrow=2,ncol=ks)
  for (i in 1:2){
    cons = a[i]*k_mat^alpha + (1-delta)*k_mat - t(k_mat)
    ret = cons^(1-sigma)/(1-sigma)
    ret[cons<0] = -Inf
    p = matrix(rep(pmat[i,], ks), ks, byrow = T)
    value_mat =  ret + beta*p%*%VF_guess
    VF[i,] = apply(value_mat,1, max) 
    
    PF[i,] = ktod[apply(value_mat, 1, which.max),]
  }
  dis = abs(max(VF - VF_guess))
  VF_guess = VF
}
Sys.time()
                            ########## plot value function
plot = data.frame(t(VF))
plot$k = ktod$linspace.klow..khigh..ks.
names(plot)[1:2] = c('high','low')
ggplot() + 
  geom_line(data = plot[,c(1,3)], aes(x=k, y=high), color='red') +
  ggtitle("Value Function") +
  xlab("k")+
  ylab("value function") + 
  geom_line(data = plot[,c(2,3)], aes(x=k, y=low), color='blue') 
     
                            ########## plot policy function
plot1 = data.frame(t(PF))
plot1$k = ktod$linspace.klow..khigh..ks.
names(plot1)[1:2] = c('high','low')
ggplot() + 
  geom_line(data = plot1[,c(1,3)], aes(x=k, y=high), color='black') +
  ggtitle("Policy Function") +
  xlab("k")+
  ylab("k'") + 
  geom_line(data = plot1[,c(2,3)], aes(x=k, y=low), color='brown') 

                              ########## plot saving
plot2 = plot1
plot2$highsave = plot2$high - (1-delta)*plot2$k
plot2$lowsave = plot2$low - (1-delta)*plot2$k
ggplot() + 
  geom_line(data = plot2[,3:4], aes(x=k, y=highsave), color='green') +
  ggtitle("Savings") +
  xlab("k")+
  ylab("saving") + 
  geom_line(data = plot2[,c(3,5)], aes(x=k, y=lowsave), color='orange') 

#################################### Q4 ######################################
pih = 74/97
pil = 1-pih
pmat = matrix(c(0.977,1-0.926,1-0.977,0.926),nrow=2,ncol=2)

########### creating the sequence of k for 1000 periods, I used the steady state k as the k_0
kseq = function(k_0,days){
  kregh = lm(plot1$high ~ plot1$k)
  summary(kregh)
  d1 = kregh$coefficients[1]
  d2 = kregh$coefficients[2]
  k = NA
  for (i in 1:days){
    k[i] = d1 + k_0*d2
    k_0 = k[i]
  }
  return(k)
}
ksq = kseq(kbar,1000)

########## drawing random {A_t} for 1000 periods, checking std of outputs and calibrating Ah
Ah = 1.0065
Al = (1-pih*Ah)/pil
a = c(Ah, Al)
A = NA
Aseq = function(A_0,days){
  A = A_0
  set.seed(12345)
  for (i in 2:days){
    if (A[i-1] == Ah){
      A[i] = sample(a,1,replace = FALSE,pmat[1,])
    }else{
      A[i] = sample(a,1,replace = FALSE,pmat[2,])
    } 
  }
  return(A)
}
Asq = Aseq(Ah,1000)
ysq = Asq*ksq^alpha
std(log(ysq))
            ###### Ah should be around 1.0065 in order to keep the std of output around 1.8%
