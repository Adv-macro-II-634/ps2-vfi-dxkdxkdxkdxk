install.packages('pracma')
library(pracma)
install.packages('zoo')
library(zoo)
library(ggplot2)
############################ Q1 - Q3
alpha = .35
beta = .99
sigma = 2
delta = .025
abar = 1
rbar = 1/beta + delta - 1
kbar = (alpha*abar/(rbar))**(1/(1-alpha))

klow = 0.01*kbar
khigh = 1.99*kbar
ks = 51
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

plot = data.frame(t(VF))
plot$k = ktod$linspace.klow..khigh..ks.
names(plot)[1:2] = c('high','low')
ggplot() + 
  geom_line(data = plot[,c(1,3)], aes(x=k, y=high), color='red') +
  ggtitle("Value Function") +
  xlab("k")+
  ylab("value function") + 
  geom_line(data = plot[,c(2,3)], aes(x=k, y=low), color='blue') 
  
plot = data.frame(t(PF))
plot$k = ktod$linspace.klow..khigh..ks.

names(plot)[1:2] = c('high','low')
ggplot() + 
  geom_line(data = plot[,c(1,3)], aes(x=k, y=high), color='red') +
  ggtitle("Value Function") +
  xlab("k")+
  ylab("value function") + 
  geom_line(data = plot[,c(2,3)], aes(x=k, y=low), color='blue') 

################ Q4

