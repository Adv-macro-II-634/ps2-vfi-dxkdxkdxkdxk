install.packages('pracma')
library(pracma)
install.packages('zoo')
library(zoo)

alpha = .35
beta = .99
sigma = 2
delta = .025
abar = 1
rbar = 1/beta + delta - 1
kbar = (alpha*abar/(rbar))**(1/(1-alpha))

klow = 0.01*kbar
khigh = 1.99*kbar
ks = 100
kgrid = data.frame(linspace(klow, khigh, ks))
k_mat = data.frame(rep(kgrid,ks))

a = c(1.1, 0.678)
pmat = matrix(c(0.977,1-0.926,1-0.977,0.926),nrow=2,ncol=2)

dis = 1
tol = 1e-06
value_mat = matrix(0,nrow=2,ncol=ks)
VF_guess = matrix(0,nrow=2,ncol=ks)
VF = matrix(0,nrow=2,ncol=ks)
PF = matrix(0,nrow=2,ncol=ks)

Sys.time()
while (dis>tol){
  for (i in 1:2){
    value_mat = (a[i]*k_mat^alpha + (1-delta)*k_mat - k_mat)^(1-sigma)/(1-sigma) +
                    beta*pmat[i,]*VF_guess
    VF[i,] = sapply(value_mat, max) 
    PF[i,] = kgrid[apply(value_mat, 2, which.max),]
  }
}
Sys.time()
