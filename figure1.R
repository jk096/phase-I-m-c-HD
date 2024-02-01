rm(list=ls())
library(Matrix)
library(mnormt)
library(mvtnorm)
library(DepthProc)
library(PEtests)
library(parallel)
library(foreach)
library(doParallel)


####0.####
n=100
vv=20
B=5000
SIM=1000
alpha=0.05
delt=c(0.15,0.3,0.5,1,1.5)
#eta=c(0.01,0.05,0.1,0.15,0.2)
#eta=c(0.01,0.05,0.2,0.5,0.8)
eta=c(0.01,0.05,0.1,0.5,0.8)
####1.####
#1.1
#function：independent
cov_f1=function(p){
    omiga=diag(1,p)
    return(omiga)
}
#function：correlated
cov_f2=function(p){
    omiga=matrix(rep(0,p^2),p,p)
    for (i in 1:p) {
        for (j in 1:p) {
            omiga[i,j]=0.6^abs(i-j)
        }
    }
    return(omiga)
}

#1.1.2
U_f=function(p,eta,delt){
    A=matrix(0,p,p)
    pp=ceiling(p*eta)
    for (i in 1:pp) {
        for (j in 1:pp) {
            A[i,j]=delt
            
        }
        
    }
    return(A)
}


matrix_sqrt=function(A){
    eigen_decomp <- eigen(A)
    
    D <- diag(eigen_decomp$values)
    V <- eigen_decomp$vectors
    D_sqrt <- diag(sqrt(eigen_decomp$values))
    
    A_sqrt <- V %*% D_sqrt %*% t(V)
    return(A_sqrt)
    
}

#1.1.3
#normal
generate_f=function(n,p,miu,SIG){
    SIG_sqrt=matrix_sqrt(SIG)
    dat=matrix(NA,n,p)
    for (u in 1:n) {
        z=rnorm(p,0,1)
        dat[u,]=SIG_sqrt%*%z+miu
    }
    return(dat)
}
#t
generate_f=function(n,p,miu,SIG){
    SIG_sqrt=matrix_sqrt(SIG)
    dat=matrix(NA,n,p)
    for (u in 1:n) {
        z=rt(p,df = 4)/sqrt(2)
        dat[u,]=SIG_sqrt%*%z+miu
    }
    return(dat)
}
#gamma
generate_f=function(n,p,miu,SIG){
    SIG_sqrt=matrix_sqrt(SIG)
    dat=matrix(NA,n,p)
    for (u in 1:n) {
        z=rgamma(p,shape = 4,scale = 0.5)-4*0.5
        dat[u,]=SIG_sqrt%*%z+miu
    }
    return(dat)
}


# multivariate phase I mean&covariance depth-based——Li(2014)
#function:calculate SQ(n1)
SQ_f=function(x,n1,n){
    #depth
    D=sapply(1:n, function(j) depth(x[j,], x[1:n1,],method = 'Euclidean'))
    Di=D[1:n1]
    Dj=D[(n1+1):n]
    #R
    R=sapply(1:length(Dj), function(j) sum(Di < Dj[j]) + 1/2 * sum(Di == Dj[j]))
    #Q
    Q=sum(R)
    #SQ(n1)
    SQ=(n1 * (n - n1) / 2 - Q) / sqrt(n1 * (n-n1) * (n+1) / 12)
    return(SQ)
}

####2.simulation  function####
icteststa_f=function(n,p,miu1,SIG1,vv){
    #IC data
    x=generate_f(n,p,miu1,SIG1)
    
    #HD
    #T_meanandcovc=rep(NA,n-vv)
    Tchi_meanandcovc=Tcau_meanandcovc=SQi=rep(0,n-vv)
    for (l in vv:(n-vv)) {
        x_1=x[1:l,];x_2=x[(l+1):n,]
        
        tm=meantest(x_1,x_2,method='cq',delta=NULL)$stat
        tc=covtest(x_1,x_2,method='lc',delta=NULL)$stat
        Tcau_meanandcovc[l]=tm^2+tc^2
        Tchi_meanandcovc[l]=simultest.pe.chisq(x_1,x_2)$stat
        
        
        SQi[l]=SQ_f(x,l,n)
    }
    Tchi_meanandcov=max(Tchi_meanandcovc)
    Tcau_meanandcov=max(Tcau_meanandcovc)
    
    T_mult=max(SQi)
    
    return(c(Tchi_meanandcov,Tcau_meanandcov,T_mult))
}


#function:detection power and change point estimation
simu_f=function(n1,n2,p,miu1,miu2,SIG1,SIG2,vv,CL,SIM,tau){
    n=n1+n2
    
    CLchi_meanandcov=CL[1]
    CLcau_meanandcov=CL[2]
    CL_mult=CL[3]
    
    Tchi_meanandcovS=Tcau_meanandcovS=T_multS=rep(NA,SIM)
    tauchi_meanandcovS=taucau_meanandcovS=tau_multS=rep(NA,SIM)
    for (s in 1:SIM) {
        #IC & OC data
        x1=generate_f(n1,p,miu1,SIG1)
        x2=generate_f(n2,p,miu2,SIG2)
        x=rbind(x1,x2)
        
        #HD
        Tchi_meanandcovc=Tcau_meanandcovc=SQi=rep(0,n-vv)
        for (l in vv:(n-vv)) {
            x_1=x[1:l,];x_2=x[(l+1):n,]
            tm=meantest(x_1,x_2,method='cq',delta=NULL)$stat
            tc=covtest(x_1,x_2,method='lc',delta=NULL)$stat
            Tcau_meanandcovc[l]=tm^2+tc^2
            Tchi_meanandcovc[l]=simultest.pe.chisq(x_1,x_2)$stat
            SQi[l]=SQ_f(x,l,n)
        }
        
        Tchi_meanandcovS[s]=max(Tchi_meanandcovc)
        tauchi_meanandcovS[s]=which.max(Tchi_meanandcovc)
        Tcau_meanandcovS[s]=max(Tcau_meanandcovc)
        taucau_meanandcovS[s]=which.max(Tcau_meanandcovc)
        
        T_multS[s]=max(SQi)
        tau_multS[s]=which.max(SQi)
        
    }
    #detection power
    dpower1=sum(Tchi_meanandcovS>=CLchi_meanandcov)/SIM   #HD——chisq
    dpower2=sum(Tcau_meanandcovS>=CLcau_meanandcov)/SIM
    dpower3=sum(T_multS>=CL_mult)/SIM               #multivariate
    
    
    ##Post-signal diagnostic
    # τ's estimation
    
    ocid1=which(Tchi_meanandcovS>=CLchi_meanandcov)   #HD——chisq
    ocid2=which(Tcau_meanandcovS>=CLcau_meanandcov)
    ocid3=which(T_multS>=CL_mult)               #multivariate
    
    
    #HD——chisq
    if (length(ocid1)==0) {
        psd1_1=Inf
        psd1_2=psd1_3=-Inf
    }else{
        
        tau_heat1=tauchi_meanandcovS[ocid1]
        
        tau_dif1=abs(tau_heat1-tau)
        #Post-signal diagnostic
        psd1_1=mean(tau_dif1)      #mean bias
        psd1_2=sum(tau_dif1<=3)/SIM#pro of detecting the change point within 3 steps
        psd1_3=sum(tau_dif1<=5)/SIM#pro of detecting the change point within 5 steps
    }
    
    if (length(ocid2)==0) {
        psd2_1=Inf
        psd2_2=psd2_3=-Inf
    }else{
        
        tau_heat2=taucau_meanandcovS[ocid2]
        
        tau_dif2=abs(tau_heat2-tau)
        #Post-signal diagnostic
        psd2_1=mean(tau_dif2)      #mean bias
        psd2_2=sum(tau_dif2<=3)/SIM#pro of detecting the change point within 3 steps
        psd2_3=sum(tau_dif2<=5)/SIM#pro of detecting the change point within 5 steps
    }
    #multivariate
    if (length(ocid3)==0) {
        psd3_1=Inf
        psd3_2=psd3_3=-Inf
    }else{
        tau_heat3=tau_multS[ocid3]                      
        
        tau_dif3=abs(tau_heat3-tau)
        #Post-signal diagnostic
        psd3_1=mean(tau_dif3)
        psd3_2=sum(tau_dif3<=3)/SIM
        psd3_3=sum(tau_dif3<=5)/SIM
        
    }
    
    return(c(dpower1,dpower2,dpower3,
             psd1_1,psd1_2,psd1_3,
             psd2_1,psd2_2,psd2_3,
             psd3_1,psd3_2,psd3_3))
    
}

tau=n1=50;n2=n-n1




####3 simulation-p50####
p=50

#3.1 independent（cov_f1）
#3.1.1  calculate control limit
miu1=rep(0,p);SIG1=cov_f1(p)
t0=Sys.time();t0
icteststa=foreach(i=1:B,.combine = 'cbind')%dopar%icteststa_f(n,p,miu1,SIG1,vv)
CLchi_meanandcov=sort(icteststa[1,])[(1-alpha)*B]
CLcau_meanandcov=sort(icteststa[2,])[(1-alpha)*B]
CL_mult=sort(icteststa[3,])[(1-alpha)*B]
CL=c(CLchi_meanandcov,CLcau_meanandcov,CL_mult);CL
t1=Sys.time();t=t1-t0;t

#3.1.2 detection power and Post-signal diagnostic
#
miu2=matrix(NA,length(eta)*length(delt),p)
SIG2=array(NA,dim = c(p,p,length(eta)*length(delt)))
for (i in 1:length(eta)) {
    etaa=ceiling(eta[i]*p)
    for (j in 1:length(delt)) {
        miu2[(i-1)*length(eta)+j,]=c(rep(delt[j],etaa),rep(0,p-etaa))
        SIG2[,,(i-1)*length(eta)+j]=SIG1+U_f(p,eta[i],delt[j])
    }
    
}
t0=Sys.time();t0
result1=foreach(i=1:nrow(miu2),.combine = 'cbind')%dopar%
    simu_f(n1,n2,p,miu1,miu2[i,],SIG1,SIG2[,,i],vv,CL,SIM,tau)
result1_1=cbind(result1,CL)
t1=Sys.time();t=t1-t0;t
t

