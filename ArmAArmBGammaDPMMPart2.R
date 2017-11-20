#----------------------------------------------------
#------- Analysis of Arm A and Arm B Data -----------
#--------- Part 2 Posterior Inference ---------------
#----------------------------------------------------



###### Packages Required: ####
###### MCMC Pack ####
###### mvtnorm #####


#---------------------------------------------------------------#
#---------------------------------------------------------------#
#------------------ Predefined Functions -----------------------#
#---------------------------------------------------------------#
#---------------------------------------------------------------#
########################################
#### Mean Residual Life Function #######
########################################
### t is the survival time
### p is the vector of weights
### a is the vector of shape parameters
### a is the vector of rate parameters



mrl.fun = function(t,p,a,b){
  S = pgamma(t, a, b, lower.tail = FALSE)
  
  #part1 = exp(log(p) + a*log(t) -(t*b) + (a-1)*log(b))
  #part2 = gamma(a) 
  part1 = p*dgamma(t, (a+1), rate = b)
  part2 = a/(b^2)
  part3 = (p*S*a/b) - (p*S*t)
  temp1 = pmin(part1, 1e308) 
  temp2 = pmin(part2, 1e308)
  temp = sum((temp1*temp2 + part3))
  return(temp)
}




########################################
######## Survival Function #############
########################################
### t is the survival time
### p is the vector of weights
### a is the vector of shape parameters
### a is the vector of rate parameters

surv.fun = function(t,p,a,b){
  temp = sum(p*(pgamma(t, a,b, lower.tail = FALSE)))
  return(temp)
}

########################################
######### Density Function #############
########################################
### t is the survival time
### p is the vector of weights
### a is the vector of shape parameters
### a is the vector of rate parameters

den.fun = function(t,p,a,b){
  temp = sum(p*(dgamma(t, a,b)))
  return(temp)
}


#---------------------------------------------------------------#
################################################################
######## Posterior Predictive with Interval Estimates ##########
################################################################
#---------------------------------------------------------------#
y.grid=c(seq(0.000001,2000,2))
k =length(y.grid)
r.1A = r.2A = r.3A =r.4A=r.1B = r.2B = r.3B =r.4B= matrix(0, nrow = EF, ncol = k)
F_A=f_A=S_A=H_A=M_A=F_B=f_B=S_B=H_B=M_B=matrix(0,nrow=3,ncol=k)
QuantA = QuantB= matrix(0, nrow = 3, ncol = EF)
M.2A=M.8A=M.05A=M.95A=M.1A=M.9A=numeric(k)
M.2B=M.8B=M.05B=M.95B=M.1B=M.9B=numeric(k)
for(i in 1:EF){
  
  r.1A[i,] = apply(t(y.grid),2, den.fun, p=pthA[,(i)], a = exp(thetathA[,(i)]), 
                  b=exp(phithA[,(i)] ) )  
  r.1B[i,] = apply(t(y.grid),2, den.fun, p=pthB[,(i)], a = exp(thetathB[,(i)]), 
                   b=exp(phithB[,(i)] ) )  
  r.2A[i,] = apply(t(y.grid),2, surv.fun, p=pthA[,(i)], a = exp(thetathA[,(i)]), 
                  b=exp(phithA[,(i)] ) ) 
  r.2B[i,] = apply(t(y.grid),2, surv.fun, p=pthB[,(i)], a = exp(thetathB[,(i)]), 
                  b=exp(phithB[,(i)] ) ) 
  r.3A[i,] = apply(t(y.grid),2, mrl.fun, p=pthA[,(i)], a = exp(thetathA[,(i)]), 
                  b=exp(phithA[,(i)] ) )/r.2A[i,]
  r.3B[i,] = apply(t(y.grid),2, mrl.fun, p=pthB[,(i)], a = exp(thetathB[,(i)]), 
                  b=exp(phithB[,(i)] ) )/r.2B[i,]
  r.4A[i,] =   r.1A[i,]/r.2A[i,]
  r.4B[i,] =   r.1B[i,]/r.2B[i,]
  
  #----- Obtaining the Posterior Densities of the Quantiles     
  #QuantA[1,i] = y.grid[max(which(1-r.2A[i,]  <=0.25))]
  #QuantA[2,i] = y.grid[max(which(1-r.2A[i,]  <=0.5))]
  #QuantA[3,i] = y.grid[max(which(1-r.2A[i,]  <=0.9))]
  #QuantB[1,i] = y.grid[max(which(1-r.2B[i,]  <=0.25))]
  #QuantB[2,i] = y.grid[max(which(1-r.2B[i,]  <=0.5))]
  #QuantB[3,i] = y.grid[max(which(1-r.2B[i,]  <=0.9))]
  
  if(i %% 100==0){
    cat(paste0("iteration: ", i, "\n"))}}
{
f_A[2,]=apply(r.1A, 2, mean)
f_A[1,]=apply(r.1A, 2, quantile, probs = 0.025)
f_A[3,]=apply(r.1A, 2, quantile, probs = 0.975)

f_B[2,]=apply(r.1B, 2, mean)
f_B[1,]=apply(r.1B, 2, quantile, probs = 0.025)
f_B[3,]=apply(r.1B, 2, quantile, probs = 0.975)


S_A[2,]=apply(r.2A, 2, mean)
S_A[1,]=apply(r.2A, 2, quantile, probs = 0.025)
S_A[3,]=apply(r.2A, 2, quantile, probs = 0.975)


S_B[2,]=apply(r.2B, 2, mean)
S_B[1,]=apply(r.2B, 2, quantile, probs = 0.025)
S_B[3,]=apply(r.2B, 2, quantile, probs = 0.975)

H_A[2,]=apply(r.4A, 2, quantile, probs = 0.5)
H_A[1,]=apply(r.4A, 2, quantile, probs = 0.025)
H_A[3,]=apply(r.4A, 2, quantile, probs = 0.975)

H_B[2,]=apply(r.4B, 2, quantile, probs = 0.5)
H_B[1,]=apply(r.4B, 2, quantile, probs = 0.025)
H_B[3,]=apply(r.4B, 2, quantile, probs = 0.975)

F_A[2,]=1-S_A[2,]
F_A[1,]=1-S_A[1,]
F_A[3,]=1-S_A[3,]

F_B[2,]=1-S_B[2,]
F_B[1,]=1-S_B[1,]
F_B[3,]=1-S_B[3,]

M.1A=apply(r.3A, 2, quantile, probs = 0.1)
M.9A=apply(r.3A, 2, quantile, probs = 0.9)
M.2A=apply(r.3A, 2, quantile, probs = 0.2)
M.8A=apply(r.3A, 2, quantile, probs = 0.8)
M.05A=apply(r.3A, 2, quantile, probs = 0.05)
M.95A=apply(r.3A, 2, quantile, probs = 0.95)
M_A[2,]=apply(r.3A, 2, quantile, probs = 0.5)
M_A[1,]=apply(r.3A, 2, quantile, probs = 0.025)
M_A[3,]=apply(r.3A, 2, quantile, probs = 0.975)

M.1B=apply(r.3B, 2, quantile, probs = 0.1)
M.9B=apply(r.3B, 2, quantile, probs = 0.9)
M.2B=apply(r.3B, 2, quantile, probs = 0.2)
M.8B=apply(r.3B, 2, quantile, probs = 0.8)
M.05B=apply(r.3B, 2, quantile, probs = 0.05)
M.95B=apply(r.3B, 2, quantile, probs = 0.95)
M_B[2,]=apply(r.3B, 2, quantile, probs = 0.5)
M_B[1,]=apply(r.3B, 2, quantile, probs = 0.025)
M_B[3,]=apply(r.3B, 2, quantile, probs = 0.975)
}
Diff_AB = r.3A - r.3B


count.p = 0
prob.p=numeric(k)
for(j in 1:k){
  count.p = 0
  for(i in 1:EF){
    
    if(Diff_AB[i,j] > 0){
      count.p =count.p+1}
  } 
  
  prob.p[j] =  count.p/EF}


plot(y.grid, M_A[2,],xlim = c(0,2000),ylim = c(0,1500),type="l", lty = 2, lwd=2,
       xlab = "Time (days)", ylab="", main = "Mean Residual Life",cex.lab=1,cex.axis=1)
lines(y.grid, M_B[2,],lwd=2)
#lines(y.grid, M_A[1,])
#lines(y.grid, M_A[3,])

plot(y.grid, M_B[2,],xlim = c(0,2000),ylim = c(0,3500),type="l",
     main ="", xlab = "Time (days)", ylab = "",cex.lab=1,cex.axis=1)
lines(y.grid, M_B[1,])
lines(y.grid, M_B[3,])



plot(y.grid, S_A[1,], type ="l")
lines(y.grid, S_A[2,])
lines(y.grid, S_A[3,])

plot(y.grid, S_B[1,], type ="l")
lines(y.grid, S_B[2,])
lines(y.grid, S_B[3,])


plot(y.grid, f_A[1,], type ="l",ylim=c(0,.5))
lines(y.grid, f_A[2,])
lines(y.grid, f_A[3,])


plot(y.grid, f_B[1,], type ="l",ylim=c(0,.5))
lines(y.grid, f_B[2,])
lines(y.grid, f_B[3,])



#---------------------------------------------------------------#
################################################################
####################### Prior Predictive #######################
################################################################
#---------------------------------------------------------------#


N = 50
I1 = 2000
a.muA=matrix(c(2.5,-3),ncol=1)
a.muB=matrix(c(2.6,-2.9),ncol=1)
B.mu=matrix(c(.21,0,0,.21),ncol=2)
a.Sig=4
B.Sig=matrix(c(.21,0,0,.21),ncol=2)
a.alpha=3
b.alpha=1

y.grid=c(seq(0.000001,2000,2))
k =length(y.grid)

r.1pA = r.2pA = r.3pA =r.4pA=r.1pB = r.2pB = r.3pB =r.4pB= matrix(0, nrow = I1, ncol = k)


for(i in 1:I1){
  alphapA = rgamma(1, a.alpha, rate = b.alpha)
  alphapB = rgamma(1, a.alpha, rate = b.alpha)
  mupA = rmvnorm(1, a.muA,B.mu)
  SigpA= riwish(a.Sig,B.Sig)
  mupB = rmvnorm(1, a.muB,B.mu)
  SigpB= riwish(a.Sig,B.Sig)
  thetaA = rmvnorm(N, mupA, SigpA)
  thetaB = rmvnorm(N, mupB, SigpB)
  
  pVA=SB(N,alphapA, rep(0,N))
  prcA=pVA[1:N]
  
  pVB=SB(N,alphapB, rep(0,N))
  prcB=pVB[1:N]
  
  
  r.1pA[i,] = apply(t(y.grid),2, den.fun, p=prcA, a =  exp(thetaA[,1]), 
                   b=exp(thetaA[,2]) )  
  r.1pB[i,] = apply(t(y.grid),2, den.fun, p=prcB, a =  exp(thetaB[,1]), 
                   b=exp(thetaB[,2])  )  
  r.2pA[i,] = apply(t(y.grid),2, surv.fun, p=prcA, a = exp(thetaA[,1]), 
                   b=exp(thetaA[,2])  ) 
  r.2pB[i,] = apply(t(y.grid),2, surv.fun, p=prcB, a = exp(thetaB[,1]), 
                   b=exp(thetaB[,2])  ) 
  r.3pA[i,] = apply(t(y.grid),2, mrl.fun, p=prcA, a =  exp(thetaA[,1]), 
                   b=exp(thetaA[,2])  )/r.2A[i,]
  r.3pB[i,] = apply(t(y.grid),2, mrl.fun, p=prcB, a =  exp(thetaB[,1]), 
                   b=exp(thetaB[,2])  )/r.2B[i,]
  r.4pA[i,] =   r.1pA[i,]/r.2pA[i,]
  r.4pB[i,] =   r.1pB[i,]/r.2pB[i,]
  
 
  if(i %% 100==0){
    cat(paste0("iteration: ", i, "\n"))}}

Diffp_AB = r.3pA - r.3pB
count.prior = 0
prob.prior=numeric(k)
for(j in 1:k){
  count.prior = 0
  for(i in 1:I1){
    
    if(Diffp_AB[i,j] > 0){
      count.prior =count.prior+1}
  } 
  
  prob.prior[j] =  count.prior/I1}


#---------------------------------------------------------------#
################################################################
####################### Plot Differences #######################
################################################################
#---------------------------------------------------------------#

plot(y.grid, prob.p, xlim = c(0,2000),ylim = c(0.4,1),type="l", lty = 1, lwd=2,
     xlab = "Time (days)", ylab="Probability", main = "",cex.lab=1,cex.axis=1)

lines(y.grid, prob.prior,lty=2,lwd =2)

#time 0 is index 1
#time 100 is index 51
#time 250 is index 126
#time 500 is index 251
#time 800 is index 401
#time 1600 is index 801

plot(density(Diff_AB[,1]),xlim = c(-3000,3000),ylim =c(0,.002),
        xlab = "Difference in MRL time (days)", ylab="Density", main="Difference at time 0",lwd =2 )
abline(v = 0, lwd =2, lty =2)

plot(density(Diff_AB[,51]),xlim = c(-3000,3000),ylim =c(0,.002),
     xlab = "Difference in MRL time (days)", ylab="Density", main="Difference at time 100",lwd =2 )
abline(v = 0, lwd =2, lty =2)



plot(density(Diff_AB[,126]),xlim = c(-3000,3000),ylim =c(0,.002),
     xlab = "Difference in MRL time (days)", ylab="Density", main="Difference at time 250",lwd =2 )
abline(v = 0, lwd =2, lty =2)



plot(density(Diff_AB[,251]),xlim = c(-3000,3000),ylim =c(0,.002),
     xlab = "Difference in Time (days)", ylab="Density", main="Difference at time 500",lwd =2 )
abline(v = 0, lwd =2, lty =2)



plot(density(Diff_AB[,401]),xlim = c(-3000,3000),ylim =c(0,.002),
     xlab = "Difference in MRL time (days)", ylab="Density", main="Difference at time 800",lwd =2 )
abline(v = 0, lwd =2, lty =2)



plot(density(Diff_AB[,801]),xlim = c(-3000,3000),ylim =c(0,.002),
     xlab = "Difference in MRL time (days)", ylab="Density", main="Difference at time 1600",lwd =2 )
abline(v = 0, lwd =2, lty =2)

