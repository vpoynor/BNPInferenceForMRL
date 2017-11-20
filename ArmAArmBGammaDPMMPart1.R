#----------------------------------------------------
#------- Analysis of Arm A and Arm B Data -----------
#------------------ Part 1 MCMC ---------------------
#----------------------------------------------------



###### Packages Required: ####
###### MCMC Pack ####
###### mvtnorm #####

#---------------------------------------------------------------#
#---------------------------------------------------------------#
#------------------ Predefined Functions -----------------------#
#---------------------------------------------------------------#
#---------------------------------------------------------------#


############################# 
####### Distinct L's ########
#############################
### This function takes in the vector (L)
### that indicates which component each survival 
### time (may be right censored) belongs, and produces 
### the vector of the distinct components (L.star).

distinct=function(L){
  N=length(L)
  L.star=numeric()
  L.star[1]=L[1]    
  for(r in 2:N){               
    k=length(L.star)
    q=1                            
    for(s in 1:k){
      if(L[r]!=L.star[s]) {q=q+1}}  
    if(q!=k){L.star[(k+1)]=L[r]}}   # if the rth component of L 
    #does not match any of the current values in L.star, 
    #then it is a new distinct component
  return(L.star)
}


###################################################### 
### Indicator of what distribution to draw the Z_l ###
######################################################
### This function takes in the lth component along with 
### the vector of distinct components, and outputs a 0 if l 
### is not a distinct component and a 1 is l is a distinct component.

z.draw=function(l,L.star){
  n.star=length(L.star)
  q=0
  j=1
  while(l != L.star[j] && j <= n.star){  
    q=j
    j=j+1}
  if(q == n.star){q=0}   #l is not a member of L.star
  else{q=1}  #l is a member of L.star
  return(q)
}


########################### 
### Number of repeats M ###
###########################
### This function takes in the vector (L) that indicates 
### which component each survival time (may be right censored) 
### belongs along with the total number of components (N).  
### The function provides the vector of frequencies of each 
### distinct component for which each data was assigned. 

reps=function(L,N){
  M=c(rep(0,N))
  n=length(L)
  for(l in 1:N){
    for(i in 1:n){
      if(l == L[i]){
        M[l]=M[l] + 1}}}
  
  return(M)} 


############################### 
###### Stick-Breaking #########
###############################
### This function takes in the number of components (N), 
### the current value of alpha, and the current vector of 
### frequencies of the distinct components (M).  
### The function produces the vector of probabilities associated 
### with the N components via the Stick-breaking process.  
### The function also provides the latent parameters drawn from the 
### Beta distribution (representing where the stick is being broken) 
### used in computing the vector of probabilities.

SB=function(N,alpha,M){
  p=V=numeric()
  
  V[1]=rbeta(1,(1 + M[1]),(alpha+sum(M[2:N])))
  p[1]=V[1]
  if(V[1]==1){V[1]=0.9999999999}
  
  for(i in 2:(N-1)){
    V[i]=rbeta(1, (M[i] +1), (sum(M[(i+1):N]) + alpha))
    if(V[i] <= 0){V[i] = 1e-308}
    if(V[i]==1){V[i]=0.9999999999}
    p[i]=V[i]*prod((1-V[1:(i-1)]))
  }
  
  p[N]=1-sum(p[1:(N-1)])
  if(p[N] <= 0){p[N] = 1e-308}
  return(c(p,V))
}


###########################################################
##### Obtaining probabilities (p.t) for sampling L's ######
###########################################################
### This function takes in the probability vector (p), the data vector, 
### the current shape parameter (gam), the current scale parameter (sig), 
### and the vector (d) indicating right censoring-1 and observed values-0.  
### The function provides the probability vector of the components now 
### incorporating the data.

probs.L.distn=function(p,y,theta,phi,d){
  
  N=length(p)
  p.1=numeric()
  
  for(i in 1:N){
    
    if(d==0){
      p.1[i]=p[i]*dgamma(y,exp(theta[i]),exp(phi[i]))}
    if(d==1){
      p.1[i]=p[i]*(1-pgamma(y,exp(theta[i]),exp(phi[i])))}
    
    if(p.1[i]==1){p.1[i]=0.99999999999}
    if(p.1[i]==0){p.1[i]=1e-308}}
  
  p.2=sum(p.1)
  p.t=p.1/p.2
  return(p.t)
}



####################################################################################################
####################################################################################################
########################################
###### @@MCMC Begins HERE!!!!!@@ #######
########################################
####################################################################################################
####################################################################################################

#########################
### Define parameters ###
#########################


setwd("/Users/valeriepoynor/Dropbox (CSU Fullerton)/Re__Submission")


ArmB = read.table("ArmB.txt", header = TRUE)
y=ArmB$time
n = length(y)
d=ArmB$cens
#I=1200001
#N=50



################### Arm A #########
#N = 50
#I = 1205000
#a.mu=matrix(c(2.5,-3),ncol=1)
#B.mu=matrix(c(.21,0,0,.21),ncol=2)
#a.Sig=4
#B.Sig=matrix(c(.21,0,0,.21),ncol=2)
#a.alpha=3
#b.alpha=1


################### Arm B #########
#N = 50
#I = 1205000
#a.mu=matrix(c(2.6,-2.9),ncol=1)
#B.mu=matrix(c(.21,0,0,.21),ncol=2)
#a.Sig=4
#B.Sig=matrix(c(.21,0,0,.21),ncol=2)
#a.alpha=3
#b.alpha=1



L=matrix(0,nrow=n,ncol=(I+1))
theta = phi = p = matrix(0,nrow=N,ncol=(I+1))
Sig=array(dim=c(2,2,(I+1)))
mu=matrix(0,nrow=2,ncol=(I+1))
n.star=alpha=numeric(length=(I+1))


###################
### Initialize ###
###################

tc=rep(a.mu[1],N)
pc=rep(a.mu[2],N)
mc=a.mu
Sc=B.mu
ac=1
Lc=numeric(n)
#Lc=c(rep(seq(1,30,1),2),31,32) #ArmA
#Lc=c(rep(seq(1,29,1),2),30) #ArmB
L.star=distinct(Lc)
nc=length(L.star)
accept=0

###################
### Loop Begins ###
###################

# Start the clock!
ptm <- proc.time()

for(i in 1:I){
  # Prints every 100th iteration
  if(i %% 10000==0){
    cat(paste0("iteration: ", i, "\n"))
  }
  
  
  ####################### 
  ####theta and phi ####
  #######################
  
  M=reps(Lc,N)
  
  for(l in 1:N){
    if(0==z.draw(l,L.star)){
      D=rmvnorm(1,mc, Sc)
      tc[l]=D[1]
      pc[l]=D[2]
    }else{
      
      S2=matrix(c(0.1,  0.1,  0.1,  0.5),ncol=2) 
      Pro=rmvnorm(1,c(tc[l],pc[l]),1*S2)
      Theta=c(tc[l],pc[l])
      
      k.c=k.o=0
      r.o=r.c=numeric()
      for(h in 1:n){
        if(Lc[h]==l){
          if(d[h]==0){k.o=k.o+1
          r.o[k.o]=y[h]}else{
            k.c=k.c+1
            r.c[k.c]=y[h]}
        }}
      
      if(runif(1) < exp(
        
        log(dmvnorm(Pro,mc,Sc)) 
        + sum(log(dgamma(r.o,exp(Pro[1]),exp(Pro[2])) )) 
        +  sum(log(pmax(1-(pgamma(r.c,exp(Pro[1]),exp(Pro[2]) )) , 1e-308)))
        - log(dmvnorm(Theta,mc,Sc))
        - sum(log(dgamma(r.o,exp(Theta[1]),exp(Theta[2]))))
        - sum(log(pmax(1-(pgamma(r.c,exp(Theta[1]),exp(Theta[2]))) , 1e-308)))
        
      ) ){
        
        tc[l]=Pro[1]
        pc[l]= Pro[2]
        accept=accept+1}else{
          tc[l]=Theta[1]
          pc[l]=Theta[2]}
      
    }}
  
  
  #############
  #### p's ####
  #############
  
  pV=SB(N,ac,M)
  prc=pV[1:N]
  
  ############# 
  #### L's ####
  #############
  
  for(k in 1:n){
    
    p.t=probs.L.distn(prc,y[k],tc,pc,d[k])
    
    Lc[k]=sample.int(N,1,replace=TRUE,prob=p.t)
  }
  
  ##########################################################
  ### Obtaining vector and number of distinct components ###
  ##########################################################
  
  L.star=distinct(Lc)
  
  nc=length(L.star)
  
  ###############################
  ### Prior Parameter Updates ###
  ###############################
  
  S.mu=solve(solve(B.mu) + (nc*solve(Sc)))
  m.mu=S.mu%*%(solve(B.mu)%*%a.mu + (solve(Sc)%*%matrix(c(sum(tc[L.star]),sum(pc[L.star])),ncol=1) ) )
  mc=t(rmvnorm(1,m.mu,S.mu))
  
  mm=matrix(c(tc[L.star],pc[L.star]),ncol=nc,nrow=2,byrow=TRUE) 
  
  df=nc + a.Sig
  #place1=matrix(0,nrow=2,ncol=nc) 
  #for(l in 1:nc){ place1[,l]=mm[,l] - mc}
  mat= B.Sig + (mm - c(mc))%*%t(mm-c(mc))
  Sc=riwish(df,mat)
  
  ##################
  ##### ALPHA ######
  ##################
  
  V=pV[(N+1):(2*N)]
  sh=N+a.alpha-1
  rt=-sum(log(1-V[1:(N - 1)])) + b.alpha
  ac=rgamma(1,sh,rate=rt)
  
  
  L[,(i+1)]=Lc
  theta[,(i+1)]=tc
  phi[,(i+1)]=pc
  mu[,(i+1)]=mc
  n.star[(i+1)]=nc
  Sig[,,(i+1)]=Sc
  alpha[(i+1)]=ac
  p[,(i+1)]=prc
}

# Stop the clock
proc.time() - ptm


##############################
######### barplot n.star ##########
##############################
# Plots the number of Distinct Components

r=numeric(N)
for(j in 1:N){
  for(i in 1:((I+1))){
    if(n.star[i]==j)
      r[j]=r[j]+1}
}
barplot(r, names.arg = seq(1,N,1))




#-----------------------------------------------------#
#------------------- THINNING ------------------------#
#-----------------------------------------------------#
#

EF=2000 #Effective Sample Size.  We ran the MCMC long enough to have
#2000 posterior independent samples
thetath=phith=pth=matrix(0,nrow=N,ncol = EF)
alphath = numeric(EF)
Lth= matrix(0, nrow = n, ncol =EF)
muth = matrix(0, nrow=2, ncol = EF)
n.starth = numeric(EF)
Sigth = array(dim=c(2,2,(EF)))
B=5000
th = 600
for(i in 1:EF){
  thetath[,i]= theta[,th*(i) +B]
  phith[,i]=phi[,th*(i)+B]
  pth[,i]= p[,th*(i)+B]
  alphath[i] = alpha[th*i +B]
  n.starth[i] =n.star[th*i +B]
  muth[,i] = mu[,i*th +B]
  Sigth[,,i] = Sig[,,th*i +B]
  Lth[,i] = L[,th*i +B]
}



