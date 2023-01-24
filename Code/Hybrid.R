#New version for revised manuscript (for PeerJ) (Original verion on Github)

#Hybrid Model for Sim Testing


Hybrid_model<-function(nitt,thin,burnin,chains,W_Above,W_Below,Nobs_Above,Nobs_Below,dist_Above,dist_Below2,dist_Below3,
ones.dist_A,ones.dist_B,y,z.mat,M,T,Bin.start.H,Bin.end.H,start.H,end.H,Zone.H,ZoneS.H,ZoneF.H){


 sink("Hybrid.jags")
cat("
model {

#Priors for half-normal detection function

#Half-normal priors
b.df.0_Above ~ dnorm(0,0.0001)               
b.df.0_Below ~ dnorm(0,0.0001)

              
#Constants for intergating half-normal
pi <- 3.141593
a <- 0.147


#######################################
#Priors for cmr model
#######################################

gamma1~dunif(0,1)  

for(j in 1:(T-1)){
 gamma[j]~dunif(0,1)  
   }   

b0~dnorm(0,0.0001)
b.dist~dnorm(0,0.0001)
   
b0.phi~dnorm(0,0.0001)
b.time~dnorm(0,0.0001) 


for(i in 1:Nobs_Above){

##########################
#Half Normal
#######################

phi_A[i]<- L.f0_A[i]/C
ones.dist_A[i] ~ dbern(phi_A[i])
   
mu.df_A[i] <-exp(b.df.0_Above)
                 
sig.df_A[i] <- exp(mu.df_A[i])
sig2.df_A[i] <- sig.df_A[i]*sig.df_A[i]

# estimate effective strip width...
# ...(esw, which is also the integral from 0 to W of g(y))
# argument of 'error func' in the integrand for detect func g(y)

x_A[i] <-(W_Above/sqrt(2*sig2.df_A[i]))
  
# approximation of the 'error function'
erf_A[i] <- sqrt(1-exp(-x_A[i]*x_A[i]*(4/pi + a*x_A[i]*x_A[i])/(1 +
a*x_A[i]*x_A[i])))

#effective strip width 
esw_Above[i] <- sqrt(pi * sig2.df_A[i] / 2) * erf_A[i]
     

#Ones Trick
L.f0_A[i] <- exp(-dist_Above[i]*dist_Above[i] / (2*sig2.df_A[i])) * 1/esw_Above[i]   #half-normal likelihood

  }  #End loop
 
mean.esw_Above <- mean(esw_Above[])
N.tot.Above<-Nobs_Above*W_Above/mean.esw_Above  #Visual estimate only
 
for(i in 1:Nobs_Below){

##########################
#Half Normal
#######################

phi_B[i]<- L.f0_B[i]/C
ones.dist_B[i] ~ dbern(phi_B[i])
   
mu.df_B[i] <-exp(b.df.0_Below)
                 
sig.df_B[i] <- exp(mu.df_B[i])
sig2.df_B[i] <- sig.df_B[i]*sig.df_B[i]

# estimate effective strip width...
# ...(esw, which is also the integral from 0 to W of g(y))
# argument of 'error func' in the integrand for detect func g(y)

x_B[i] <-(W_Below/sqrt(2*sig2.df_B[i]))
  
# approximation of the 'error function'
erf_B[i] <- sqrt(1-exp(-x_B[i]*x_B[i]*(4/pi + a*x_B[i]*x_B[i])/(1 +
a*x_B[i]*x_B[i])))

#effective strip width 
esw_Below[i] <- sqrt(pi * sig2.df_B[i] / 2) * erf_B[i]
     

#Ones Trick
L.f0_B[i] <- exp(-dist_Below2[i]*dist_Below2[i] / (2*sig2.df_B[i])) * 1/esw_Below[i]   #half-normal likelihood

  }  #End loop
 
mean.esw_Below <- mean(esw_Below[])
N.tot.Below_all<-Nobs_Below*W_Below/mean.esw_Below  #Acoustic estimate only

#Adjust for the smaller zone of overlap with visual platform
N.tot.Below=N.tot.Below_all*W_Above/W_Below 
 
##############################################################          
#CMR Model
############################################################
 
 
for(i in 1:M){
dist_Below3[i]~dunif(0,W_Below) #Need distribution for unobserved distances
time1[i]~dunif(1,maxTime)

R[i,1]<-(YMid.vec[1]^2+dist_Below3[i]^2)^0.5 
logit(p[i,1])<-b0+b.dist*R[i,1]
 
z[i,1]~dbern(gamma1)
mu[i,1]<-z[i,1]*p[i,1]
y[i,1]~dbern(mu[i,1])
y.rep[i,1]~dbern(mu[i,1])
recruitable[i,1]<-1
everRecruit[i]<-1-z[i,1]
time[i,1]<-time1[i]*z[i,1]

Eval[i,1]<-mu[i,1]
resid.true[i,1]<-pow(y[i,1]-Eval[i,1],2)/(Eval[i,1]+0.5)
resid.pred[i,1]<-pow(y.rep[i,1]-Eval[i,1],2)/(Eval[i,1]+0.5)

for(j in 2:T){
R[i,j]<-(YMid.vec[j]^2+dist_Below3[i]^2)^0.5 
logit(p[i,j])<-b0+b.dist*R[i,j]

#logit(p[i,j])<-b0[j]+b.dist*dist_Below3[i]+b.dist2*dist_Below3[i]^2  
logit(phi[i,j-1])<-b0.phi+b.time*time[i,j-1]  
survived[i,j]<- (phi[i,j-1]*z[i,j-1]) #*(1-everRecruit[i])+ everRecruit[i]*z[i,j-1]
recruitable[i,j]<- recruitable[i,(j-1)]*(1-z[i,j-1])
mu2[i,j]<-  survived[i,j] + gamma[j-1]*recruitable[i,j]
z[i,j]~dbern(mu2[i,j])
mu3[i,j]<-z[i,j]*p[i,j]
y[i,j]~dbern(mu3[i,j])
time[i,j]<-(time[i,j-1]+1)*z[i,j]

#For caluclating Bayein p-value
y.rep[i,j]~dbern(mu3[i,j])
Eval[i,j]<-mu3[i,j]
resid.true[i,j]<-pow(y[i,j]-Eval[i,j],2)/(Eval[i,j]+0.5)
resid.pred[i,j]<-pow(y.rep[i,j]-Eval[i,j],2)/(Eval[i,j]+0.5)
  }
}

for(i in 1:M){
forage[i,1]<-z[i,1]
surface[i,1]<-0
for(j in 2:T){
forage[i,j]<-z[i,j]*(1-z[i,j-1])
surface[i,j]<-z[i,j-1]*(1-z[i,j])
}
}

#Derived parameters
for(j in 1:T){
N[j]<-sum(z[1:M,j])
F[j]<-sum(forage[1:M,j])
S[j]<-sum(surface[1:M,j])
}

for(i in 1:M){
Nind[i]<-sum(z[i,Bin.start.H:Bin.end.H])
Nalive[i]<-1-equals(Nind[i],0)
}

Nsuper<-sum(Nalive[1:M])

totSurface<-sum(S[start.H:end.H])
totForage<-sum(F[start.H:end.H])

#Average recruiting per interval
mean.F<-totForage/Zone.H
 
#Average surfacing per interval
mean.S<-totSurface/Zone.H
 
#Percent diving per interval
p.F<-mean.F/Nsuper

#Percent surfacing per interval
p.S<-mean.S/Nsuper

#Calculate duplicates
duplicates<-N.tot.Below*p.F*ZoneF.H+N.tot.Below*p.S*ZoneS.H

#Calculate abundance estimate adjusting for duplicates   
N.total<- N.tot.Below+N.tot.Above-duplicates

#Calculate surface availability (i.e., availability bias)    
availability<-N.tot.Above/N.total

#Final abundance estimate using DS-DS Method (ignores duplicates)
N.total2<- N.tot.Below+N.tot.Above

availability2<-N.tot.Above/N.total2
  
fit.true<-sum(resid.true)
fit.pred<-sum(resid.pred)
#
 } #end model


",fill = TRUE)
sink()

# Bundle data
Model.data <- list(y = y,
C=10000,
W_Above=W_Above,
W_Below=W_Below,
ones.dist_A=ones.dist_A,
ones.dist_B=ones.dist_B,
M=M,
T=T,
z.mat=z.mat,
Bin.start.H= Bin.start.H,
Bin.end.H=Bin.end.H,
start.H=start.H,
end.H=end.H,
Zone.H=Zone.H,
ZoneS.H=ZoneS.H,
ZoneF.H=ZoneF.H,
Nobs_Above=Nobs_Above,
Nobs_Below=Nobs_Below,
dist_Above=dist_Above,
dist_Below2=dist_Below2, #For DS analysis
dist_Below3=dist_Below3, #For CMR analysis of subset
maxTime=maxTime,
YMid.vec=YMid.vec+0.249
  )
  
 
# Initial values
inits <- function(){list(b0 = runif(1, 0, 2), b.dist = runif(1, 0, 1), b0.phi = runif(1, 8, 10),
                     b.time = runif(1, -1, 0),gamma= runif(T-1, 0, 0.1), z=z.mat, sigma.e=1, gamma1= runif(1, 0, 0.8),
                     ln_b.df.0_Below=log(1.5), ln_b_Below=log(2),ln_b.df.0_Above=log(1.5), ln_b_Above=log(2))}


varsToMonitor<-c("N.total","availability","N.tot.Below","N.tot.Above","b0", "b.dist",
"b.df.0_Above","mean.esw_Above","z","duplicates","F","S","N","gamma","gamma1","totForage","totSurface",
"phi", "b0.phi", "b.time","N.total2","availability2","mean.F","mean.Surface","p.F","p.S")


#MCMC settings
niter <- nitt
nthin <- thin
nburn <-burnin
nchains <- chains



(before <- Sys.time())
hybrid<- jags.model(
    file="./Hybrid.jags",
    data = Model.data,
    inits = inits,
    n.chains = nchains,
    n.adapt = nburn
)
(mid1 <- Sys.time())

(mid2 <- Sys.time())
out <- jags.samples(
    model = hybrid,
    variable.names = varsToMonitor,
    n.iter = niter,
    thin = nthin,
    progress.bar = 'text'
  )
(after <- Sys.time())

(diff(c(before,mid1,mid2,after)))[c(1,3)]

                out

 } #end

    