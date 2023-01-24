#New version for revised manuscript (for PeerJ) (Original verion on Github)

#For New Sims


 #JAGS Code Availabiity Model
 

CMR_DS_model<-function(nitt,thin,burnin,chains,dist_Above,Nobs_Above,W_Above,dist_Below1,W_Below,
  ones.dist_A,y,z.mat,M,T,age,Bin.start,Bin.end,S.start,S.end,F.start,F.end,YMid.vec,
  maxTime){


 sink("CMR_DS.jags")
cat("
model {


#Priors for detection function

#Half-normal
b.df.0_Above ~ dnorm(0,0.0001)

#Constants for intergating half-normal
pi <- 3.141593
a <- 0.147

#Priors for cmr model
gamma1~dunif(0,1) 

for(j in 1:(T-1)){
     gamma[j]~dunif(0,1)  
   }  
   
b0~dnorm(0,0.0001)
b.dist~dnorm(0,0.0001)  
 
b0.phi~dnorm(0,0.0001)
b.time~dnorm(0,0.0001) 

          
for( i in 1:Nobs_Above){

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
          
for(i in 1:M){

dist_Below1[i]~dunif(0,W_Below) #Need distribution for unobserved distances
Time.Init[i]~dunif(1,maxTime)  #Need distribution for unobserved time 
                                #already in the foraging phase in the first bin

R[i,1]<-(YMid.vec[1]^2+dist_Below1[i]^2)^0.5 
logit(p[i,1])<-b0+b.dist*R[i,1]
z[i,1]~dbern(gamma1)
mu[i,1]<-z[i,1]*p[i,1]
y[i,1]~dbern(mu[i,1])
recruitable[i,1]<-1
everRecruit[i]<-1-z[i,1]
time[i,1]<-Time.Init[i]*z[i,1]

y.rep[i,1]~dbern(mu[i,1])
Eval[i,1]<-mu[i,1]
resid.true[i,1]<-pow(y[i,1]-Eval[i,1],2)/(Eval[i,1]+0.5)
resid.pred[i,1]<-pow(y.rep[i,1]-Eval[i,1],2)/(Eval[i,1]+0.5)

for(j in 2:T){
R[i,j]<-(YMid.vec[j]^2+dist_Below1[i]^2)^0.5 
logit(p[i,j])<-b0+b.dist*R[i,j]
logit(phi[i,j-1])<-b0.phi+b.time*time[i,j-1]
survived[i,j]<- (phi[i,j-1]*z[i,j-1])
recruitable[i,j]<- recruitable[i,(j-1)]*(1-z[i,j-1])
mu2[i,j]<-  survived[i,j] + gamma[j-1]*recruitable[i,j]
z[i,j]~dbern(mu2[i,j])
mu3[i,j]<-z[i,j]*p[i,j]
y[i,j]~dbern(mu3[i,j])
time[i,j]<-(time[i,j-1]+1)*z[i,j]

#For caluclating Bayesian p-value
y.rep[i,j]~dbern(mu3[i,j])
Eval[i,j]<-mu3[i,j]
resid.true[i,j]<-pow(y[i,j]-Eval[i,j],2)/(Eval[i,j]+0.5)
resid.pred[i,j]<-pow(y.rep[i,j]-Eval[i,j],2)/(Eval[i,j]+0.5)

  }
}

for(i in 1:M){
forage[i,1]<-z[i,1]*step(W_Above-dist_Below1[i])
surface[i,1]<-0
for(j in 2:T){
forage[i,j]<-z[i,j]*(1-z[i,j-1])*step(W_Above-dist_Below1[i])  #enter the foraging (i. e., clicking) phase
surface[i,j]<-z[i,j-1]*(1-z[i,j])*step(W_Above-dist_Below1[i])  #enter the surfacing phase
}
}

for(j in 1:T){

N[j]<-sum(z[1:M,j])
F[j]<-sum(forage[1:M,j])
S[j]<-sum(surface[1:M,j])

}

for(i in 1:M){

Nind[i]<-sum(z[i,Bin.start:Bin.end])*step(W_Above-dist_Below1[i])
Nalive[i]<-1-equals(Nind[i],0)

}


#Only count duplicates that occur in interval of overlap
totForage<-sum(F[F.start:F.end])
totSurface<-sum(S[S.start:S.end])

#Total duplicates in zone of overlap
duplicates<-totForage+totSurface

Nsuper<-sum(Nalive[1:M])

N.tot.Below<-Nsuper

N.total<-N.tot.Below+N.tot.Above-duplicates

availability<-N.tot.Above/N.total

fit.true<-sum(resid.true)
fit.pred<-sum(resid.pred)

  
 } #end model


",fill = TRUE)
sink()

#Bundle data
Model.data <- list(y = y,
dist_Below1=dist_Below1, 
dist_Above=dist_Above,
Nobs_Above=Nobs_Above, 
C=10000,
W_Above=W_Above,
W_Below=W_Below,
ones.dist_A= ones.dist_A, 
M=M,
T=T,
z.mat=z.mat, 
Bin.start=Bin.start,
Bin.end=Bin.end,
F.start=F.start,
F.end=F.end,
S.start=S.start,
S.end=S.end,
maxTime=maxTime,
#YMid.vec=YMid.vec+0.249
YMid.vec=y_bins-0.001
  )

# Initial values
inits <- function(){list(b0 =6, b.dist = -1.5, b.df.0_Above=1, 
            gamma1= runif(1, 0.6, 0.7), gamma= runif(T-1, 0, 0.1), b0.phi = runif(1, 6, 8),
              b.time = runif(1, -2, -0.5),  z=z.mat)}

  
varsToMonitor<-c("N.total","availability","N.tot.Below","N.tot.Above","b0", "b.dist",
"b.df.0_Above","mean.esw_Above","z","duplicates","F","S","N","gamma","gamma1","totForage","totSurface",
"phi", "b0.phi", "b.time","fit.true","fit.pred")


#MCMC settings
niter <- nitt
nthin <- thin
nburn <-burnin
nchains <- chains



(before <- Sys.time())
CMR<- jags.model(
    file="./CMR_DS.jags",
    data = Model.data,
    inits = inits,
    n.chains = nchains,
    n.adapt = nburn
)
(mid1 <- Sys.time())

(mid2 <- Sys.time())
out <- jags.samples(
    model = CMR,
    variable.names = varsToMonitor,
    n.iter = niter,
    thin = nthin,
    progress.bar = 'text'
  )
(after <- Sys.time())

(diff(c(before,mid1,mid2,after)))[c(1,3)]

out

 } #end
