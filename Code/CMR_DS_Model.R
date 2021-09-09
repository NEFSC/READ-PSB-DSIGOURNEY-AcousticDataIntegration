#For New Sims


 #JAGS Code Availabiity Model
 


CMR_model<-function(nitt,thin,burnin,chains,dist_Below,dist_Above,W_Above,W_Below,
  ones.dist_A,y,z.mat,M,T,age,Bin.start,Bin.end,S.start,S.end,D.start,D.end){


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
   
for(j in 2:(T-1)){
   b0[j]~dnorm(0,0.0001)
   }
   
b0[1]<-b0[2]
b0[T]<-b0[T-1]

b.dist~dnorm(0,0.0001)
b.dist2~dnorm(0,0.0001)

   
b0.phi~dnorm(0,0.0001)
b.age~dnorm(0,0.0001) 
        
for(i in 1:M){
 for(j in 1:T){
   age[i,j]~dunif(0,j)
        }
         }           
          
for( i in 1:N_Above){

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
N.tot.Above<-N_Above*W_Above/mean.esw_Above  #Visual estimate only
   
       
for(i in 1:M){

dist_Below1[i]~dunif(0,W_Below) #Need distribution for unobserved distances
 
logit(p[i,1])<-b0[1]+b.dist*dist_Below1[i]+b.dist2*dist_Below1[i]^2
z[i,1]~dbern(gamma1)
mu[i,1]<-z[i,1]*p[i,1]
y[i,1]~dbern(mu[i,1])
recruitable[i,1]<-1
everRecruit[i]<-1-z[i,1]

y.rep[i,1]~dbern(mu[i,1])
Eval[i,1]<-mu[i,1]
resid.true[i,1]<-pow(y[i,1]-Eval[i,1],2)/(Eval[i,1]+0.5)
resid.pred[i,1]<-pow(y.rep[i,1]-Eval[i,1],2)/(Eval[i,1]+0.5)

for(j in 2:T){
logit(p[i,j])<-b0[j]+b.dist*dist_Below1[i]+b.dist2*dist_Below1[i]^2
logit(phi[i,j-1])<-b0.phi+b.age*age[i,j-1]
survived[i,j]<- (phi[i,j-1]*z[i,j-1])
recruitable[i,j]<- recruitable[i,(j-1)]*(1-z[i,j-1])
mu2[i,j]<-  survived[i,j] + gamma[j-1]*recruitable[i,j]
z[i,j]~dbern(mu2[i,j])
mu3[i,j]<-z[i,j]*p[i,j]
y[i,j]~dbern(mu3[i,j])

#For caluclating Bayein p-value
y.rep[i,j]~dbern(mu3[i,j])
Eval[i,j]<-mu3[i,j]
resid.true[i,j]<-pow(y[i,j]-Eval[i,j],2)/(Eval[i,j]+0.5)
resid.pred[i,j]<-pow(y.rep[i,j]-Eval[i,j],2)/(Eval[i,j]+0.5)

  }
}

for(i in 1:M){
recruit[i,1]<-z[i,1]
emigrate[i,1]<-0
for(j in 2:T){
recruit[i,j]<-z[i,j]*(1-z[i,j-1])
emigrate[i,j]<-z[i,j-1]*(1-z[i,j])
}
}


for(j in 1:T){

N[j]<-sum(z[1:M,j])
D[j]<-sum(recruit[1:M,j])
S[j]<-sum(emigrate[1:M,j])

}

for(i in 1:M){

Nind[i]<-sum(z[i,Bin.start:Bin.end])
Nalive[i]<-1-equals(Nind[i],0)

}


#Only count duplicates that occur in interval of overlap
totDive<-sum(D[D.start:D.end])
totSurface<-sum(S[S.start:S.end])

duplicates<-totDive+totSurface

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
N_Above=N_Above, 
C=10000,
W_Above=W_Above,
W_Below=W_Below,
ones.dist_A= ones.dist_A, 
M=M,
T=T,
z.mat=z.mat, 
Bin.start=Bin.start,
Bin.end=Bin.end,
D.start=D.start,
D.end=D.end,
S.start=S.start,
S.end=S.end,
age=age
  )

# Initial values
inits <- function(){list(b0 = runif(T, -1, 0), b.dist = runif(1, -1, 0),b.dist2 = runif(1, -1, 0), b.df.0_Above=1, 
             gamma1= runif(1, 0.1, 0.9),  gamma= runif(T-1, 0, 0.1), b0.phi = runif(1, 0, 1),
              b.age = runif(1, 0, 1), z=z.mat)}


varsToMonitor<-c("N.total","availability","N.tot.Below","N.tot.Above","b0", "b.dist", 
"b.df.0_Above","mean.esw_Above","z","duplicates","D","S","N","gamma","gamma1","totDive","totSurface",
"phi", "b0.phi", "b.age","fit.true","fit.pred")


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


 
     
CMR_test<- CMR_model(nitt=nitt,thin=thin,burnin=burnin,chains=chains,
dist_Above=dist_Above,W_Above=W_Above,W_Below=W_Below,ones.dist_A=ones.dist_A,y=y,z.mat=z.mat,
M=M,T=T,age=age,Bin.start=Bin.start,Bin.end=Bin.end,S.start=S.start,S.end=S.end,D.start=D.start,
D.end=D.end)
##

 #Goodness of Fit
    plot(CMR_test$fit.pred, CMR_test$fit.true)
    abline(a=0,b=1)
    mean(CMR_test$fit.pred> CMR_test$fit.true)

                       n.sims<-length(CMR_test$N.total[1,,1])
                       
                       #Abundance Estimate
                       N.total_CMR_DS<-median(CMR_test$N.total[1,,1])
                       sd_N.total_CMR_DS<-sd(CMR_test$N.total[1,,1])
                       cv_N.total_CMR_DS<- sd_N.total_CMR_DS/N.total_CMR_DS
                       
                       sort_N.total_CMR_DS<-sort(CMR_test$N.total[1,,1]) #For calculating 95% credible intervals 
                       upper_N.total_CMR_DS<-sort_N.total_CMR_DS[round(0.975*n.sims)]
                       lower_N.total_CMR_DS<-sort_N.total_CMR_DS[round(0.025*n.sims)]
                       
                       ##Coverage (1=credible intervals include true value, 0=credible intervals do not include true value)
                       coverage_N.total_CMR_DS<-as.numeric(N.true<upper_N.total_CMR_DS)*as.numeric( N.true>lower_N.total_CMR_DS)

                      
                      #Availability Estimate 
                       availability_CMR_DS<-median(CMR_test$availability[1,,1])
                       sd_availability_CMR_DS<-sd(CMR_test$availability[1,,1])
                       cv_availability_CMR_DS<- sd_availability_CMR_DS/availability_CMR_DS
                       
                       sort_availability_CMR_DS<-sort(CMR_test$availability[1,,1]) #For calculating 95% credible intervals 
                       upper_availability_CMR_DS<-sort_availability_CMR_DS[round(0.975*n.sims)]
                       lower_availability_CMR_DS<-sort_availability_CMR_DS[round(0.025*n.sims)]
                       
                       #Coverage (1=credible intervals include true value, 0=credible intervals do not include true value)
                       coverage_availability_CMR_DS<-as.numeric(Availability<upper_availability_CMR_DS)*as.numeric(Availability>lower_availability_CMR_DS)
                       
                       
                       #Acoustic Estimate (i.e.,total below)
                       N.tot.Below_CMR_DS<-median(CMR_test$N.tot.Below[1,,1])
                       sd_N.tot.Below_CMR_DS<-sd(CMR_test$N.tot.Below[1,,1])
                       cv_N.tot.Below_CMR_DS<- sd_N.tot.Below_CMR_DS/N.tot.Below_CMR_DS
                       
                       sort_N.tot.Below_CMR_DS<-sort(CMR_test$N.tot.Below[1,,1]) #For calculating 95% credible intervals 
                       upper_N.tot.Below_CMR_DS<-sort_N.tot.Below_CMR_DS[round(0.975*n.sims)]
                       lower_N.tot.Below_CMR_DS<-sort_N.tot.Below_CMR_DS[round(0.025*n.sims)]
                       
                       #Coverage (1=credible intervals include true value, 0=credible intervals do not include true value)
                       coverage_N.tot.Below_CMR_DS<-as.numeric(N.tot.Below<upper_N.tot.Below_CMR_DS)*as.numeric(N.tot.Below>lower_N.tot.Below_CMR_DS)
                       
                       
                       #Visual Estimate (i.e.,total above)
                       N.tot.Above_CMR_DS<-median(CMR_test$N.tot.Above[1,,1])
                       sd_N.tot.Above_CMR_DS<-sd(CMR_test$N.tot.Above[1,,1])
                       cv_N.tot.Above_CMR_DS<- sd_N.tot.Above_CMR_DS/N.tot.Above_CMR_DS
                       
                       sort_N.tot.Above_CMR_DS<-sort(CMR_test$N.tot.Above[1,,1]) #For calculating 95% credible intervals 
                       upper_N.tot.Above_CMR_DS<-sort_N.tot.Above_CMR_DS[round(0.975*n.sims)]
                       lower_N.tot.Above_CMR_DS<-sort_N.tot.Above_CMR_DS[round(0.025*n.sims)]
                       
                       #Coverage (1=credible intervals include true value, 0=credible intervals do not include true value)
                       coverage_N.tot.Above_CMR_DS<-as.numeric(N.tot.Above<upper_N.tot.Above_CMR_DS)*as.numeric(N.tot.Above>lower_N.tot.Above_CMR_DS)


                     