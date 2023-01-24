

  
sink("Hybrid_spermwhale.jags")
cat("
model {


#Priors for surface DS functions

#Hazard rate for VLT 
ln_b.df.0_Above ~ dnorm(0,0.0001)T(-2,2) 
ln_b_Above ~ dnorm(0,0.0001)T(-2,2)
b.df.0_Above<-exp(ln_b.df.0_Above)
b_Above<-exp(ln_b_Above)

#Hazard rate for PAM
ln_b.df.0_Below ~ dnorm(0,0.0001)T(-2,2) 
ln_b_Below ~ dnorm(0,0.0001)T(-2,2)              
b.df.0_Below<-exp(ln_b.df.0_Below)
b_Below<-exp(ln_b_Below)

#MR Component for VLT analysis
b1.0_Above ~ dnorm(0,0.0001)  #Groups 1 and 2
b2.0_Above ~ dnorm(0,0.0001)  #Sperm whale model (Group 3)
b.dist1_Above~ dnorm(0,0.0001)  #Sperm whale model (Group 3)
b.dist2_Above~ dnorm(0,0.0001)  #Sperm whale model (Group 3)
  
pi <- 3.141593
a <- 0.147

#Informed prior availability bias based on estimate from Palka et al. (2017)          
Av~dbeta(a_Av, b_Av)

for( i in 1:Nsights_Above){


#################################
 #DS Part of MRDS
#################################   
  
mu.df_Above[i] <- b.df.0_Above

esw_pred_Above[i,1]<-0

#Numerically integrate hazard rate function
for (c in 1:(steps-1)){
 dx1_Above[i,c]<- 1-exp(-((y.step_Above[c]/mu.df_Above[i])^(-b_Above)))
 dx2_Above[i,c]<- 1-exp(-((y.step_plus1_Above[c]/mu.df_Above[i])^(-b_Above)))
 esw_pred_Above[i,c+1]<-esw_pred_Above[i,c]+0.5*step.size_Above*(dx1_Above[i,c]+dx2_Above[i,c])
     }

esw_Above[i]<-esw_pred_Above[i,steps]

#Ones Trick to sample from likelihood
L.f0_Above[i] <- (1-exp(-((dist_Above[i]/ mu.df_Above[i])^(-b_Above)))) * 1/esw_Above[i]  # hazard rate likelihood

phi_Above[i]<- L.f0_Above[i]/C
ones.dist_Above[i] ~ dbern(phi_Above[i])
   
pcd_A[i]<-esw_Above[i]/W_Above

################################   
#MR part of MRDS
################################
    
logit(p1_Above[i])<-b1.0_Above+b.dist1_Above*dist_Above[i]
logit(p2_Above[i])<-b2.0_Above+b.dist2_Above*dist_Above[i]
     
logit(p0.1[i])<-b1.0_Above
logit(p0.2[i])<-b2.0_Above
p.dot_Above[i]<-p1_Above[i]+p2_Above[i]-p1_Above[i]*p2_Above[i]

pi_Above.1[i]<-p1_Above[i]*(1-p2_Above[i])/p.dot_Above[i] #y10
pi_Above.2[i]<-(1-p1_Above[i])*p2_Above[i]/p.dot_Above[i] #y01
pi_Above.3[i]<-p1_Above[i]*p2_Above[i]/p.dot_Above[i]    #y11

Obs_Above[i,1]~dbern(pi_Above.1[i])
Obs_Above[i,2]~dbern(pi_Above.2[i])
Obs_Above[i,3]~dbern(pi_Above.3[i])
     
g_0[i]<-p0.1[i]+p0.2[i]-p0.1[i]*p0.2[i]
       
Abundance_A[i]<-gs[i]/(pcd_A[i]*g_0[i]) 
   }

mean.esw_Above <- mean(esw_Above[])
   
#Calculate g(0)
P0_Above<-mean(g_0[])
 
N_Above<-sum(Abundance_A[])
N_Above_total<- N_Above/Av #Adjust for availability bias

######################
#DS model for PAM data
######################  
for( i in 1:Nsights_Below){

##########################
#Hazard Rate
#######################

mu.df_Below[i] <- b.df.0_Below
esw_pred_Below[i,1]<-0

#Numerically integrate hazard rate function
for (c in 1:(steps-1)){
  dx1_Below[i,c]<- 1-exp(-((y.step_Below[c]/mu.df_Below[i])^(-b_Below)))
  dx2_Below[i,c]<- 1-exp(-((y.step_plus1_Below[c]/mu.df_Below[i])^(-b_Below)))
  esw_pred_Below[i,c+1]<-esw_pred_Below[i,c]+0.5*step.size_Below*(dx1_Below[i,c]+dx2_Below[i,c])
     }

esw_Below[i]<-esw_pred_Below[i,steps]

#Ones Trick to sample from likelihood
L.f0_Below[i] <- (1-exp(-((dist_Below1[i]/ mu.df_Below[i])^(-b_Below)))) * 1/esw_Below[i]  # hazard rate likelihood

phi_Below[i]<- L.f0_Below[i]/C
ones.dist_B[i] ~ dbern(phi_Below[i])
   
pcd.B[i]<-esw_Below[i]/W_Below

Abundance_B[i]<-1/pcd.B[i]*1.12
}

mean.esw_Below <- mean(esw_Below[])

N_Below<-sum(Abundance_B[])
   
 
##############################################################          
#Fit CMR model to PAM data
############################################################
  
#gamma~dunif(0,1)    
gamma1~dunif(0,1) 

for(j in 1:(T-1)){
 gamma[j]~dunif(0,1)  
   }    

b0~dnorm(0,0.0001)    
b.dist~dnorm(0,0.0001)
   
b0.phi~dnorm(0,0.0001)
b.time~dnorm(0,0.0001) 


for(i in 1:M){

dist_Below[i]~dunif(0,W_Below.CH) #Need distribution for unobserved distances
time1[i]~dunif(1,maxTime)

R[i,1]<-(YMid.vec[1]^2+dist_Below[i]^2)^0.5 
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
R[i,j]<-(YMid.vec[j]^2+dist_Below[i]^2)^0.5 
logit(p[i,j])<-b0+b.dist*R[i,j]
logit(phi[i,j-1])<-b0.phi+b.time*time[i,j-1]  
survived[i,j]<- (phi[i,j-1]*z[i,j-1]) #*(1-everRecruit[i])+ everRecruit[i]*z[i,j-1]
recruitable[i,j]<- recruitable[i,(j-1)]*(1-z[i,j-1])
mu2[i,j]<- survived[i,j] + gamma[j-1]*recruitable[i,j]
z[i,j]~dbern(mu2[i,j])
mu3[i,j]<-z[i,j]*p[i,j]
y[i,j]~dbern(mu3[i,j])
time[i,j]<-(time[i,j-1]+1)*z[i,j]
 
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


for(j in 1:T){
N[j]<-sum(z[1:M,j])
F[j]<-sum(forage[1:M,j])
S[j]<-sum(surface[1:M,j])
}

for(i in 1:M){
Nind[i]<-sum(z[i,Bin.start.H:Bin.end.H])
Nalive[i]<-1-equals(Nind[i],0)
}

#Only count duplicates that occur in interval of overlap
totF<-sum(F[start.H:end.H])
totS<-sum(S[start.H:end.H])

Nsuper<-sum(Nalive[1:M])

#Average transitioning to foraging stage per bin
mean.F<-totF/Zone.H

#Average transitioning to surfacing stage per bin
mean.S<-totS/Zone.H

#Per capita number tranistioning to foraging per bin
p.F<- mean.F/Nsuper

#Percent surface per bin
p.S<- mean.S/Nsuper

#Calculate duplicates
duplicates<-N_Below*p.F*ZoneF.H+N_Below*p.S*ZoneS.H

#Total abundance adjusted for duplicates
N.total<- N_Below+N_Above-duplicates

#Availability bias
availability<-N_Above/N.total
  
fit.true<-sum(resid.true)
fit.pred<-sum(resid.pred)

 } #end model


",fill = TRUE)
sink()

#Bundle data
Model.data <- list(y = y,
  C=10000,
  W_Above=W_Above/1000, #Convert to KM
  W_Below=W_Below/1000,  #Convert to KM
  W_Below.CH=W_Below.CH/1000,  #Convert to KM
  ones.dist_Above=ones.dist_Above,
  ones.dist_B=ones.dist_B,
  steps=steps,
  y.step_Below=y.step_Below,
  y.step_plus1_Below=y.step_plus1_Below,
  step.size_Below=step.size_Below,
  y.step_Above=y.step_Above,
  y.step_plus1_Above=y.step_plus1_Above,
  step.size_Above=step.size_Above,
  M=M,
  T=T,
  start.H=start.H, #For averaging
  end.H=end.H, #For averaging
  Bin.start.H=Bin.start.H, #For Nsuper
  Bin.end.H=Bin.end.H, #For Nsuper
  a_Av=a_Av,
  b_Av=b_Av,
  Obs_Above=VLT_MRDS[,7:9],
  dist_Above=c(VLT_MRDS[,3])/1000+0.00001,
  gs=VLT_MRDS[,6],#Group size
  Nsights_Above=Nsights_Above,
  Nsights_Below=Nsights_Below,
  Zone.H=Zone.H,
  ZoneF.H=ZoneF.H,
  ZoneS.H=ZoneS.H,
  dist_Below=dist_Below/1000, #Convert to KM
  dist_Below1=dist_Below1/1000, #Convert to KM
  YMid.vec=y_mid.vec/1000, #Convert to KM
  maxTime=35
  )

#Initial values
inits <- function(){list(b0 = runif(1, 0, 2), b.dist = runif(1, 0, 1),b0.phi = runif(1, 0, 1),
                     b.time = runif(1, -1, 0), gamma= runif(T-1, 0.01, 0.07),  gamma1= runif(1, 0.1, 0.6), z=z.mat, 
                     ln_b.df.0_Below=log(1.5), ln_b_Below=log(2),ln_b.df.0_Above=log(1.5), ln_b_Above=log(2))}

 
varsToMonitor<-c("b0", "b.dist", "b.df.0_Above", "N_Above", "N.total", "mean.esw_Above", 
                  "mean.esw_Below", "b_Below", "availability", "duplicates","b_Above", 
                  "Nsuper","F", "S", "N", "gamma","gamma1","totF", "totS",  
                  "N_Above_total", "P0_Above", "b1.0_Above","b2.0_Above","b.dist1_Above",
                  "b.dist2_Above", "b0.phi", "b.time", "N_Below","b.df.0_Below","mean.F",
                  "mean.S","p.F","p.S", "fit.true", "fit.pred")

    # MCMC settings
niter <- nitt
nthin <- thin
nburn <-burnin
nchains <- chains



(before <- Sys.time())
hybrid<- jags.model(
    file="./Hybrid_spermwhale.jags",
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

             
  
        



 


 