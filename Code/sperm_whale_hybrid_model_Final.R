 

 sink("SA.jags")
cat("
model {

#Priors for surface detection function

#Hazard rate

#Visual
ln_b.df.0_Above ~ dnorm(0,0.0001)T(-2,2) 
ln_b_Above ~ dnorm(0,0.0001)T(-2,2)
               
#Acoustic
ln_b.df.0_Below ~ dnorm(0,0.0001)T(-2,2) 
ln_b_Below ~ dnorm(0,0.0001)T(-2,2)
b.df.0_Above<-exp(ln_b.df.0_Above)
b_Above<-exp(ln_b_Above)
b.df.0_Below<-exp(ln_b.df.0_Below)
b_Below<-exp(ln_b_Below)

       
#MR Component
b1.0_Above ~ dnorm(0,0.0001)  #Groups 1 and 2
b2.0_Above ~ dnorm(0,0.0001)  #Sperm whale model (Group 3)
b.dist1_Above~ dnorm(0,0.0001)  #Sperm whale model (Group 3)
b.dist2_Above~ dnorm(0,0.0001)  #Sperm whale model (Group 3)

#Priors for CMR analysis  
gamma~dunif(0,1)    
gamma1~dunif(0,1)  

for(j in 2:(T-1)){
 b0[j]~dnorm(0,0.0001)
   }
   
b0[1]<-b0[2]
b0[T]<-b0[T-1]

b.dist~dnorm(0,0.0001)
b.dist2~dnorm(0,0.0001)
   
b0.phi~dnorm(0,0.0001)
b.age~dnorm(0,0.0001) 
b.phi.dist~dnorm(0,0.0001) 
b.phi.dist2~dnorm(0,0.0001) 

#Priors for unknown ages
for(i in 1:M){
 for(j in 1:T){
  age[i,j]~dunif(0,j)
        }
         }
         
#Priors for unknown perpendicular distances         
for(i in 1:M){
   dist_Below[i]~dunif(0,W_Below.CH) #Need distribution for unobserved distances
   }
   
   
#Informative Prior for surface availability          
Av~dbeta(a_Av, b_Av)

###############################
#Loop through visual observation
#################################

for( i in 1:Nsights_Above){

#DS Part of MRDS

mu.df_Above[i] <- b.df.0_Above
esw_pred_Above[i,1]<-0

#Numerically integrate hazard rate function
for (c in 1:(steps-1)){
 dx1_Above[i,c]<- 1-exp(-((y.step_Above[c]/mu.df_Above[i])^(-b_Above)))
 dx2_Above[i,c]<- 1-exp(-((y.step_plus1_Above[c]/mu.df_Above[i])^(-b_Above)))
 esw_pred_Above[i,c+1]<-esw_pred_Above[i,c]+0.5*step.size_Above*(dx1_Above[i,c]+dx2_Above[i,c])
     }

esw_Above[i]<-esw_pred_Above[i,steps]

f0_Above[i] <- 1/esw_Above[i]

#Ones Trick
L.f0_Above[i] <- (1-exp(-((dist_Above[i]/ mu.df_Above[i])^(-b_Above)))) * 1/esw_Above[i]  # hazard rate likelihood

phi_Above[i]<- L.f0_Above[i]/C
ones.dist_Above[i] ~ dbern(phi_Above[i])
   
pcd_A[i]<-esw_Above[i]/W_Above

 
#MR part of MRDS

    
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
     
p_0[i]<-p0.1[i]+p0.2[i]-p0.1[i]*p0.2[i]
     
Abundance_A[i]<-gs[i]/(pcd_A[i]*p_0[i]) 
        
}

mean.esw_Above <- mean(esw_Above[])
   
#Calculate p(0)
P0_Above<-mean(p_0[])
 
N_Above<-sum(Abundance_A[])
N_Above_total<- N_Above/Av
  

###############################
#Loop through acoustic observations
#################################

for( i in 1:Nsights_Below){

mu.df_Below[i] <- b.df.0_Below

#Numerically integrate hazard rate function
esw_pred_Below[i,1]<-0

for (c in 1:(steps-1)){
 dx1_Below[i,c]<- 1-exp(-((y.step_Below[c]/mu.df_Below[i])^(-b_Below)))
 dx2_Below[i,c]<- 1-exp(-((y.step_plus1_Below[c]/mu.df_Below[i])^(-b_Below)))
 esw_pred_Below[i,c+1]<-esw_pred_Below[i,c]+0.5*step.size_Below*(dx1_Below[i,c]+dx2_Below[i,c])
     }

esw_Below[i]<-esw_pred_Below[i,steps]


#Ones Trick
L.f0_Below[i] <- (1-exp(-((dist_Below1[i]/ mu.df_Below[i])^(-b_Below)))) * 1/esw_Below[i]  # hazard rate likelihood

phi_Below[i]<- L.f0_Below[i]/C
ones.dist_B[i] ~ dbern(phi_Below[i])
   
pcd.B[i]<-esw_Below[i]/W_Below

Abundance_B[i]<-1/pcd.B[i] 
}

mean.esw_Below <- mean(esw_Below[])

N_Below<-sum(Abundance_B[])
   

       
##############################################################          
#CMR analysis of click trains
############################################################


for(i in 1:M){
logit(p[i,1])<-b0[1]+b.dist*dist_Below[i]+b.dist2*dist_Below[i]^2   
z[i,1]~dbern(gamma1)
mu[i,1]<-z[i,1]*p[i,1]
y[i,1]~dbern(mu[i,1])
y.rep[i,1]~dbern(mu[i,1])
recruitable[i,1]<-1
everRecruit[i]<-1-z[i,1]

Eval[i,1]<-mu[i,1]
resid.true[i,1]<-pow(y[i,1]-Eval[i,1],2)/(Eval[i,1]+0.5)
resid.pred[i,1]<-pow(y.rep[i,1]-Eval[i,1],2)/(Eval[i,1]+0.5)

for(j in 2:T){
logit(p[i,j])<-b0[j]+b.dist*dist_Below[i]+b.dist2*dist_Below[i]^2  #+b.dist2*dist_Below[i]^2+e[i]
logit(phi[i,j-1])<-b0.phi+b.age*age[i,j-1]  #+b.phi.dist*dist_Below[i]+b.phi.dist2*dist_Below[i]^2
survived[i,j]<- (phi[i,j-1]*z[i,j-1]) #*(1-everRecruit[i])+ everRecruit[i]*z[i,j-1]
recruitable[i,j]<- recruitable[i,(j-1)]*(1-z[i,j-1])
mu2[i,j]<-  survived[i,j] + gamma*recruitable[i,j]
z[i,j]~dbern(mu2[i,j])
mu3[i,j]<-z[i,j]*p[i,j]
y[i,j]~dbern(mu3[i,j])
y.rep[i,j]~dbern(mu3[i,j]) #For Bayesian G.O.F
 
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
D[j]<-sum(recruit[1:M,j])  #Enter dive state
S[j]<-sum(emigrate[1:M,j]) #Enter surface state
}



for(i in 1:M){
Nind[i]<-sum(z[i,1:end.interval])
Nalive[i]<-1-equals(Nind[i],0)
}

#Only count duplicates that occur in interval of overlap
totD<-sum(D[D.start.bin:D.end.bin])
totS<-sum(S[S.start.bin:S.end.bin])

Nsuper<-sum(Nalive[1:M])

#Average entering dive state per distance bin
meanD<-totD/Window

#Average entering surface state per distance bin
meanS<-totS/Window

#Per capita entering dive state per distance bin
pD<- meanD/Nsuper

#Percent emigrate per interval
pS<- meanS/Nsuper

#Calculate duplicates
duplicatesD<-N_Below*pD*WindowD #Initially Above then Below
duplicatesS<-N_Below*pS*WindowS #Initially Below then Above

duplicates<-duplicatesD+duplicatesS #

#Total abundance adjusted for duplicates
N_Total<- N_Below+N_Above-duplicates

#Availability bias
availability<-N_Above/N_Total

fit.true<-sum(resid.true)
fit.pred<-sum(resid.pred)
#
 } #end model


",fill = TRUE)
sink()

# Bundle data
Model.data <- list(y = y,
  C=10000,
  W_Above=W_Above/1000,
  W_Below=W_Below/1000,
  W_Below.CH=W_Below.CH/1000,
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
  z.mat=z.mat,
  D.start.bin=D.start.bin,
  D.end.interval=D.end.bin,
  S.start.bin=B.start.bin,
  S.end.bin=S.end.bin,
  start.bin=start.bin,
  end.bin=end.bin,
  a_Av=a_Av,
  b_Av=b_Av,
  Obs_Above=Dat_Sightings[,7:9],
  dist_Above=c(Dat_Sightings[,3])/1000+0.00001, #Convert to KM
  gs=Dat_Sightings[,6],
  Nsights_Above=Nsights_Above,
  Nsights_Below=Nsights_Below,
  Window=Window,
   WindowD=WindowD,
   WindowS=WindowS,
  age=age,
  dist_Below=dist_Below/1000,  #Convert to KM
  dist_Below1=dist_Below1/1000 #Convert to KM
  )

# Initial values
inits <- function(){list(b0 = runif(T, 0, 2), b.dist = runif(1, 0, 1),b.dist2 = runif(1, -1, 0),  b0.phi = runif(1, 0, 1),
                     b.age = runif(1, -1, 0),gamma= runif(1, 0.01, 0.07),  gamma1= runif(1, 0.1, 0.6),z=z.mat, sigma.e=1,
                      b.phi.dist = runif(1, -1, 1),b.phi.dist2 = runif(1, -1, 1),ln_b.df.0_Below=log(1.5), ln_b_Below=log(2),
                      ln_b.df.0_Above=log(1.5), ln_b_Above=log(2))}


varsToMonitor<-c("b0", "b.dist","b.dist2", "b.df.0_Above", "N_Above", "N_Total", "mean.esw_Above", "mean.esw_Below",      
                 "b_Below", "P_Below", "z","availability", "duplicates","b_Above", "Nsuper", "D", "S", "N", "phi", "p", "gamma","gamma1",
                 "totD", "totS", "N_Above_total", "P0.1_Above", "P0.2_Above", "P0_Above", "b1.0_Above","b2.0_Above","b.dist1_Above",
                  "b.dist2_Above", "b0.phi", "b.age", "fit.true", "fit.pred", "N_Below","b.df.0_Below",
                  "meanD","meanS","pD","pS","duplicatesD","duplicatesS","Av")


    # MCMC settings
niter <- nitt
nthin <- thin
nburn <-burnin
nchains <- chains



(before <- Sys.time())
HMM<- jags.model(
    file="./SA.jags",
    data = Model.data,
    inits = inits,
    n.chains = nchains,
    n.adapt = nburn
)
(mid1 <- Sys.time())

(mid2 <- Sys.time())
out <- jags.samples(
    model = HMM,
    variable.names = varsToMonitor,
    n.iter = niter,
    thin = nthin,
    progress.bar = 'text'
  )
(after <- Sys.time())

(diff(c(before,mid1,mid2,after)))[c(1,3)]

             
  

 


 