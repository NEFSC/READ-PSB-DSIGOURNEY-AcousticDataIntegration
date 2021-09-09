

 #JAGS Code Availabiity Model

  Hybrid_model<-function(nitt,thin,burnin,chains,W_Above,W_Below,ones.dist_Above, steps,y.step_Below,y.step_plus1_Below,step.size_Below,
  y.step_Above,y.step_plus1_Above,step.size_Above,y,z.mat,M,T,start.interval,end.interval,a_Av,b_Av,Nsights_Above,age,Y_I,
  dist_Below,Nsights_Below,dist_Below1,ones.dist_B,W_Below.CH,Window,Window.B,Window.E){


 sink("SA.jags")
cat("
model {


#Priors for hazard rate (DS model)

 #Visual data
  ln_b.df.0_Above ~ dnorm(0,0.0001)T(-2,2) #Scale parameter
  ln_b_Above ~ dnorm(0,0.0001)T(-2,2) #Shape parameter
  
  b.df.0_Above<-exp(ln_b.df.0_Above)
  b_Above<-exp(ln_b_Above)
  
  #PAM data             
  ln_b.df.0_Below ~ dnorm(0,0.0001)T(-2,2) 
  ln_b_Below ~ dnorm(0,0.0001)T(-2,2)
  
  b.df.0_Below<-exp(ln_b.df.0_Below)
  b_Below<-exp(ln_b_Below)

       
#Priors for MR Component of MRDS model for visual data
   b1.0_Above ~ dnorm(0,0.0001)  #Groups 1 and 2
   b2.0_Above ~ dnorm(0,0.0001)  #Sperm whale model (Group 3)
   b.dist1_Above~ dnorm(0,0.0001)  #Sperm whale model (Group 3)
   b.dist2_Above~ dnorm(0,0.0001)  #Sperm whale model (Group 3)
  

          
 #Informed prior for availability bias         
  Av~dbeta(a_Av, b_Av)

for( i in 1:Nsights_Above){


#################################
 #DS Part of MRDS  (Visual data)
#################################   
  
 #Use hazard rate funtion
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

 #Ones Trick to fit hazard rate function
   L.f0_Above[i] <- (1-exp(-((dist_Above[i]/ mu.df_Above[i])^(-b_Above)))) * 1/esw_Above[i]  # hazard rate likelihood


   phi_Above[i]<- L.f0_Above[i]/C
   ones.dist_Above[i] ~ dbern(phi_Above[i])
   
   pcd_A[i]<-esw_Above[i]/W_Above

################################   
    #MR part of MRDS
################################
    
   logit(p1_Above[i])<-b1.0_Above+b.dist1_Above*dist_Above[i] #Team 1           
   logit(p2_Above[i])<-b2.0_Above+b.dist2_Above*dist_Above[i] #Team 2
   
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
     
     
   Abundance_A[i]<-gs[i]/(pcd_A[i]*p_0[i])*Y_I_A[i]
        
}

  mean.esw_Above <- mean(esw_Above[])
   

  N_Above<-sum(Abundance_A[]) #Visual abundance estimate
  N_Above_total<- N_Above/Av  #Adjust visual abundance estimate for surface availabilty
  
    
#Fit DS model to PAM data
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

  f0_Below[i] <- 1/esw_Below[i]

 #Ones Trick
   L.f0_Below[i] <- (1-exp(-((dist_Below1[i]/ mu.df_Below[i])^(-b_Below)))) * 1/esw_Below[i]  # hazard rate likelihood


   phi_Below[i]<- L.f0_Below[i]/C
   ones.dist_B[i] ~ dbern(phi_Below[i])
   
   pcd.B[i]<-esw_Below[i]/W_Below

  Abundance_B[i]<-1/pcd.B[i]*Y_I_B[i] #Abundace estimate fro PAM data
}

  mean.esw_Below <- mean(esw_Below[])
  N_Below<-Nsights_Below*W_Below/mean.esw_Below  #Acoustic estimate only
  
       
##############################################################          
#Fit CMR model to PAM capture histories
############################################################

   gamma1~dunif(0,1)  

 for(j in 2:(T-1)){
    b0[j]~dnorm(0,0.0001)
   }
   
  for(j in 1:(T-1)){
    gamma[j]~dunif(0,1)  
   }
  
 #Constraints  
   b0[1]<-b0[2]
   b0[T]<-b0[T-1]

   b.dist~dnorm(0,0.0001)
   b.dist2~dnorm(0,0.0001)
   
    b0.phi~dnorm(0,0.0001)
    b.age~dnorm(0,0.0001) 
    b.phi.dist~dnorm(0,0.0001) 
    b.phi.dist2~dnorm(0,0.0001) 

  #Distribution for unobserved agaes
   for(i in 1:M){
     for(j in 1:T){
        age[i,j]~dunif(0,j)
        }
         }


  for(i in 1:M){
    dist_Below[i]~dunif(0,W_Below.CH) #Distribution for unobserved distances
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
   mu2[i,j]<-  survived[i,j] + gamma[j-1]*recruitable[i,j]
   z[i,j]~dbern(mu2[i,j])
   mu3[i,j]<-z[i,j]*p[i,j]
   y[i,j]~dbern(mu3[i,j])
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
   B[j]<-sum(recruit[1:M,j])
   Em[j]<-sum(emigrate[1:M,j])
     }



  for(i in 1:M){
   Nind[i]<-sum(z[i,1:end.interval])
   Nalive[i]<-1-equals(Nind[i],0)
  }


#Only count duplicates that occur in interval of overlap
    totEm<-sum(Em[start.interval:end.interval])
    totRec<-sum(B[start.interval:end.interval])



Nsuper<-sum(Nalive[1:M])

#Average recruit per interval
mean.Rec<-totRec/Window

#Average emigrant per interval
mean.Em<-totEm/Window

#Percent recruit per interval
p.Rec<- mean.Rec/Nsuper

#Percent emigrate per interval
p.Em<- mean.Em/Nsuper

#Duplicates<-N_Below*gamma*Window
  duplicates.Rec<-N_Below*p.Rec*Window.B #Initially above then below
  duplicates.Em<-N_Below*p.Em*Window.E #Initially Below then Above

#Total duplicates
  duplicates<-duplicates.Rec+duplicates.Em

#Total abundance adjusted for duplicates
  N.total<- N_Below+N_Above-duplicates

#Availability bias
availability<-N_Above/N.total
availability2<-N_Above/N_Below

 #For Bayesian G.O.F.  
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
  #em.start.interval=em.start.interval,
#  em.end.interval=em.end.interval,
#  B.start.interval=B.start.interval,
#  B.end.interval=B.end.interval,
  start.interval=start.interval,
  end.interval=end.interval,
  a_Av=a_Av,
  b_Av=b_Av,
  Obs_Above=Dat_Sightings[,7:9],
  dist_Above=c(Dat_Sightings[,3])/1000+0.00001,
  gs=Dat_Sightings[,6],
  Nsights_Above=Nsights_Above,
  Nsights_Below=Nsights_Below,
  Window=Window,
   Window.B=Window.B,
   Window.E=Window.E,
  age=age,
  Y_I_A=Y_I_A,
  Y_I_B=Y_I_B,
  dist_Below=dist_Below/1000,
  dist_Below1=dist_Below1/1000
  )

# Initial values
  inits <- function(){list(b0 = runif(T, 0, 2), b.dist = runif(1, 0, 1),b.dist2 = runif(1, -1, 0),  b0.phi = runif(1, 0, 1),
                     b.age = runif(1, -1, 0),gamma= runif(T-1, 0.01, 0.07),  gamma1= runif(1, 0.1, 0.6),z=z.mat, sigma.e=1,
                      b.phi.dist = runif(1, -1, 1),b.phi.dist2 = runif(1, -1, 1),ln_b.df.0_Below=log(1.5), ln_b_Below=log(2),
                      ln_b.df.0_Above=log(1.5), ln_b_Above=log(2))}


  varsToMonitor<-c("b0", "b.dist","b.dist2", "b.df.0_Above", "N_Above", "N.total", "mean.esw_Above", "mean.esw_Below",      
                 "b_Below", "P_Below", "z","availability", "duplicates","b_Above", "Nsuper", "B", "N", "phi", "p", "gamma","gamma1","Em",
                 "totEm", "totRec", "N_Above_total", "P0.1_Above", "P0.2_Above", "P0_Above", "b1.0_Above","b2.0_Above","b.dist1_Above",
                  "b.dist2_Above", "b0.phi", "b.age", "sigma.e", "fit.true", "fit.pred", "b.phi.dist", "b.phi.dist2","N_Below","b.df.0_Below",
                  "mean.Rec","mean.Em","p.Rec","p.Em","duplicates.Rec","duplicates.Em","availability2")


    #MCMC settings
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

                out
              

 } #end

