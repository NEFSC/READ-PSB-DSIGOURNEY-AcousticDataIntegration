#######################
#Script to run simulations decsribed in Sigourney et al.
###############################################


library(plyr)
library(dplyr)
library(tidyr)
library(coda)
library(rjags)
library(reshape)
library(here)
library(furrr)
  
#For parallel processing using furr package  
plan("multisession") 

#Find working directories
here()

#Choose number of simulations
sims=100


#Set mcmc parameters
burnin=5000
nitt=15000
chains=1
thin=15
  
#Total number of all whales
N.true<-400

#Choose subset of click trains for Hybrid Method
n.sub=35

#set maxTime
maxTime=35
  
#Parameters to simulate detection for PAM data 
b0_b=7
b.dist_b=-1.0
  
#Half-normal parameters for VLT and PAM detection
#VLT (A for Above)
b0.A=1
sigma.A<-exp(b0.A)

#PAM (B for Below)
b0.B=1  
sigma.B<-exp(b0.B)

#Truncation distances 
W_Below<-7 #Truncation for acoustic/below (currently 7 KM)
W_Above<-5 #Truncation for visual/above (currently 5 KM)


#Pick maximum forward distance
initial.y <- W_Below-0.001 #-0.1 #Forward distance at start of simulation(adjusting for error)
 
#
y_max=W_Below
y_min=-5
y_width<-0.25  #width of forward distance bins
  
y_bins<-rev(seq(y_min,y_max,y_width))

n_all<-length(y_bins)
n_bins<-length(y_bins)-1

Y1<-y_bins[1:n_bins]
Y2<-y_bins[2:n_all]
YMid.vec=(Y1+Y2)/2

##########################
#For CMR-DS analysis
########################

#Zone of overlap with VLT (starts at y= 5 KM)
Bin.start<-9   #
Bin.end<-28

#Zone of overlap for surfacing        
S.start<-Bin.start+1
S.end<-Bin.end-9

#Zone of overlap for foraging       
F.start<-Bin.start+4  
F.end<-Bin.end         

########################
#For Hybrid analysis
#########################
Bin.start.H<-Bin.start
Bin.end.H<-Bin.end

#Choose bins to calculate mean Fj and Sj       
start.H<-7
end.H<-20
Zone.H<-end.H+1-start.H
         
ZoneS.H=S.end-S.start
ZoneF.H=F.end-F.start


#Source functions 
source(here::here("Code","simulate_data.R"))
source(here::here("Code","sum_sims.R"))
source(here::here("Code","make_inputs.R"))
source(here::here("Code","CMR_DS.R"))
source(here::here("Code","Hybrid.R"))


Run.sims<-function(s){

#simulate data       
sims<-simulate_data(N.all=N.true,W_Below=W_Below,
           W_Above=W_Above,b0_b=b0_b,b.dist_b=b.dist_b,sigma.A=sigma.A,
           sigma.B=sigma.B)

Final_Dat<-sims$Final_Dat
Below_Dat<-sims$Below_Dat
Above.obs<-sims$Above.obs
Nobs_Above<-sims$Nobs_Above
Below.obs<-sims$Below.obs
Nobs_Below<-sims$Nobs_Below

#summarize simulations
sum.sims<-sum_sims(Final_Dat=Final_Dat,W_Above=W_Above)
N.true.Below<-sum.sims$N.true.Below
N.true.Above<-sum.sims$N.true.Above
Availability<-sum.sims$Availability
Dups_click<-sum.sims$Dups_click
 
#Organize data for input into models
input_data<-make_inputs(Below_Dat=Below_Dat,Above.obs=Above.obs,
             Nobs_Above=Nobs_Above,Below.obs=Below.obs,Nobs_Below=Nobs_Below)


             
y<-input_data$y
M<-input_data$M
T<-input_data$T
z.mat<-input_data$z.mat
y.sub<-input_data$y.sub
M.sub<-input_data$M.sub
T.sub<-input_data$T.sub
z.mat.sub<-input_data$z.mat.sub
dist_Above<-input_data$dist_Above
dist_Below1<-input_data$dist_Below1
dist_Below2<-input_data$dist_Below2
dist_Below3<-input_data$dist_Below3
ones.dist_A<-input_data$ones.dist_A
ones.dist_B<-input_data$ones.dist_B

#######################
#Run CMR_DS Model  
 #######################
CMR_DS<- CMR_DS_model(nitt=nitt,thin=thin,burnin=burnin,chains=chains,
dist_Above=dist_Above,Nobs_Above=Nobs_Above,W_Above=W_Above,
dist_Below1=dist_Below1,W_Below=W_Below,ones.dist_A=ones.dist_A,
y=y,z.mat=z.mat,M=M,T=T,Bin.start=Bin.start,
Bin.end=Bin.end,S.start=S.start,S.end=S.end,F.start=F.start,
F.end=F.end,YMid.vec=YMid.vec,maxTime=maxTime)      
                 
#Abundance Estimate
N.total_CMR_DS<-median(CMR_DS$N.total[1,,1])
sd_N.total_CMR_DS<-sd(CMR_DS$N.total[1,,1])
cv_N.total_CMR_DS<- sd_N.total_CMR_DS/N.total_CMR_DS
                      
#Availability Estimate 
availability_CMR_DS<-median(CMR_DS$availability[1,,1])
sd_availability_CMR_DS<-sd(CMR_DS$availability[1,,1])
cv_availability_CMR_DS<- sd_availability_CMR_DS/availability_CMR_DS               
                                   
#Duplicates
Dups_CMR_DS<-median(CMR_DS$duplicates[1,,1])
sd_Dups_CMR_DS<-sd(CMR_DS$duplicates[1,,1])
cv_Dups_CMR_DS<- sd_Dups_CMR_DS/Dups_CMR_DS
              
#Acoustic Estimate (i.e.,total below)
N.tot.Below_CMR_DS<-median(CMR_DS$N.tot.Below[1,,1])
sd_N.tot.Below_CMR_DS<-sd(CMR_DS$N.tot.Below[1,,1])
cv_N.tot.Below_CMR_DS<- sd_N.tot.Below_CMR_DS/N.tot.Below_CMR_DS

######################
#Run Hybrid Model (and DS-DS Model)
######################  
Hybrid<-Hybrid_model(nitt=nitt,thin=thin,burnin=burnin,chains=chains,W_Above=W_Above,W_Below=W_Below,Nobs_Above=Nobs_Above,
Nobs_Below=Nobs_Below,dist_Above=dist_Above,dist_Below2=dist_Below2,dist_Below3=dist_Below3,ones.dist_A=ones.dist_A,
ones.dist_B=ones.dist_B,y=y.sub,z.mat=z.mat.sub,M=M.sub,T=T.sub,Bin.start.H=Bin.start.H,Bin.end.H=Bin.end.H,start.H=start.H,
end.H=end.H,Zone.H=Zone.H,ZoneS.H=ZoneS.H,ZoneF.H=ZoneF.H)   
                
#Abundance Estimate
N.total_Hybrid<-median(Hybrid$N.total[1,,1])
sd_N.total_Hybrid<-sd(Hybrid$N.total[1,,1])
cv_N.total_Hybrid<- sd_N.total_Hybrid/N.total_Hybrid
         
#Availability Estimate 
availability_Hybrid<-median(Hybrid$availability[1,,1])
sd_availability_Hybrid<-sd(Hybrid$availability[1,,1])
cv_availability_Hybrid<- sd_availability_Hybrid/availability_Hybrid
             
#Duplicates
Dups_Hybrid<-median(Hybrid$duplicates[1,,1])
sd_Dups_Hybrid<-sd(Hybrid$duplicates[1,,1])
cv_Dups_Hybrid<- sd_Dups_Hybrid/Dups_Hybrid
                
#Acoustic Estimate (i.e.,total below)
N.tot.Below_Hybrid<-median(Hybrid$N.tot.Below[1,,1])
sd_N.tot.Below_Hybrid<-sd(Hybrid$N.tot.Below[1,,1])
cv_N.tot.Below_Hybrid<- sd_N.tot.Below_Hybrid/N.tot.Below_Hybrid
               
#Visual Estimate (i.e.,total above)
N.tot.Above_Hybrid<-median(Hybrid$N.tot.Above[1,,1])
sd_N.tot.Above_Hybrid<-sd(Hybrid$N.tot.Above[1,,1])
cv_N.tot.Above_Hybrid<- sd_N.tot.Above_Hybrid/N.tot.Above_Hybrid
                  
#Abundance Estimate DS-DS Method
N.total_DS_DS<-median(Hybrid$N.total2[1,,1])
sd_N.total_DS_DS<-sd(Hybrid$N.total2[1,,1])
cv_N.total_DS_DS<- sd_N.total_DS_DS/N.total_DS_DS
              
#Availability Estimate 
availability_DS_DS<-median(Hybrid$availability2[1,,1])
sd_availability_DS_DS<-sd(Hybrid$availability2[1,,1])
cv_availability_DS_DS<- sd_availability_DS_DS/availability_DS_DS             
                   
                      
summary_data<-matrix(NA,nrow=1,ncol=33)
summary_data<-as.data.frame(summary_data)
    colnames(summary_data)<-c("Sim","W_Below","W_Above","y_max","y_min","y_width","N.true", "N.tot.Below", "N.tot.Above", 
     "Availability", "Duplicates", "N.Total_CMR_DS", "CV_N.Total_CMR_DS", "Availability_CMR_DS","CV_Availability_CMR_DS",
     "Duplicates_CMR_DS","CV_Duplicate_CMR_DS", "N.Below_CMR_DS","CV_N.Below_CMR_DS","N.Total_Hybrid","CV_N.Total_Hybrid",
     "Availability_Hybrid", "CV_Availability_Hybrid", "Duplicates_Hybrid","CV_Duplicates_Hybrid","N.Below_Hybrid",
     "CV_N.Below_Hybrid","N.Above_Hybrid","CV_N.Above_Hybrid","N.Total_DS_DS", "CV_N.Total_DS_DS","Availability_DS_DS",
     "CV_Availability_DS_DS")

#Summarize results for simulations               
 
  summary_data[1,1]<-s #+(sims*run-sims)
  summary_data[1,2]<-W_Below
  summary_data[1,3]<-W_Above
  summary_data[1,4]<-y_max
  summary_data[1,5]<-y_min
  summary_data[1,6]<-y_width
  summary_data[1,7]<-N.true 
  summary_data[1,8]<-N.true.Below
  summary_data[1,9]<-N.true.Above  
  summary_data[1,10]<-Availability #True surface availability  
  summary_data[1,11]<-Dups_click   #True duplicates                                     
  
  summary_data[1,12]<-N.total_CMR_DS  
  summary_data[1,13]<-cv_N.total_CMR_DS 
   
  summary_data[1,14]<-availability_CMR_DS  
  summary_data[1,15]<-cv_availability_CMR_DS     
  
  summary_data[1,16]<-Dups_CMR_DS  
  summary_data[1,17]<-cv_Dups_CMR_DS  
  
  summary_data[1,18]<-N.tot.Below_CMR_DS  
  summary_data[1,19]<-cv_N.tot.Below_CMR_DS
     
  summary_data[1,20]<-N.total_Hybrid  
  summary_data[1,21]<-cv_N.total_Hybrid 
   
  summary_data[1,22]<-availability_Hybrid  
  summary_data[1,23]<-cv_availability_Hybrid  
  
  summary_data[1,24]<-Dups_Hybrid  
  summary_data[1,25]<-cv_Dups_Hybrid    
  
  summary_data[1,26]<-N.tot.Below_Hybrid  
  summary_data[1,27]<-cv_N.tot.Below_Hybrid
  
  summary_data[1,28]<-N.tot.Above_Hybrid  
  summary_data[1,29]<-cv_N.tot.Above_Hybrid
   
  summary_data[1,30]<-N.total_DS_DS
  summary_data[1,31]<-cv_N.total_DS_DS 
   
  summary_data[1,32]<-availability_DS_DS  
  summary_data[1,33]<-cv_availability_DS_DS  

   return(summary_data)      
            }
            
  new_summary<-future_map_dfr(1:sims,~Run.sims(.))
           
   write.csv(new_summary, "summary_sims.csv")    
   
                             
   
           
       
             
                                     
