

#Load libraries

library(rjags)
library(reshape)
library(dplyr)
library(here)


#Set mcmc parameters
nitt=15000
burnin=10000
chains=1                  
thin=15
     

########################
#Inputs for Hybrid analysis
#########################

Bin.start.H<-1
Bin.end.H<-25

#Choose bins to calculate mean Fj and Sj       
start.H<-15
end.H<-29
Zone.H<-end.H+1-start.H
#Zone.H<-25
#Zone.H<-Bin.end.H+1-Bin.start.H
         
ZoneF.H=22
ZoneS.H=10

#Set up forward distance bins
Y_max<-10000 
y_max<- -Y_max 
y_min<-10000

y_width<-309  #width of forward distance bins (meters)
y_bins<-seq(y_max,y_min ,y_width)
y_bins<-y_bins*-1
y_mid.vec_all<-y_bins-y_width/2   
  
#############################
#Import and format datasets
############################

###################################################
#Import capture history (CH) matrix
CH_mat_all<-read.csv(here::here("Data", "CH_Matrix.csv")) 
dist_Below<-CH_mat_all$dist_below #For modelling detetcion within CMR model
ncols<-dim(CH_mat_all)[2]
CH_mat<-CH_mat_all[,4:ncols]

#Truncate total columns in CH matrix
Col.start<-18
Col.end<-48
CH_mat<-CH_mat[,Col.start:Col.end]
ch.indicator<-which(rowSums(CH_mat)>0) #Was individual detected within truncated CH matrix
CH_mat<-CH_mat[ch.indicator,]
#Truncate distance bins
y_mid.vec<-y_mid.vec_all[Col.start:Col.end]

dist_Below<-dist_Below[ch.indicator]
W_Below.CH<-max(dist_Below,na.rm=TRUE)
    
#Data augmentation step
y<-as.matrix(CH_mat)
nz<-60
y<-rbind(y,matrix(rep(0,ncol(y)*nz),ncol=ncol(y),byrow=TRUE))
M<-nrow(y)
T<-ncol(y)
z.mat<-matrix(1,nrow=dim(y)[1], ncol=dim(y)[2])

dist_Below<-c(dist_Below,rep(NA,nz))
############################################################


############################################################
#Load and format sperm whale VLT data for MRDS analysis
VLT_data<-read.csv(here::here("Data","sperm_whale_VLT_data.csv"))

#Format data for input into a Bayesian MRDS analysis   
source(here::here("Code","format_MRDS_data.R"))

Nsights_Above<-dim(VLT_MRDS)[1]
ones.dist_Above<-rep(1,Nsights_Above)

#Truncation distance for VLT 
W_Above=7600

#For integrating hazard rate model
steps<-50
y.step_all_Above<-seq(0.00001,W_Above/1000, length=steps)
y.step_Above<-y.step_all_Above[1:(steps-1)]
y.step_plus1_Above<-y.step_all_Above[2:steps]
step.size_Above<-(W_Above/1000)/(steps-1)

#Informed prior for availability bias   
Av_est=0.613 
Av_CV=0.247  

c_Av=((Av_est*(1-Av_est))/((Av_est*Av_CV)^2))-1
a_Av=c_Av*Av_est
b_Av=c_Av*(1-Av_est)  

###########################################################


##########################################################
#Load and format sperm whale PAM data for DS analysis
summary_mat1<-read.csv(here::here("Data", "sperm_whale_PAMLT_data.csv"))
summary_mat1<-subset(summary_mat1, meanX<W_Below)
dist_Below1<-summary_mat1$meanX
Nsights_Below<-dim(summary_mat1)[1]
ones.dist_B<-rep(1,Nsights_Below)  
 
W_Below=W_Above #Same as VLT 
   
#For integrating hazard rate model
steps<-50
y.step_all_Below<-seq(0.00001,W_Below/1000, length=steps)
y.step_Below<-y.step_all_Below[1:(steps-1)]
y.step_plus1_Below<-y.step_all_Below[2:steps]
step.size_Below<-(W_Below/1000)/(steps-1)

###############################################################

    
#Run Hybrid analysis and summarize results           
source(here::here("Code","Hybrid_sperm_whale_analysis.R")

#####################
#Summarize output
#####################

#Abundance estimate Hybrid
abun.H<-quantile(out$N.total[1,,1],0.5)
abun.H.upper<-quantile(out$N.total[1,,1],0.975)
abun.H.lower<-quantile(out$N.total[1,,1],0.025)

#Aavailability bias estimate Hybrid
avail.H<-quantile(out$availability[1,,1],0.5)
avail.H.upper<-quantile(out$availability[1,,1],0.975)
avail.H.lower<-quantile(out$availability[1,,1],0.025)

