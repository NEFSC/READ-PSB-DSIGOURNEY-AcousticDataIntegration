

library(rjags)
library(dplyr)
library(here)


#setwd
setwd("C:/Acoustic Analysis/HMM Simulations Test/Final SPWH Analysis")


#Truncate CH matrix
CH.start<-18
CH.end<-48

#Start and and bins  for window of overlap  
start.bin<-1
end.bin<-25

Window=end.bin-(start.bin-1)

#Bins for calculating the average number entering the dive state  
D.start.bin<-7
D.end.bin<-16

#Bins for calculating the average number entering the surfacing state  
S.start.bin<-7
S.end.bin<-16

#Size of Window where dupilcates can occur
Window.D=22
Window.S=10

#Truncation distances for Distance Analysis (meters)         
W_Above=7600
W_Below=W_Above

##############################################    
#Import Data

#Read in capture history (CH) )matrix
All_dat<-read.csv("Final_CH_Matrix.csv") 
ncols<-dim(All_dat)[2]
detected_mat<-All_dat[,4:ncols] #CH Matrix
dist_Below<-All_dat$dist_below #Visual Perpendicular distances for MRDS
  
#Import visual data     
Vis_data<-read.csv("Dat_Sightings_SPWH.csv")
Nsights_Above<-dim(Vis_data)[1] 
Dat_Sightings<-Vis_data
   
#Import passive acoustic data
summary_mat1<-read.csv("Final_Acoustic_Input_For_CDS_Analysis_test1.csv")
summary_mat1<-subset(summary_mat1, meanX<W_Below)
dist_Below1<-summary_mat1$meanX  #Acoustic perpendicular distances for DS analysis
Nsights_Below<-dim(summary_mat1)[1]

 
############################################################
#For numerical integration of hazard rate function
steps<-50

#For passive acoustic data
y.step_all_Below<-seq(0.00001,W_Below/1000, length=steps)
y.step_Below<-y.step_all_Below[1:(steps-1)]
y.step_plus1_Below<-y.step_all_Below[2:steps]
step.size_Below<-(W_Below/1000)/(steps-1)

#For visual data
y.step_all_Above<-seq(0.00001,W_Above/1000, length=steps)
y.step_Above<-y.step_all_Above[1:(steps-1)]
y.step_plus1_Above<-y.step_all_Above[2:steps]
step.size_Above<-(W_Above/1000)/(steps-1)
#################################################################


##############################################################
#Make age matrix

get.first<-function(x)min(which(x!=0))
get.last<-function(x)max(which(x!=0))
  
f<-apply(detected_mat,1,get.first)
g<-apply(detected_mat,1,get.last)
n<-nrow(detected_mat)
t<-ncol(detected_mat)
         
age<-matrix(0,nrow=n,ncol=t) 
 for (i in 1:n){
   for (j in f[i]:t){
     age[i,j]<-j-f[i]+1
            }
             }

zz<-matrix(0,nrow=n,ncol=t) 
 for (i in 1:n){
    for (j in f[i]:g[i]){
       zz[i,j]<-1
            }
            }
            
########################
#Subset ch matrix

detected_mat<-detected_mat[,CH.start:CH.end]
ch.indicator<-which(rowSums(detected_mat)>0)
detected_mat<-detected_mat[ch.indicator,] 
age<-age[ch.indicator,CH.start:CH.end]
dist_Below<-dist_Below[ch.indicator]

W_Below.CH<-max(dist_Below,na.rm=TRUE)  #Max truncation for CMR data

###############################################################
#For Data Augmentation  step
y<-as.matrix(detected_mat)
nz<-60  #number of extra rows to augment
y<-rbind(y,matrix(rep(0,ncol(y)*nz),ncol=ncol(y),byrow=TRUE))
M<-nrow(y)
T<-ncol(y)

z.mat<-matrix(1,nrow=dim(y)[1], ncol=dim(y)[2])

dist_Below<-c(dist_Below,rep(NA,nz))
age<-rbind(age,matrix(rep(NA,ncol(age)*nz),ncol=ncol(age),byrow=TRUE))
    
#For JAGS "ones" trick to fit hazard rate likelihood 
ones.dist_Above<-array(NA,Nsights_Above)

for (i in 1:Nsights_Above){
     ones.dist_Above[i]<-1
     }
     
ones.dist_B<-array(NA,Nsights_Below)  

for (i in 1:Nsights_Below){
     ones.dist_B[i]<-1
     }
    

     
#Set up informative Beta prior for availability bias
Av_est=0.613 #estimate of surface availability for sperm whales from Palka et al. (2017)
Av_CV=0.247  #sperm whale (shipboard)
c_Av=((Av_est*(1-Av_est))/((Av_est*Av_CV)^2))-1
a_Av=c_Av*Av_est
b_Av=c_Av*(1-Av_est)  
           
W_Below.CH<-max(dist_Below,na.rm=TRUE)
 
#MCMC inputs          
nitt=30000
burnin=20000
chains=2                  
thin=3
 


#Run Model 
 
source("Hybrid_model_for_sperm_whales.R")
  
#save( Output_All_dot, file="test3_Final_Redo_3.RData")
  
   
   
   