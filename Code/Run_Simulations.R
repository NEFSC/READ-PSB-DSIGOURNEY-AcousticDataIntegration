



#setwd("C:/Acoustic Analysis/Final Report Summer 2020/Code")

library(plyr)
library(dplyr)
library(tidyr)
library(coda)
library(rjags)
library(reshape)
library(here)

#Find working directories
here()

#Choose number of simulations
sims=1
  
#Total number of whales
N.true<-600
N=N.true 

#Choose subset of click trains for Hybrid Method
n.sub=35
  
#Parameters to simulate detection for acoustic (below) data 
b0_b=6
b.dist_b=-1.5
  
#Half-normal parametrs for visual and below detection
b0.A=1
sigma.A<-exp(b0.A)

b0.B=1  #for testing Hybrid Model
sigma.B<-exp(b0.B)

#Truncation distances 
W_Below<-6 #Truncation for acoustic/below (currently 6 KM)
W_Above<-6 #Truncation for visual/above (currently 6 KM)


#Pick maximum forward distance
initial.y <- W_Below-0.1 #Forward distance at start of simulation(adjusting for error)
 
#
y_max=W_Below
y_min=-5
y_width<-0.5  #width of forward distance bins
  
y_bins<-rev(seq(y_min,y_max,y_width))

n_bins<-length(y_bins)-1

#For CMR-DS analysis
Bin.start<-6
Bin.end<-17
       
S.start<-Bin.start+1
S.end<-9
        
D.start<-Bin.start+1
D.end<-Bin.end

#For Hybrid analysis
Bin.start.H<-6
Bin.end.H<-17
       
start.H<-7
end.H<-20
Window.H<-end.H+1-start.H
         
WindowS.H=3
WindowD.H=11
    
summary_data<-matrix(NA,nrow=sims,ncol=40)
summary_data<-as.data.frame(summary_data)
    colnames(summary_data)<-c("Sim","W_Below","W_Above","y_max","y_min","y_width","N.true", "N.tot.Below", "N.tot.Above", 
     "Availability", "N.Total_CMR_DS", "CV_N.Total_CMR_DS","Cov_N.Total_CMR_DS", "Availability_CMR_DS","CV_Availability_CMR_DS",
     "Cov_Availability_CMR_DS", "N.tot.Below_CMR_DS","CV_N.tot.Below_CMR_DS","Cov_N.tot.Below_CMR_DS", 
     "N.tot.Above_CMR_DS","CV_N.tot.Above_CMR_DS","Cov_N.tot.Above_CMR_DS","N.Total_Hybrid","CV_N.Total_Hybrid", 
     "Cov_N.Total_Hybrid","Availability_Hybrid", "CV_Availability_Hybrid", "Cov_Availability_Hybrid","N.tot.Below_Hybrid","CV_N.Below_Hybrid",
     "Cov_N.Below_Hybrid","N.Above_Hybrid","CV_N.Above_Hybrid","Cov_N.Above_Hybrid", "N.Total_DS_DS", "CV_N.Total_DS_DS","Cov_N.Total_DS_DS",
     "Availability_DS_DS","CV_Availability_DS_DS","Cov_Availability_DS_DS")

#Run simulations and save results 
for (s in 1:sims){
       
source("Simulate_Data.R") #Simulate acoustic nd visual datasets
   
source("Organize_Simulation_Data.R") #Organize data forinout into model 
   
source("CMR_DS_Model.R")   #Analyze data with CMR model  
   
source("Hybrid_Model.R") #Analyze data with Hybrid model
   
                  
#Summarize results for simulation s               
  summary_data[s,1]<-s
  summary_data[s,2]<-W_Below
  summary_data[s,3]<-W_Above
  summary_data[s,4]<-y_max
  summary_data[s,5]<-y_min
  summary_data[s,6]<-y_width
  summary_data[s,7]<-N.true 
  summary_data[s,8]<-N.tot.Below
  summary_data[s,9]<-N.tot.Above  
  summary_data[s,10]<-Availability #True surface availability                                       
  
  summary_data[s,11]<-N.total_CMR_DS  
  summary_data[s,12]<-cv_N.total_CMR_DS 
  summary_data[s,13]<-coverage_N.total_CMR_DS   
   
  summary_data[s,14]<-availability_CMR_DS  
  summary_data[s,15]<-cv_availability_CMR_DS  
  summary_data[s,16]<-coverage_availability_CMR_DS    
  
  summary_data[s,17]<-N.tot.Below_CMR_DS  
  summary_data[s,18]<-cv_N.tot.Below_CMR_DS
  summary_data[s,19]<-coverage_N.tot.Below_CMR_DS
  
  summary_data[s,20]<-N.tot.Above_CMR_DS  
  summary_data[s,21]<-cv_N.tot.Above_CMR_DS
  summary_data[s,22]<-coverage_N.tot.Above_CMR_DS
   
  summary_data[s,23]<-N.total_Hybrid  
  summary_data[s,24]<-cv_N.total_Hybrid 
  summary_data[s,25]<-coverage_N.total_Hybrid   
   
  summary_data[s,26]<-availability_Hybrid  
  summary_data[s,27]<-cv_availability_Hybrid  
  summary_data[s,28]<-coverage_availability_Hybrid    
  
  summary_data[s,29]<-N.tot.Below_Hybrid  
  summary_data[s,30]<-cv_N.tot.Below_Hybrid
  summary_data[s,31]<-coverage_N.tot.Below_Hybrid
  
  summary_data[s,32]<-N.tot.Above_Hybrid  
  summary_data[s,33]<-cv_N.tot.Above_Hybrid
  summary_data[s,34]<-coverage_N.tot.Above_Hybrid
   
  summary_data[s,35]<-N.total_DS_DS
  summary_data[s,36]<-cv_N.total_DS_DS 
  summary_data[s,37]<-coverage_N.total_DS_DS   
   
  summary_data[s,38]<-availability_DS_DS  
  summary_data[s,39]<-cv_availability_DS_DS  
  summary_data[s,40]<-coverage_availability_DS_DS
      
            }
            
           
   #       write.csv(summary_data, file=paste0("summary_sims_","comparison","_",Sim.Scenario,"_",run,"_6000_gamma_t",".csv",sep=""))    
   
                             
   
           
       
             
                                     
