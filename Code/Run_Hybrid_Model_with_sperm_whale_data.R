   
 
library(rjags)
library(reshape)
  
#Set working dir.  
setwd("C:/Acoustic Analysis")

#define truncation distances
 W_Above=7600 #Define truncation distance for visual data
 W_Below=W_Above
 
 
#define zone of overlap
 start.interval<-10
 end.interval<-25
 Window=end.interval-(start.interval-1)
 #Window.B=14
# Window.E=1
  Window.B=22
 Window.E=7
           

##############################################    
#Import Data

  #Import visual data     
     Vis_data<-read.csv("Dat_Sightings_SPWH.csv")
     Nsights_Above<-dim(Vis_data)[1] 
     Dat_Sightings<-Vis_data
   
 #Import passive acoustic data
  summary_mat1<-read.csv("Final_Acoustic_Input_For_CDS_Analysis_test1.csv")
  summary_mat1<-subset(summary_mat1, meanX<W_Below)
  dist_Below1<-summary_mat1$meanX
  Nsights_Below<-dim(summary_mat1)[1]
  Y_I_B<-as.numeric(summary_mat1$maxY>0)*as.numeric(summary_mat1$minY<8000)
 
 
     
  #Import subset of capture histories from  
     CH<-read.csv("detected_mat_spwh.csv")
     ncols.CH<-dim(CH)[2]
      dist_Below<-CH$Perp.Dist
      detected_mat<-CH[,3:ncols.CH]
     
     #NEED TO IMPORT CORESPNDING Perd Dist
     


############################################################
#For numerical integration of hazard rate function
 steps<-50

    #For passive acoustic data
     y.step_all_Below<-seq(0.00001,W_Below/1000, length=steps)
     y.step_Below<-y.step_all_Below[1:(steps-1)]
     y.step_plus1_Below<-y.step_all_Below[2:steps]
     step.size_Below<-(W_Below/1000)/(steps-1)

     #For passive acoustic data
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
    
    
 #Indicator for which visual sighting are inside zone of overlap (forward truncation distance)
 #    Y_I<-as.numeric(Dat_Sightings$Y<abs(y_max))
   Y_I_A<-as.numeric(Dat_Sightings$Y<8000)
     
#Set up informative Beta prior for availability bias
  Av_est=0.613 #estimate of surface availability for sperm whales from Palka et al. (2017)
  Av_CV=0.247  #sperm whale (shipboard)
  c_Av=((Av_est*(1-Av_est))/((Av_est*Av_CV)^2))-1
  a_Av=c_Av*Av_est
  b_Av=c_Av*(1-Av_est)  
           
         
           W_Below.CH<-max(dist_Below,na.rm=TRUE)
 
 #MCMC inputs          
 nitt=10000
 burnin=5000
 chains=1                  
 thin=10

#Call JAGS function
 source("Hybrid_model_for_sperm_whales.R")
 
#Run Model 
 Output_All<-Hybrid_model(nitt=nitt,thin=thin,burnin=burnin,chains=chains,
    W_Above=W_Above, ones.dist_Above= ones.dist_Above,W_Below=W_Below,
     steps=steps,y.step_Below=y.step_Below,y.step_plus1_Below=y.step_plus1_Below,
     step.size_Below=step.size_Below,y.step_Above=y.step_Above,y.step_plus1_Above=y.step_plus1_Above,
     step.size_Above=step.size_Above,y=y,z.mat=z.mat,M=M,T=T,start.interval=start.interval,
     end.interval=end.interval, a_Av=a_Av,b_Av=b_Av,Nsights_Above=Nsights_Above, age=age, Y_I=Y_I,
     dist_Below=dist_Below,Nsights_Below=Nsights_Below, dist_Below1=dist_Below1, ones.dist_B= ones.dist_B,
     W_Below.CH=W_Below.CH,Window=Window,Window.B=Window.B,Window.E=Window.E)
     
 
par(mfrow=c(1,2))   
hist(Output_All$N.total[1,,1], main="Abundance Estimate Method 1",xlab="Abundance")
abline(v=median(Output_All$N.total[1,,1]), col='black', lty=3,lwd=4)
hist(Output_All$N_Above_total[1,,1], main="Abundance Estimate CDS",xlab="Abundance")
abline(v=median(Output_All$N.total[1,,1]), col='black', lty=3,lwd=4)

hist(Output_All$availability[1,,1], main="Availability Estimate Method 1",xlab="Availability")
av1<-median(Output_All$availability[1,,1])
abline(v=av1, col='black', lty=3,lwd=4)
abline(v=0.61, col='red', lwd=3)
  
