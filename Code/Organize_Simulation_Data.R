   

#All_dat<-read.csv("Simulated_Below_for_analysis.csv")


# All_dat<-subset(All_dat, Y<(y_max*-1) & X<x_max_trunc & Y>(y_min*-1)) 
   
     
ID_vec<-unique(Below_Dat$NewEventID)
N.whales<-length(ID_vec)
             
summary_mat<-matrix(NA,nrow= N.whales,ncol=4)

for (i in 1:N.whales){
        
cWhale<-ID_vec[i]
cData<-subset(Below_Dat, NewEventID==cWhale)
clicks<-dim(cData)[1]
meanX<-mean(as.numeric(cData$X))           
minY<-min(as.numeric(cData$Y))
maxY<-max(as.numeric(cData$Y))
             
summary_mat[i,1]<-cWhale
summary_mat[i,2]<-meanX
summary_mat[i,3]<-minY
summary_mat[i,4]<-maxY     
             }
             
summary_mat_all<-as.data.frame(summary_mat)
colnames(summary_mat_all)<-c("ID",  "meanX",  "minY", "maxY")
         
#write.csv(summary_mat_all, "summary_mat_all_sims.csv")
     
  
summary_mat<-summary_mat_all
ID_vec2<-unique(summary_mat$ID)
N.whales<-length(ID_vec2)    
            
#Prepare data for model         

Data_mat_final<-matrix(NA,nrow=1,ncol=7)

for (c in 1:N.whales){  
cWhale<-ID_vec2[c]
cSummary.Dat<-subset(summary_mat,ID==cWhale)
x<-cSummary.Dat$meanX
cData<-subset(Below_Dat, NewEventID==cWhale)
N.rows<-dim(cData)[1]   
y<-as.numeric(cData$Y)
  
y_mat<-matrix(NA,nrow=N.rows,ncol=n_bins)
Bin<-matrix(NA,nrow=N.rows,ncol=1)
  
#Assign distance bins
for (i in 1:N.rows){
  for (j in 1:n_bins){
       y_mat[i,j]<-as.numeric(y[i]<=y_bins[j])
      }
      Bin[i,1]<-sum(y_mat[i,])
      }
      
y<-as.matrix(y)
y_final<-cbind(y,Bin)
y_final<-as.data.frame(y_final)
colnames(y_final)<-c("Y", "Bin")
Data_mat<-matrix(NA,nrow=n_bins, ncol=7)

for (j in 1:n_bins){
 
cMat<-subset(y_final, Bin==j)
midY<-(y_bins[j]+y_bins[j+1])/2
 
Data_mat[j,1]<-cWhale #Whale number
Data_mat[j,2]<-j #Bin number
Data_mat[j,3]<-x
Data_mat[j,4]<-midY
Data_mat[j,5]<-y_bins[j+1]
Data_mat[j,6]<-y_bins[j] 
Data_mat[j,7]<-as.numeric(dim(cMat)[1]>0)

      }
     
Data_mat_final<-rbind(Data_mat_final,Data_mat)

      }
      
Data_mat_final<-Data_mat_final[2:dim(Data_mat_final)[1],]
Data_mat_final<-as.data.frame(Data_mat_final)
colnames(Data_mat_final)<-c("Whale","Bin","Distance","Mid.Y","Y_min","Y_max","Detection.Indicator")
 # write.csv(Data_mat_final,"Raw_Data_formatted.csv")
     
#Make CH matrix  
wide_detected <- reshape(Data_mat_final, v.names = "Detection.Indicator", idvar = "Whale",
                 timevar = "Bin", direction = "wide")
  
ncols<-dim(wide_detected)[2]
detected_mat<-as.matrix(wide_detected[,6:(ncols)]) 
#write.csv(detected_mat,"detected.csv")

#Subsample Observed Data 
samp<-(1:N.whales)
subsample<-sample(samp,n.sub, replace=FALSE)
wide_detected.sub<-wide_detected[subsample,]
detected_mat.sub<-as.matrix(wide_detected.sub[,6:(ncols)]) 

#Get perpendicular distances for visual analysis
dist_Above<-Above.obs$distance #For Hybrid Model

#Get peerpendicaulr distnaces for acoustic analysis
dist_Below1<-wide_detected$Distance #For CMR Model 
dist_Below2<-Below.obs$distance #For Hybrid Model (Line transect analysis)
dist_Below3<-wide_detected.sub$Distance #For Hybrid Model (CMR analysis)


#For "Ones"  trick
ones.dist_A<-array(1,N_Above)
ones.dist_B<-array(1,N_Below)

#Data Augmentation CMR-DS Method
get.first<-function(x)min(which(x!=0))
get.last<-function(x)max(which(x!=0))
f<-apply(detected_mat,1,get.first)
g<-apply(detected_mat,1,get.last)
n<-nrow(detected_mat)
t<-ncol(detected_mat)
         
y<-as.matrix(detected_mat)
nz<-120
y<-rbind(y,matrix(rep(0,ncol(y)*nz),ncol=ncol(y),byrow=TRUE))
M<-nrow(y)
T<-ncol(y)

#Make z matrix        
z.mat<-matrix(1,nrow=dim(y)[1], ncol=dim(y)[2])


#Make age matrix  
age<-matrix(0,nrow=n,ncol=t) 
for (i in 1:n){
    for (j in f[i]:t){
     age[i,j]<-j-f[i]+1
            }
            }
age<-rbind(age,matrix(rep(NA,ncol(age)*nz),ncol=ncol(age),byrow=TRUE)) 
    

dist_Below1<-c(dist_Below1,rep(NA,nz))


#Data Augmentation Hybrid method
f.sub<-apply(detected_mat.sub,1,get.first)
g.sub<-apply(detected_mat.sub,1,get.last)
n.sub<-nrow(detected_mat.sub)
t.sub<-ncol(detected_mat.sub)

y.sub<-as.matrix(detected_mat.sub)
nz.sub<-50
y.sub<-rbind(y.sub,matrix(rep(0,ncol(y.sub)*nz.sub),ncol=ncol(y.sub),byrow=TRUE))
M.sub<-nrow(y.sub)
T.sub<-ncol(y.sub)

#Make z matrix        
z.mat.sub<-matrix(1,nrow=dim(y.sub)[1], ncol=dim(y.sub)[2])

#Make age matrix  
age.sub<-matrix(0,nrow=n.sub,ncol=t.sub) 
for (i in 1:n.sub){
    for (j in f.sub[i]:t.sub){
     age.sub[i,j]<-j-f.sub[i]+1
            }
            }
age.sub<-rbind(age.sub,matrix(rep(NA,ncol(age.sub)*nz.sub),ncol=ncol(age.sub),byrow=TRUE)) 
    

dist_Below3<-c(dist_Below3,rep(NA,nz.sub))



#MCMC parameters
nitt=3000
burnin=3000
chains=1
thin=3
        
 

