

#Format VLT data for MRDS analysis

N_total<-dim(VLT_data)[1]
             
Dat1<-matrix(NA,N_total,6)
Dat2<-matrix(NA,N_total,6)

Observ_t<-VLT_data$team   #1 for upper/front and 2 for lower/back
Detect_t<-VLT_data$detected
Dist_t<-VLT_data$distance  #Convert to kilometers
Beaufort_t<-VLT_data$beaufort
Size_t<-VLT_data$size
Y_t<-VLT_data$Y

for (i in 1:N_total) {

#Dat1
Dat1[i,1]=Observ_t[i]
Dat1[i,2]=Detect_t[i]
Dat1[i,3]=Dist_t[i]
Dat1[i,4]=Y_t[i]
Dat1[i,5]=Beaufort_t[i]
Dat1[i,6]=Size_t[i]

#Dat2
Dat2[i,1]=Observ_t[i]
Dat2[i,2]=Detect_t[i]
Dat2[i,3]=Dist_t[i]
Dat2[i,4]=Y_t[i]
Dat2[i,5]=Beaufort_t[i]
Dat2[i,6]=Size_t[i]

                        }

Dat1<-as.data.frame(Dat1)
colnames(Dat1)=c("Observer", "Detected", "Distance", "Y", "Beaufort", "Group.Size")

Dat1<-subset(Dat1,Observer==1)
Obs1<-Dat1$Detected

Dat2<-as.data.frame(Dat2)
colnames(Dat2)=c("Observer", "Detected", "Distance", "Y", "Beaufort", "Group.Size")

Dat2<-subset(Dat2,Observer==2)
Obs2<-as.numeric(Dat2$Detected)

N_sights<-length(Obs1)
y10<-array(NA,N_sights)
y01<-array(NA,N_sights)
y11<-array(NA,N_sights)

for (i in 1:N_sights){
y10[i]<-(1*(Obs1[i]>0))*(1*(Obs2[i]<1))
y01[i]<-(1*(Obs1[i]<1))*(1*(Obs2[i]>0))
y11[i]<-(1*(Obs1[i]>0))*(1*(Obs2[i]>0))
                     }

VLT_MRDS<-cbind(Dat1,y10,y01,y11)
