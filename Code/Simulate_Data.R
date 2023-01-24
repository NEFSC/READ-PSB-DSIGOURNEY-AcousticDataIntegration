#New version for revised manuscript (for PeerJ) (Original verion on Github)

simulate_data<-function(N.all,W_Below,W_Above,b0_b,b.dist_b,
                   sigma.A,sigma.B){ 
                   
#Starting perpendicular distance to trackline
perp.dist<-runif(N.all,0,W_Below)
tTime<-60 #Transect time (total time of obsrvation process for each whale)
speed<-0.25

#Dive Parameters (From Watwood et al. 2006)
DR<-1.2*60/1000 #Descent rate (convert to kilometers per minute)
     
#Blind spot
B.angle=0
B.rad=B.angle/57.2958
      
Final_Dat<-matrix(NA,nrow=1,ncol=30)


 for (i in 1:N.all){
      
#perpendicular distance whale i
cDist<-perp.dist[i]
      
#Assign dive parameters 
Dive.Time<-rnorm(1,55,2) #Total dive cycle (including time at surface)

   p1<-0.17  #17 percent of time at surface
   p2<-0.06  #6 percent of time in non-clicking descent stae
   p3<-0.11  #11 percent of time in clicking descent state
   p4<-0.50  #50 percent of time in clicking foraging state
   p5<-0.16  #16 percent of time in clicking ascent state
   
Surface.time<-p1*Dive.Time
Descent.time<-(p2+p3)*Dive.Time
Ascent.time<-p5*Dive.Time

sTime<-runif(1,0,Dive.Time) #Time of current dive cycle (where is the whale in its current dive cycle)
                      
#Cumulative Amount of time in each dive state
  Time1<-Dive.Time*p1
  Time2<-Dive.Time*p2+Dive.Time*p1
  Time3<-Dive.Time*p3+Dive.Time*p2+Dive.Time*p1
  Time4<-Dive.Time*p4+Dive.Time*p3+Dive.Time*p2+Dive.Time*p1
  Time5<-Dive.Time*p5+Dive.Time*p4+Dive.Time*p3+Dive.Time*p2+Dive.Time*p1
  Time6<-Dive.Time

#Simulate dive process for whale i
cMat<-matrix(NA, nrow=tTime, ncol=26)

for (j in 1:(tTime)) {
                           
#Current Dive State
cTime<-((sTime+(j-1))%%Dive.Time)#current time mod maxTime

#Assign a dive state based on where the whale is in its dive cycle        
State1<-as.numeric(cTime<=Time1)
State2<-as.numeric(cTime>Time1)*as.numeric(cTime<=Time2)
State3<-as.numeric(cTime>Time2)*as.numeric(cTime<=Time3)
State4<-as.numeric(cTime>Time3)*as.numeric(cTime<=Time4)
State5<-as.numeric(cTime>Time4)*as.numeric(cTime<=Time5)
State6<-as.numeric(cTime>Time5)*as.numeric(cTime<=Time6)
             
#State1 is the surface state (i.e. available to visual team)

Above<-State1
SI<-Above
           
#Clicking state          

Clicking<-State3+State4 
Below<-State2+State3+State4+State5
DI<-State2+State3
BI<-State4 
AI<-State5
             
#Calculate maxDeopth bassed on descent time and descent rate
maxDepth<-Descent.time*DR
AR<-maxDepth/Ascent.time

       
             
#Calculate current depth and forward distance at time j
cDepth<-(maxDepth*BI+(cTime-Surface.time)*DR*DI+(maxDepth-(cTime-(Dive.Time-Ascent.time))*AR)*AI)
cY<-initial.y-speed*(j-1)           
cSR<-(cDist^2+cY^2+cDepth^2)^0.5 #Calculate current slant range
cAngle.R<-(atan(cDist/cY)*as.numeric(cY>0))+((atan(abs(cY)/cDist)+90*0.017453294)*as.numeric(cY<=0)) #Calculate current bearing in radians
cAngle.D<- cAngle.R/0.017453294 #Convert radians to degrees
cY.calc<-cos(cAngle.R)*cSR #Y calculated from observed bearing and current slant range (assumes 2D) 
cX.calc<-sin(cAngle.R)*cSR  #X calculated from observed bearing and current slant range (assumes 2D) 
                        
#Determine if whale is in a blindspot (i .e., cannot be detected by acoustic array) 
#X.blind<-cY*tan(B.rad)
#blindspot<-as.numeric(cDist<X.blind) #blindspot=1 if whale is in blindspot 
blindspot<-0  #turn blindspot off for now

#Simulate Observation Proces            
Obs.below<-as.numeric(runif(1)<plogis(b0_b+b.dist_b*cSR))*Clicking*(1- blindspot)     
  
#Record simulated data  
cMat[j,1]<-i
cMat[j,2]<-cDist 
cMat[j,3]<-cY 
cMat[j,4]<-cDepth 
cMat[j,5]<-cSR
cMat[j,6]<-cX.calc
cMat[j,7]<-cY.calc
cMat[j,8]<-cAngle.R
cMat[j,9]<-cAngle.D
cMat[j,10]<-Above
cMat[j,11]<-Clicking
cMat[j,12]<-Below
cMat[j,13]<-Obs.below
cMat[j,14]<-cTime
cMat[j,15]<-Dive.Time
cMat[j,16]<-cTime/Dive.Time
cMat[j,17]<-maxDepth
cMat[j,18]<-SI
cMat[j,19]<-DI
cMat[j,20]<-BI
cMat[j,21]<-AI                       
cMat[j,22]<-Time1
cMat[j,23]<-j
cMat[j,24]<-blindspot
cMat[j,25]<-(State1*1)+(State2*2)+(State3*3)+(State4*4)+(State5*5)+(State6*6)
cMat[j,26]<-(Above+Clicking)*(1-blindspot)  #In an observable state
                   
              }
              
#Keep track of tranistions into (-1) and out of (1) the Diving state
trans.D<-array(NA,tTime)
trans.D[1]<-0
trans.D[2:(tTime)]<-cMat[2:(tTime),12] -cMat[1:(tTime)-1,12]

              
#Keep track of tranistions into (-1) and out of (1) the foraging (i.e., clicking) state
trans.F<-array(NA,tTime)
trans.F[1]<-0
trans.F[2:(tTime)]<-cMat[2:(tTime),11]-cMat[1:(tTime)-1,11] 


Below_vec<-cMat[,12]
Ever_Below<-array(NA,tTime)                                    
Ever_Below[1]<-Below_vec[1]
newID<-array(NA,tTime)
newID[1]<-i
             
  for (j in 2:(tTime)) {
    Ever_Below[j]<-Ever_Below[j-1]+as.numeric(trans.F[j]<0)
    newID[j]<-i+N.all*as.numeric(Ever_Below[j]>1)#All new IDs > N.all (true number of all whales available to acoutic platform)
                       }

cMat_All<-cbind(cMat,trans.D,trans.F,Ever_Below,newID)
Final_Dat<-rbind(Final_Dat,cMat_All)
                      
                      }
                                           
n.tot.rows<-dim(Final_Dat)[1] #Get rid of extra row   
Final_Dat<-Final_Dat[2:n.tot.rows,]
                   
Final_Dat<-as.data.frame(Final_Dat)
colnames(Final_Dat)<-c("EventID", "X.True", "Y.True", "Depth","RadialDist", "X", "Y", "BearingR", "BearingD", 
                 "Above.Indicator", "Click.Indicator", "Below.Indicator", "obs.Below", "Current.Dive.Time", 
                 "Dive.Cycle.Time", "Percent.Dive.Time",  "maxDepth",  "Surface.Indicator", "Diving.Indicator", "Bottom.Indicator", 
                 "Ascending.Indicator", "Surf.Time", "Time", "Blind", "State", "Observable", "trans.D","trans.F",
                  "Ever_Below", "NewEventID") 
                 
                 
#Export simuated data                      
#write.csv(Final_Dat, "Sim.Real.Data.csv") 
             
#Subset observations by towed array (i.e., data for analysis)
Below_Dat<-subset(Final_Dat, obs.Below==1 & Y<W_Below) 
#write.csv(Below_Dat, "Simulated_Below_for_analysis.csv")
                 
                 
###########################################  
#Simulate detection from half-normal
##########################################

#Visual (Above))
Above_Dat_All<-subset(Final_Dat, Above.Indicator==1 & X.True<W_Above & Y.True>0 & Y.True<W_Above)  
Above_Dat<-Above_Dat_All%>%dplyr::group_by(EventID)%>%summarise(distance = first(X.True))
N.tot.Above<-dim(Above_Dat)[1]  


#Half-Normal  detection
distance.A<-Above_Dat$distance
p.det.A<-exp(-distance.A^2/(2*sigma.A))
obs.A<-as.numeric(runif(N.tot.Above)<p.det.A)
simData.A<-data.frame(distance=distance.A,obs=obs.A)
Above.obs<-subset(simData.A, obs==1)
Nobs_Above<-dim(Above.obs)[1]
p.above<-Nobs_Above/N.tot.Above
  
#Acoustic (Below) (for testing hybrid model)  
Below_Dat_All<-subset(Final_Dat, Below.Indicator==1 & Y.True>0)  
Below_Dat2<-Below_Dat_All%>%dplyr::group_by(EventID)%>%summarise(distance = first(X.True))
N.tot.Below<-dim(Below_Dat2)[1]  
  
#Half-Normal  detection
distance.B<-Below_Dat2$distance
p.det.B<-exp(-distance.B^2/(2*sigma.B))
obs.B<-as.numeric(runif(N.tot.Below)<p.det.B)
simData.B<-data.frame(distance=distance.B,obs=obs.B)
Below.obs<-subset(simData.B, obs==1)
Nobs_Below<-dim(Below.obs)[1]
p.below<-Nobs_Below/N.tot.Below
                
 sims<-list(Final_Dat=Final_Dat,Below_Dat=Below_Dat,Above.obs=Above.obs,Nobs_Above=Nobs_Above,
              Below.obs=Below.obs,Nobs_Below=Nobs_Below)        
 return(sims)
     }    
       

