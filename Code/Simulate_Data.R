
 
#Starting perpendicular distance to trackline
perp.dist<-runif(N,0,W_Below)
tTime<-60 #Transect time (total time of obsrvation process for each whale)
speed<-0.309 #Speed of ship in kilometers per minute
 

#Dive Parameters (From Watwood et al. 2006)
DR<-1.2*60/1000 #Descent rate (convert to kilometers per minute)
         
#Blind spot
B.angle=0
B.rad=B.angle/57.2958
      
Final_Dat<-matrix(NA,nrow=1,ncol=27)

 for (i in 1:N.true){
      
#perpendicular distance whale i
cDist<-perp.dist[i]
      
Dive.Time<-rnorm(1,55,5) #Randomly assign a total time for a complete dive cycle
   
p1<-0.17  #17 percent of time at surface
p2<-0.06  #6 percent of time in non-clicking descent stae
p3<-0.11 #11 percent of time in clicking descent stae
p4<-0.50 #x percent of time in clicking foraging stae
p5<-0.16    #16 percent of time in non-clicking ascent state
   
   
Surface.time<-p1*Dive.Time
Descent.time<-(p2+p3)*Dive.Time
Ascent.time<-p5*Dive.Time
         
maxTime<-Dive.Time  #Total time of dive cycle
sTime<-runif(1,1,maxTime) #Time of current dive cycle (where is the animal in it's current dive cycle)
            
#Cumulative Amount of time in each dive state
Time1<-maxTime*p1
Time2<-maxTime*p2+maxTime*p1
Time3<-maxTime*p3+maxTime*p2+maxTime*p1
Time4<-maxTime*p4+maxTime*p3+maxTime*p2+maxTime*p1
Time5<-maxTime*p5+maxTime*p4+maxTime*p3+maxTime*p2+maxTime*p1
Time6<-maxTime

#Simulate dive process for whale i
cMat<-matrix(NA, nrow=tTime, ncol=24)

for (j in 1:(tTime)) {
                           
#Current Dive State
cTime<-((sTime+(j-1))%%maxTime)#current time mod maxTime

        
#Assign a dive state based on where the whale is in its dive cycle  
State1<-as.numeric(cTime<=Time1)#Surface
State2<-as.numeric(cTime>Time1)*as.numeric(cTime<=Time2)  #Non-clicking dive
State3<-as.numeric(cTime>Time2)*as.numeric(cTime<=Time3) #Clicking dive
State4<-as.numeric(cTime>Time3)*as.numeric(cTime<=Time4)  #Clicking and foraging 
State5<-as.numeric(cTime>Time4)*as.numeric(cTime<=Time5)
State6<-as.numeric(cTime>Time5)*as.numeric(cTime<=Time6)
           
     
#Indicator variables for dive state (i.e., above or below and clicking or non-clicking)
Above<-State1 #Surface state (availabe to visual team)
Clicking<-(State3+State4) #Clicking state (availabe to visual team)
Below<-(State2+State3+State4+State5)
S_I<-Above 
D_I<-(State2+State3) #Diving state
B_I<-State4 #Bottom foraging state
A_I<-State5 #Surfacing state

#Calculate maximum diving depth and ascent rate               
maxDepth<-Descent.time*DR
AR<-maxDepth/Ascent.time #Ascent rate
        
#Calculate current depth and forward distance at time j
cDepth<-(maxDepth*B_I+(cTime-Surface.time)*DR*D_I+(maxDepth-(cTime-(maxTime-Ascent.time))*AR)*A_I)
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
cMat[j,14]<-maxTime
cMat[j,15]<-maxDepth
cMat[j,16]<-S_I
cMat[j,17]<-D_I
cMat[j,18]<-B_I
cMat[j,19]<-A_I                       
cMat[j,20]<-Time1
cMat[j,21]<-j
cMat[j,22]<-blindspot
cMat[j,23]<-(State1*1)+(State2*2)+(State3*3)+(State4*4)+(State5*5)+(State6*6)
cMat[j,24]<-(Above+Clicking)*(1-blindspot)  #In an observable state
                   
              }
              
#Assign a NewID if whale i surfaces and then dives again and is detected (becomes a new click train)
transition<-array(NA,tTime)#Keep track of transitions 
transition[1]<-0
transition[2:(tTime)]<-cMat[1:(tTime)-1,11]-cMat[2:(tTime),11] 

Below_vec<-cMat[,12]
Ever_Below<-array(NA,tTime)                                    
Ever_Below[1]<-Below_vec[1]
newID<-array(NA,tTime)
newID[1]<-i
             
  for (j in 2:(tTime)) {
    Ever_Below[j]<-Ever_Below[j-1]+as.numeric(transition[j]<0)
    newID[j]<-i+N.true*as.numeric(Ever_Below[j]>1)#All new IDs > N.true (true number of whales)
                       }

cMat_All<-cbind(cMat,transition,Ever_Below,newID)
Final_Dat<-rbind(Final_Dat,cMat_All)
                      
                      }
                      
                      
n.tot.rows<-dim(Final_Dat)[1] #Get rid of extra row   
Final_Dat<-Final_Dat[2:n.tot.rows,]
                   
Final_Dat<-as.data.frame(Final_Dat)
colnames(Final_Dat)<-c("EventID", "X.True", "Y.True", "Depth","RadialDist", "X", "Y", "BearingR", "BearingD", 
                 "Above.Indicator", "Click.Indicator", "Below.Indicator", "obs.Below",  
                 "Dive.Cycle.Time", "maxDepth",  "Surface.Indicator", "Diving.Indicator", "Bottom.Indicator", 
                 "Ascending.Indicator", "Surf.Time", "Time", "Blind", "State", "Observable", "transition", "Ever_Below", "NewEventID") 
                 
                 
#Export simuated data                      
#write.csv(Final_Dat, "Sim.Real.Data.csv") 
             
#Subset observations by towed array (i.e., data for analysis)
Below_Dat<-subset(Final_Dat, obs.Below==1) 
#write.csv(Below_Dat, "Simulated_Below_for_analysis.csv")
                 
                 
###########################################  
#Simulate detection from half-normal
##########################################

#Visual (Above))
Above_Dat_All<-subset(Final_Dat, Above.Indicator==1 & Y.True>0)  
Above_Dat<-Above_Dat_All%>%dplyr::group_by(EventID)%>%summarise(distance = first(X.True))
N.tot.Above<-dim(Above_Dat)[1]  


#Half-Normal  detection
distance.A<-Above_Dat$distance
p.det.A<-exp(-distance.A^2/(2*sigma.A))
obs.A<-as.numeric(runif(N.tot.Above)<p.det.A)
simData.A<-data.frame(distance=distance.A,obs=obs.A)
Above.obs<-subset(simData.A, obs==1)
N_Above<-dim(Above.obs)[1]
p.above<-N_Above/N.tot.Above
  
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
N_Below<-dim(Below.obs)[1]
p.below<-N_Below/N.tot.Below
                
         
#Calculate true availability bias for visual team
Availability<-N.tot.Above/N.true #True availability    
         
       
        
 
             
       