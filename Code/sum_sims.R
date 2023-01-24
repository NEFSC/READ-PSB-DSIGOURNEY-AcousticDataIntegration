

#########################
 # SUMMARIZE SIMULATION OUTPUT
###############################


##############
 #TRUNCATION
##############
 sum_sims<-function(Final_Dat,W_Above){
#Subset out data past the visual observation window (Anything < forward distance =0)
trunc_Dat<-subset(Final_Dat,X.True<W_Above & Y.True>0 & Y.True<W_Above)
#trunc_Dat<-subset(Final_Dat,Y.True>0)
ID_vec_All<-unique(trunc_Dat$EventID)
N.true<-length(ID_vec_All) #True number of whales in zone of overlap



N.true.Below<-trunc_Dat%>%select(Below.Indicator,EventID)%>%filter(Below.Indicator==1)%>%
         distinct()%>%summarise(sum(Below.Indicator))%>%as.numeric()

N.true.Above<-trunc_Dat%>%select(Above.Indicator,EventID)%>%filter(Above.Indicator==1)%>%
         distinct()%>%summarise(sum(Above.Indicator))%>%as.numeric()
          
          

Availability<-N.true.Above/N.true #True surface availability

                     
Double_Dives<-matrix(NA,nrow= N.true,ncol=3)
Obs_Mat<-matrix(NA,nrow= N.true,ncol=2)
Dups_Mat<-matrix(NA,nrow= N.true,ncol=4)

  for (i in 1:N.true){
    cID<-ID_vec_All[i]
    cDat<-subset(trunc_Dat, EventID==cID)
    Trans_vec.D<-as.numeric(cDat$transition<0)
    Trans_vec.S<-as.numeric(cDat$transition>0)
    N.dives<-sum(Trans_vec.D)
    N.surface<-sum(Trans_vec.S)
    C_I<-as.numeric(sum(cDat$Click.Indicator)==0)
    B1_I<-as.numeric(cDat$Click.Indicator[1]==1) #Indicator that animal was below and clicking in first interval
    DD_I<-B1_I*as.numeric(N.surface>0)*as.numeric(N.dives>0)+as.numeric(N.dives>2)
    Double_Dives[i,1]<-as.numeric(DD_I>0)
    Double_Dives[i,2]<-N.dives
    Double_Dives[i,3]<-C_I
                        
    obs<-cDat$Observable
    n.intervals<-length(obs)
    unobservable<-as.numeric(sum(obs)<1)
    n.unobs<-length(which(obs<1))
    percent.unobs<-n.unobs/n.intervals
    Obs_Mat[i,1]<-unobservable
    Obs_Mat[i,2]<-percent.unobs
                      
    Above_I<-cDat$Above.Indicator
    Below_I<-cDat$Below.Indicator
    Click_I<-cDat$Click.Indicator
    A<-as.numeric(sum(Above_I)>0)
    B<-as.numeric(sum(Below_I)>0)
    C<-as.numeric(sum(Click_I)>0)
                      
    Dups_Mat[i,1]<- cID
    Dups_Mat[i,2]<-as.numeric((A+B)>1) #Indiviudal was both above and in some below state
    Dups_Mat[i,3]<-as.numeric((A+C)>1) #Indiviudal was both above and in a clicking state
    Dups_Mat[i,4]<-as.numeric(A<1)*as.numeric(sum(obs)<n.intervals) #Individuals both observable and unobservable but never Above
                      }
                      
    DD<-sum(Double_Dives[,1])
    unobservable<-sum(Obs_Mat[,1])
    max.time.obs<-max(Obs_Mat[,2])
    min.time.obs<-min(Obs_Mat[,2])
    med.time.obs<-median(Obs_Mat[,2])

    Dups_all<-sum(Dups_Mat[,2])
    Dups_click<-sum(Dups_Mat[,3])
    Dups_fake<-sum(Dups_Mat[,4])
    
    n_bins.trunc<-Bin.end-Bin.start+1 #Number of bins in the truncated zone of overlap                                         
    sum_mat<-matrix(NA,nrow= n_bins.trunc,ncol=7)
  for (i in 1:(n_bins.trunc)){
    j<-(Bin.start-1)+i
    y1<-y_bins[j]
    y2<-y_bins[j+1]
    cDat1<-subset(trunc_Dat,Y.True>y2 & Y.True<y1)
    cDat2<-subset(cDat1,Below.Indicator==1)
    cDat3<-subset(cDat1,Click.Indicator==1)
    cDat4<-subset(cDat1,Above.Indicator==1)
    cDat5<-subset(cDat1,trans.F==1)
    cDat6<-subset(cDat1,trans.F==-1)
                          
    N_1<-length(unique(cDat1$EventID))
    N_2<-length(unique(cDat2$EventID))
    N_3<-length(unique(cDat3$EventID))
    N_4<-length(unique(cDat4$EventID))
    N_5<-length(unique(cDat5$EventID))
    N_6<-length(unique(cDat6$EventID))
                            
    sum_mat[i,1]<-N_1
    sum_mat[i,2]<-N_2
    sum_mat[i,3]<-N_3
    sum_mat[i,4]<-N_4
    sum_mat[i,5]<-N_5
    sum_mat[i,6]<-N_6
    sum_mat[i,7]<-N_2/N_1
                         }
                         
    output<-list(N.true.Below=N.true.Below,N.true.Above=N.true.Above,
    Availability=Availability,Dups_click=Dups_click,sum_mat=sum_mat)
    
    return(output)                         
            }                 
                             
                             
                             