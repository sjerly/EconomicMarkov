library(dplyr)
library(data.table)
library(deSolve)
library(ggplot2)
library(kableExtra)
library(facetscales)
library(gridExtra)


##############Validation Model######################

#Monthly Transition Matrix
moy<-read.csv("./MOY4.csv")
moymk<-as.matrix(moy[,-1])
moym<-moymk

#End of Year Transition Matrix
eoy<-read.csv("./EOY4.csv")
eoymk<-as.matrix(eoy[,-1])
eoym<-eoymk

#Population Growth
added<-read.csv("./Added12.csv")
addedmk<-as.matrix(added[1,])
addedm<-addedmk

#Population Counts for Variance
counts<-read.csv("./Year Counts.csv")
countsm<-as.matrix(counts[,2:3])

#Getting Row Names from Template
Val=read.csv("./Validation Values.csv")
init<-as.matrix((Val[1,-1]))
init<-t(as.matrix(c(439, 450, 0, 519, 377,  248, 0, 275,
                    63,  146, 0, 75,  1701, 122, 0, 484,
                    1447, 35, 0, 9,   91,   28, 0, 30),nrow=1))




out6<-data.frame(A6=numeric(),VS6=numeric(),overhead6=numeric(), recertcost6=numeric(), cmrecertcost6=numeric(),  PCost6=numeric(), Totcost6=numeric(), drop6=numeric(), inelig6=numeric(), Cpc6=numeric(),
                 A12=numeric(),VS12=numeric(), overhead12=numeric(), recertcost12=numeric(), cmrecertcost12=numeric(), PCost12=numeric(), Totcost12=numeric(), drop12=numeric(), inelig12=numeric(), Cpc12=numeric(), N12=numeric(),
                 RelEnroll=numeric(), RelVs=numeric(), RelCost=numeric(), PctInelg=numeric(), relrecertcost=numeric(), relcmrecertcost=numeric(),
                 relTotcost=numeric(), reldrop6=numeric(), relinelig6=numeric(), RelCpc6=numeric())


for (q in 1:10000) { 
  
  #Adding Uncertainty
  Va<-rbeta(1,0.90*3006+1,0.10*3006+1)
  Vd<-rbeta(1,0.69*1336+1,0.31*1336+1)
  Vi<-rbeta(1,0.81*896+1,0.19*896+1)
  
  co<-rnorm(1,40160.5,40160.5/12/1.96)
  cr0<-rnorm(1,98,98/12/1.96)
  cr1<-rnorm(1,202,202/12/1.96)
  calt<-rnorm(1,2188.8,2188.8/12/1.96)
  
  c1<-rgamma(1,shape=353^2/692,rate=353/692)
  c2<-rgamma(1,shape=414^2/736,rate=414/736)  
  c3<-rgamma(1,shape=2749^2/1483,rate=2749/1483)  
  c4<-rgamma(1,shape=413^2/717,rate=413/717)
  c5<-rgamma(1,shape=383^2/796,rate=383/796)  
  c6<-rgamma(1,shape=2789^2/1969,rate=2789/1969)  

 
  for (dc in 1:24){
    for (dr in 1:24) {
      moym[dr,dc]<-ifelse(moymk[dr,dc]==0,0,rbeta(1,moymk[dr,dc]*countsm[dr]+1,(1-moymk[dr,dc])*countsm[dr]+1,ncp=0))
      eoym[dr,dc]<-ifelse(eoymk[dr,dc]==0,0,rbeta(1,eoymk[dr,dc]*countsm[dr]+1,(1-eoymk[dr,dc])*countsm[dr]+1,ncp=0))
    }
      addedm[1,dc]<-rpois(1, addedmk[1,dc])
    
    }
diag(moym)<-1-(rowSums(moym)-diag(moym))
diag(eoym)<-1-(rowSums(eoym)-diag(eoym))
  
  
  #Initialize Values
  upd<-init
  store6<-init
  drop6<-0
  inelig6<-0
  
  for (i in 1:59) {
    #Number of Drops  
    #Number of Drops  
    drop6<-drop6+ifelse(i %in% 1:6,0,
                        ifelse(i %in% c(13,25), sum(t(as.matrix(upd[c(1,5,9,13,17,21)]))%*%eoym[c(1,5,9,13,17,21),c(2,6,10,14,18,22)]),
                               sum(t(as.matrix(upd[c(1,5,9,13,17,21)]))%*%moym[c(1,5,9,13,17,21),c(2,6,10,14,18,22)])))
    
    #Number Ruled Ineligible 
    inelig6<-inelig6+ifelse(i %in% c(13,25), sum(t(as.matrix(upd[c(1,5,9,13,17,21)]))%*%eoym[c(1,5,9,13,17,21),c(4,8,12,16,20,24)]),
                                   sum(t(as.matrix(upd[c(1,5,9,13,17,21)]))%*%moym[c(1,5,9,13,17,21),c(4,8,12,16,20,24)]))
    
    ifelse(i %in% c(13,25),upd<-upd%*%eoym,upd<-upd%*%moym)  
    upd<-upd+addedm
    store6<-rbind(store6,upd) 
  }
  
  A6<-rowSums(store6[,c(1,5,9,13,17,21,3,7,11,15,19,23)])
  D6<-rowSums(store6[,c(2,6,10,14,18,22)])
  I6<-rowSums(store6[,c(4,8,12,19,20,24)])
  





#12 Month Recert Model
###Translating Matrix to 12 Month

moym12<-moym

  #Move People Move to N at the rate they formerly went to I
  moym12[,c(3,7,11,15,19,23)]<-moym12[,c(4,8,12,16,20,24)]
  
  #Drop at Half the Rate 
    for (r in 1:24) {
      for (c in 1:24) {
        if(c %in% c(2,6,10,14,18,22) & r %in% c(1,5,9,13,17,21)) {
          moym12[r,c]<-moym12[r,c]/2}}}
  
  #Those People Stay in ADAP
   
    for (r in 1:24) {
      for (c in 1:24) {
         if(c %in% c(1,5,9,13,17,21) & r %in% c(1,5,9,13,17,21)) {
            moym12[r,c]<-moym12[r,c]+moym12[r,c+1]}}}
  
  
  #People only get to inelig through N    
  
    for (r in 1:24) {
      for (c in 1:24) {
        if(c %in% c(4,8,12,16,20,24) & r != c) {
           moym12[r,c]<-0}}}
  
  #People Return to ADAP the same rate from N as they did from I
  
    for (r in 1:24) {
      for (c in 1:24) {
        if(c %in% c(1,5,9,13,17,21) & r %in% c(3,7,11,15,19,23)) {
           moym12[r,c]<-moym12[r+1,c]}}}
          
  
  #Now 1/6 of people who are in N go to I every month      
        
    for (r in 1:24) {
      for (c in 1:24) {
        if(c %in% c(4,8,12,16,20,24) & r %in% c(3,7,11,15,19,23)) {
            moym12[r,c]<-moym12[r+1,c-1]*1/6
            moym12[r,c-1]<-moym12[r+1,c-1]*5/6
            moym12[r+1,c-1]<-0}}}
  
  
  #Resetting Diagonal so that populations stay fixed
  
      for (r in 1:24) {
        for (c in 1:24) {
          if(c %in% c(3,7,11,15,19,23) & r==c) {
              moym12[r,c+1]<-0
              moym12[r,c]<-1-(sum(moym12[r,])-moym12[r,c])
              moym12[r,c+1]<-moym12[r,c]*1/6
              moym12[r,c]<- moym12[r,c]*5/6}}}
              diag(moym12)<-1-(rowSums(moym12)-diag(moym12))

###Same thing for EOY
  
  eoym12<-eoym
  
  #Move People Move to N at the rate they formerly went to I
  eoym12[,c(3,7,11,15,19,23)]<-eoym12[,c(4,8,12,16,20,24)]
  
  #Drop at Half the Rate 
  for (r in 1:24) {
    for (c in 1:24) {
      if(c %in% c(2,6,10,14,18,22) & r %in% c(1,5,9,13,17,21)) {
        eoym12[r,c]<-eoym12[r,c]/2}}}
  
  #Those People Stay in ADAP
  
  for (r in 1:24) {
    for (c in 1:24) {
      if(c %in% c(1,5,9,13,17,21) & r %in% c(1,5,9,13,17,21)) {
        eoym12[r,c]<-eoym12[r,c]+eoym12[r,c+1]}}}
  
  
  #People only get to inelig through N    
  
  for (r in 1:24) {
    for (c in 1:24) {
      if(c %in% c(4,8,12,16,20,24) & r != c) {
        eoym12[r,c]<-0}}}
  
  #People Return to ADAP the same rate from N as they did from I
  
  for (r in 1:24) {
    for (c in 1:24) {
      if(c %in% c(1,5,9,13,17,21) & r %in% c(3,7,11,15,19,23)) {
        eoym12[r,c]<-eoym12[r+1,c]}}}
  
  
  #Now 1/6 of people who are in N go to I every month      
  
  for (r in 1:24) {
    for (c in 1:24) {
      if(c %in% c(4,8,12,16,20,24) & r %in% c(3,7,11,15,19,23)) {
        eoym12[r,c]<-eoym12[r+1,c-1]*1/6
        eoym12[r,c-1]<-eoym12[r+1,c-1]*5/6
        eoym12[r+1,c-1]<-0}}}
  
  
  #Resetting Diagonal so that populations stay fixed
  
  for (r in 1:24) {
    for (c in 1:24) {
      if(c %in% c(3,7,11,15,19,23) & r==c) {
        eoym12[r,c+1]<-0
        eoym12[r,c]<-1-(sum(eoym12[r,])-eoym12[r,c])
        eoym12[r,c+1]<-eoym12[r,c]*1/6
        eoym12[r,c]<- eoym12[r,c]*5/6}}}  
        diag(eoym12)<-1-(rowSums(eoym12)-diag(eoym12))

  

####For First 6 Months nobody is removed from the program 
  
  f6mm12<-moym12

  for (r in 1:24) {
    for (c in 1:24) {
      if(c %in% c(2,6,10,14,18,22,
                  4,8,12,16,20,24) & r!=c) {
          f6mm12[r,c]<-0}}}  
          diag(f6mm12)<-1-(rowSums(f6mm12)-diag(f6mm12))   
  


 #Initialize Values
  upd<-init
  store12<-init
  drop12<-0
  inelig12<-0

for (i in 1:59) {
  
  #Number of Drops  
  drop12<-drop12+ifelse(i %in% 1:6,0,
                        ifelse(i %in% c(13,25), sum(t(as.matrix(upd[c(1,5,9,13,17,21)]))%*%eoym12[c(1,5,9,13,17,21),c(2,6,10,14,18,22)]),
                               sum(t(as.matrix(upd[c(1,5,9,13,17,21)]))%*%moym12[c(1,5,9,13,17,21),c(2,6,10,14,18,22)])))
  
  #Number Ruled Ineligible 
  inelig12<-inelig12+ifelse(i %in% 1:6,0,
                            ifelse(i %in% c(13,25), sum(t(as.matrix(upd[c(3,7,11,15,19,23)]))%*%eoym12[c(3,7,11,15,19,23),c(4,8,12,16,20,24)]),
                                   sum(t(as.matrix(upd[c(3,7,11,15,19,23)]))%*%moym12[c(3,7,11,15,19,23),c(4,8,12,16,20,24)])))
  
  
  ifelse(i %in% 1:6,upd<-upd%*%f6mm12,ifelse(i %in% c(13,25),upd<-upd%*%eoym12,upd<-upd%*%moym12))  
  upd<-upd+addedm
  store12<-rbind(store12,upd) 
}

A12<-rowSums(store12[,c(1,5,9,13,17,21,3,7,11,15,19,23)])
D12<-rowSums(store12[,c(2,6,10,14,18,22)])
I12<-rowSums(store12[,c(4,8,12,19,20,24)])
N12<-rowSums(store12[,c(3,7,11,15,19,23)])  




Metrics6<-as.data.frame(cbind(A6,I6,D6))%>%rename(A=A6,I=I6,D=D6)%>%
  mutate(VS=(A*Va+D*Vd+(I)*Vi), 
         overhead=co,
         servcost=c1*store6[,c(1,3)]+c2*store6[,c(5,7)]+c3*store6[,c(9,11)]+c4*store6[,c(13,15)]+
           c5*store6[,c(17,19)]+c6*store6[,c(21,23)],
         recertcost=cr0*A/6,
         cmrecertcost=((cr1)*(store6[,c(13,15)]+store6[,c(17,19)]+store6[,c(21,23)]))/6,
         soccost=calt*I+D,
         Totcost=servcost+overhead+recertcost+cmrecertcost+soccost,
         PCost=servcost+overhead+recertcost+cmrecertcost)

A6a<-sum(Metrics6$A)/60
VS6<-Metrics6[60,"VS"]
overhead6<-sum(Metrics6$overhead)/5
recertcost6<-sum(Metrics6$recertcost)/5
cmrecertcost6<-sum(Metrics6$cmrecertcost)/5
Totcost6<-sum(Metrics6$Totcost)/5
PCost6<-sum(Metrics6$PCost)/5
Cpc6<-PCost6/A6a*12

Metrics12<-as.data.frame(cbind(A12,I12,D12))%>%rename(A=A12,I=I12,D=D12)%>%
  mutate(VS=(A*Va+D*Vd+(I)*Vi), 
         overhead=co,
         servcost=c1*store12[,c(1,3)]+c2*store12[,c(5,7)]+c3*store12[,c(9,11)]+c4*store12[,c(13,15)]+
           c5*store12[,c(17,19)]+c6*store12[,c(21,23)],
         recertcost=cr0*A/12,
         cmrecertcost=((cr1)*(store12[,c(13,15)]+store12[,c(17,19)]+store12[,c(21,23)]))/12,
         soccost=calt*I+D,
         Totcost=servcost+overhead+recertcost+cmrecertcost+soccost,
         PCost=servcost+overhead+recertcost+cmrecertcost)


A12a<-sum(Metrics12$A)/60
N12a<-sum(N12)
VS12<-Metrics12[60,"VS"]
overhead12<-sum(Metrics6$overhead)/5
recertcost12<-sum(Metrics12$recertcost)/5
cmrecertcost12<-sum(Metrics12$cmrecertcost)/5
Totcost12<-sum(Metrics12$Totcost)/5
PCost12<-sum(Metrics12$PCost)/5
Cpc12<-PCost12/A12a*12

out6[q,]<-data.frame(A6a, VS6, overhead6, recertcost6, cmrecertcost6, PCost6, Totcost6, drop6, inelig6, Cpc6,
                     A12a,VS12,overhead12,recertcost12,cmrecertcost12,PCost12,Totcost12,drop12,inelig12,Cpc12, N12a,
                     A12a/A6a, VS12/VS6, PCost12/PCost6, N12a/A12a/60,
                     recertcost6/recertcost12, cmrecertcost6/cmrecertcost12, Totcost6/Totcost12, drop6/drop12, 
                     inelig6/inelig12, Cpc6/Cpc12)

print(q)
}

x<-rbind(sapply(out6, function(x) quantile(x,0.975)), sapply(out6, function(x) quantile(x,0.025)))


out6$budgetimp<-out6$PCost6-out6$PCost12
out6$VSImp<-out6$VS12-out6$VS6

hist(out6$budgetimp)
hist(out6$VSImp)
hist(out6$A12/60-out6$A6/60)








