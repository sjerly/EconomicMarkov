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


Va<-0.90
Vd<-0.69
Vi<-0.81
co<-481926/12
cr0<-98
cr1<-202
calt<-26266/12
c1<-353
c2<-414
c3<-2749
c4<-413
c5<-383
c6<-2789

out6<-data.frame(A6=numeric(),VS6=numeric(),overhead6=numeric(), recertcost6=numeric(), cmrecertcost6=numeric(),  PCost6=numeric(), Totcost6=numeric(), drop6=numeric(), inelig6=numeric(), Cpc6=numeric(),
                 A12=numeric(),VS12=numeric(), overhead12=numeric(), recertcost12=numeric(), cmrecertcost12=numeric(), PCost12=numeric(), Totcost12=numeric(), drop12=numeric(), inelig12=numeric(), Cpc12=numeric(), N12=numeric(),
                 RelEnroll=numeric(), RelVs=numeric(), RelCost=numeric(), PctInelg=numeric(), relrecertcost=numeric(), relcmrecertcost=numeric(),
                 relTotcost=numeric(), reldrop6=numeric(), relinelig6=numeric(), RelCpc6=numeric())

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
    
    #Matrix Multiplication for Model
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
  
  #Matrix Multiplication for Model
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





