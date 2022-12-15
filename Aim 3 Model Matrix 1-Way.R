  #Adding Uncertainty
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

moym<-moymk
eoym<-eoymk
addedm<-addedmk

#for (dr in 1:24) {
#  for (dc in 1:24) {
#    if(dc %in% c(1,5,9,13,17,21) & dr %in% c(4,8,12,16,20,24)) {
#      moym[dr,dc]<-ifelse(moymk[dr,dc]==0,0,qbeta(0.975,moymk[dr,dc]*countsm[dr]+1,(1-moymk[dr,dc])*countsm[dr]+1,ncp=0))
#      eoym[dr,dc]<-ifelse(eoymk[dr,dc]==0,0,qbeta(0.975,eoymk[dr,dc]*countsm[dr]+1,(1-eoymk[dr,dc])*countsm[dr]+1,ncp=0))
#      }}}
  

#for (dc in 1:24) {
#addedm[1,dc]<-qpois(0.025, addedmk[1,dc])
#}

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
  
  
  #Now 1/12 of people who are in N go to I every month      
  
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
  
  
  #Now 1/12 of people who are in N go to I every month      
  
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
  
  
  
  Metrics6<-as.data.frame(cbind(A6,I6,D6))%>%dplyr::rename(A=A6,I=I6,D=D6)%>%
    mutate(VS=(A*Va+D*Vd+(I)*Vi), 
           overhead=co,
           servcost=c1*store6[,c(1,3)]+c2*store6[,c(5,7)]+c3*store6[,c(9,11)]+c4*store6[,c(13,15)]+
             c5*store6[,c(17,19)]+c6*store6[,c(21,23)],
           recertcost=cr0*A/6,
           cmrecertcost=((cr1)*(store6[,c(13,15)]+store6[,c(17,19)]+store6[,c(21,23)]))/6,
           soccost=calt*I+D,
           Totcost=servcost+overhead+recertcost+cmrecertcost+soccost,
           PCost=servcost+overhead+recertcost+cmrecertcost)

  Metrics12<-as.data.frame(cbind(A12,I12,D12,N12))%>%dplyr::rename(A=A12,I=I12,D=D12,N=N12)%>%
    mutate(VS=(A*Va+D*Vd+(I)*Vi), 
           overhead=co,
           servcost=c1*store12[,c(1,3)]+c2*store12[,c(5,7)]+c3*store12[,c(9,11)]+c4*store12[,c(13,15)]+
             c5*store12[,c(17,19)]+c6*store12[,c(21,23)],
           recertcost=cr0*A/12,
           cmrecertcost=((cr1)*(store12[,c(13,15)]+store12[,c(17,19)]+store12[,c(21,23)]))/12,
           soccost=calt*I+D,
           Totcost=servcost+overhead+recertcost+cmrecertcost+soccost,
           PCost=servcost+overhead+recertcost+cmrecertcost)
  

 
 
 
 ##### Graph #############
 
 library(ggplot2)
 library(plyr)
 library(dplyr)
 library(tidyverse)
 
 Parameter<-c(
   "Recertification Staff Cost",
   "Cost of CM Assistance",
   "Service Cost: Uninsured",
   "Service Cost: Public Insurance",
   "Service Cost: Private Insurance",
   "Disenrollment Rate",
   "Ineligibility Rate",
   "Rate of Client Return from Disenrollment",
   "Rate of Client Return from Ineligibility",
   "Number of New ADAP Clients Per Month"
 )
 
 options(scipen=999)
 # this is throwing some warnings in my computer, but it is reading the data frame correctly
 df <- '
Lower_Bound	Upper_Bound	UL_Difference
2486615	2621310	134695
2498345	2609580	111235
2513491	2595031	81540
2457155	2657570	200415
2466078	2648414	182336
2260019	2873382	613363
2350403	2805796	455393
2875801	2116520 759281
2674273 2471037	203236
2129110	3108062	978952





' %>% read_table2()
 df<-cbind(Parameter,df)
 colnames(df)[2]<-"Lower Bound"
 colnames(df)[3]<-"Upper Bound" 
 
 # original value of output
 base.value <- 2553962
 
 # get order of parameters according to size of intervals
 # (I use this to define the ordering of the factors which I then use to define the positions in the plot)
 order.parameters <- df %>% arrange(UL_Difference) %>%
   mutate(Parameter=factor(x=Parameter, levels=Parameter)) %>%
   select(Parameter) %>% unlist() %>% levels()
 
 # width of columns in plot (value between 0 and 1)
 width <- 0.95
 
 # get data frame in shape for ggplot and geom_rect
 df.2 <- df %>% 
   # gather columns Lower_Bound and Upper_Bound into a single column using gather
   gather(key='type', value='output.value', "Lower Bound":"Upper Bound") %>%
   # just reordering columns
   select(Parameter, type, output.value, UL_Difference) %>%
   # create the columns for geom_rect
   mutate(Parameter=factor(Parameter, levels=order.parameters),
          ymin=pmin(output.value, base.value),
          ymax=pmax(output.value, base.value),
          xmin=as.numeric(Parameter)-width/2,
          xmax=as.numeric(Parameter)+width/2)
 
 # create plot
 # (use scale_x_continuous to change labels in y axis to name of parameters)
 png(width = 960, height = 540)
 ggplot() + 
   geom_rect(data = df.2, 
             aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, fill=type)) +
   theme_bw() + 
   theme(axis.text=element_text(color="black"), axis.title.y=element_blank(), legend.position = 'bottom',
         legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   geom_hline(yintercept = base.value) +
   scale_x_continuous(breaks = c(1:length(order.parameters)), 
                      labels = order.parameters) +
   scale_y_continuous(labels=scales::dollar_format())+
   coord_flip()
 dev.off()
 
 
 
 
 
 
 #Last Statistics
 
 
 Metrics6<-as.data.frame(cbind(A6,I6,D6))%>%dplyr::rename(A=A6,I=I6,D=D6)%>%
   mutate(VS=(A*Va+D*Vd+(I)*Vi), 
          overhead=co,
          servcost=c1*store6[,c(1,3)]+c2*store6[,c(5,7)]+c3*store6[,c(9,11)]+c4*store6[,c(13,15)]+
            c5*store6[,c(17,19)]+c6*store6[,c(21,23)],
          recertcost=cr0*A/6,
          cmrecertcost=((cr1)*(store6[,c(13,15)]+store6[,c(17,19)]+store6[,c(21,23)]))/6,
          soccost=calt*I+D,
          Totcost=servcost+overhead+recertcost+cmrecertcost+soccost,
          PCost=servcost+overhead+recertcost+cmrecertcost,
          yr=ceiling(row_number()/12))%>%group_by(yr)%>%
   summarise(overhead=sum(overhead), servcost=sum(servcost), recertcost=sum(recertcost), 
             cmrecertcost=sum(cmrecertcost), Totcost=sum(Totcost), PCost=sum(PCost), VS=mean(VS),
             A=mean(A))%>%ungroup
 
 
 Metrics12<-as.data.frame(cbind(A12,I12,D12,N12))%>%dplyr::rename(A=A12,I=I12,D=D12,N=N12)%>%
   mutate(VS=(A*Va+D*Vd+(I)*Vi), 
          overhead=co,
          servcost=c1*store12[,c(1,3)]+c2*store12[,c(5,7)]+c3*store12[,c(9,11)]+c4*store12[,c(13,15)]+
            c5*store12[,c(17,19)]+c6*store12[,c(21,23)],
          recertcost=cr0*A/12,
          cmrecertcost=((cr1)*(store12[,c(13,15)]+store12[,c(17,19)]+store12[,c(21,23)]))/12,
          soccost=calt*I+D,
          Totcost=servcost+overhead+recertcost+cmrecertcost+soccost,
          PCost=servcost+overhead+recertcost+cmrecertcost,
          yr=ceiling(row_number()/12))%>%group_by(yr)%>%
   summarise(overhead=sum(overhead), servcost=sum(servcost), recertcost=sum(recertcost), 
             cmrecertcost=sum(cmrecertcost), Totcost=sum(Totcost), PCost=sum(PCost), VS=mean(VS),
             A=mean(A), N=mean(N))%>%ungroup
 
 
 
 
 
 
 
