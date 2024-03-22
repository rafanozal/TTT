
K<-1:5; N<-1000
models<-c(1,2,3)
param<-expand.grid(K,models)
# 
# # 

M<-200
powerH1<-matrix(NA,1,2)
for (num in 1:5)
{
tarea<-param[num,1]
model<-param[num,2]

fil.res<-paste('Power_H1_',model,'_tarea',tarea,'_N_',N,'.txt',sep="")
res<-as.matrix(read.table(fil.res,h=T),M,2)
powerH1<-rbind(powerH1,res)
}
colnames(powerH1)<-c( "i" ,"test.power")
powerH1<-powerH1[-1,]
#####
powerH2<-matrix(NA,1,2)
for (num in 6:10)
{
  tarea<-param[num,1]
  model<-param[num,2]
  
  fil.res<-paste('Power_H1_',model,'_tarea',tarea,'_N_',N,'.txt',sep="")
  res<-as.matrix(read.table(fil.res,h=T),M,2)
  powerH2<-rbind(powerH2,res)
}
colnames(powerH2)<-c( "i" ,"test.power")
powerH2<-powerH2[-1,]

#####
powerH3<-matrix(NA,1,2)
for (num in 11:15)
{
  tarea<-param[num,1]
  model<-param[num,2]
  
  fil.res<-paste('Power_H1_',model,'_tarea',tarea,'_N_',N,'.txt',sep="")
  res<-as.matrix(read.table(fil.res,h=T),M,2)
  powerH3<-rbind(powerH3,res)
}
colnames(powerH3)<-c( "i" ,"test.power")
powerH3<-powerH3[-1,]


write.table(powerH1,'Power_H1_IFR_N_1000.txt')
write.table(powerH2,'Power_H1_BFR_N_1000.txt')
write.table(powerH3,'Power_H1_UFR_N_1000.txt')