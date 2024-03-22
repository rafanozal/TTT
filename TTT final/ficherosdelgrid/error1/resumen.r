
K<-1:5;
Ns<-c(50,100,500,1000)
param<-expand.grid(K,Ns)
# 
# # 

M<-200
error50<-matrix(NA,1,3)
for (num in 1:5)
{
tarea<-param[num,1]
N<-param[num,2]

fil.res<-paste('Error_TypeI_N_',N,'_tarea_',tarea,'.txt',sep="")
res<-as.matrix(read.table(fil.res,h=T),M,3)
error50<-rbind(error50,res)
}
colnames(error50)<-c( "i" ,"SiZer-0","SiZer-2")
error50<-error50[-1,]
#####
error100<-matrix(NA,1,3)
for (num in 6:10)
{
  tarea<-param[num,1]
  N<-param[num,2]
  
  fil.res<-paste('Error_TypeI_N_',N,'_tarea_',tarea,'.txt',sep="")
  res<-as.matrix(read.table(fil.res,h=T),M,3)
  error100<-rbind(error100,res)
}
colnames(error100)<-c( "i" ,"SiZer-0","SiZer-2")
error100<-error100[-1,]

#####
error500<-matrix(NA,1,3)
for (num in 11:15)
{
  tarea<-param[num,1]
  N<-param[num,2]
  
  fil.res<-paste('Error_TypeI_N_',N,'_tarea_',tarea,'.txt',sep="")
  res<-as.matrix(read.table(fil.res,h=T),M,3)
  error500<-rbind(error500,res)
}
colnames(error500)<-c( "i" ,"SiZer-0","SiZer-2")
error500<-error500[-1,]
#####################

error1000<-matrix(NA,1,3)
for (num in 16:20)
{
  tarea<-param[num,1]
  N<-param[num,2]
  
  fil.res<-paste('Error_TypeI_N_',N,'_tarea_',tarea,'.txt',sep="")
  res<-as.matrix(read.table(fil.res,h=T),M,3)
  error1000<-rbind(error1000,res)
}
colnames(error1000)<-c( "i" ,"SiZer-0","SiZer-2")
error1000<-error1000[-1,]

write.table(error50,'Error_TypeI_N_50.txt')
write.table(error100,'Error_TypeI_N_100.txt')
write.table(error500,'Error_TypeI_N_500.txt')
write.table(error1000,'Error_TypeI_N_1000.txt')