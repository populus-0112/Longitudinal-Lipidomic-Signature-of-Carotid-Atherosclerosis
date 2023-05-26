plaque=read.table('cs_4.csv',sep=',',header=T)
dim(plaque)
# change column names
col_name <- scan('cs_4.csv', what = "", sep = ",", 
                 nlines = 1, quiet = TRUE, skip = 0, strip.white = TRUE)
colnames(plaque)=col_name

# atherosclerosis incidence 
table(plaque[,'ATHER4'],useNA='always')

plaque=plaque[!is.na(plaque[,'ATHER4']),] 

plaque$SEX<-relevel(as.factor(plaque$SEX), ref = "M")
plaque$S4SMOKE <- gsub("N","2",plaque$S4SMOKE)
plaque$S4SMOKE <- gsub("E","2",plaque$S4SMOKE)
plaque$S4SMOKE <- gsub("Y","1",plaque$S4SMOKE)
plaque$S4SMOKE<-relevel(as.factor(plaque$S4SMOKE), ref = 2)
plaque$S4ADADM2020<-relevel(as.factor(plaque$S4ADADM2020), ref = "NFG")
table(plaque[,'S4HTN2'],useNA='always')
plaque=plaque[plaque[,'S4HTN2']!='',]
plaque$S4HTN2<-relevel(as.factor(plaque$S4HTN2), ref = "N")

plaque=plaque[!is.na(plaque[,'S4BMI']),]
plaque=plaque[!is.na(plaque[,'S4AGE']),]
plaque=plaque[!is.na(plaque[,'S4GFR_CKD_EPI']),] 

transform<-function(var){
  plaque[,var]<-(plaque[,var]-mean(plaque[,var],na.rm=T))/sd(plaque[,var],na.rm=T)
}

plaque[,'S4AGE']=transform('S4AGE')
plaque[,'S4BMI']=transform('S4BMI')
plaque[,'S4GFR_CKD_EPI']=transform('S4GFR_CKD_EPI') 

nm='AC(10:0)'
if (nm%in%colnames(plaque)){
  start_lipid=which(colnames(plaque)=='AC(10:0)')
} else {
  start_lipid=which(colnames(plaque)=='AC.10.0.')
}

##standardize the lipids
for(i in start_lipid:(start_lipid+1541)){
  plaque[,i]<-scale(plaque[,i],scale=TRUE,center=TRUE)
}

name=colnames(plaque)[start_lipid:(start_lipid+1541)]
n_lipids=1542
lipids_select<-plaque[,colnames(plaque)%in%name] 

####################################
beta<-rep(0,n_lipids)
se<-rep(0,n_lipids)
OR<-rep(0,n_lipids)
upper_CI<-rep(0,n_lipids)
lower_CI<-rep(0,n_lipids)
pvals<-rep(0,n_lipids)
library("geepack")
for(i in 1:n_lipids){
  lipid<-lipids_select[,i]   
  fit<-geeglm(ATHER4~lipid+S4AGE+SEX+S4SMOKE+S4BMI+S4ADADM2020+S4HTN2+S4GFR_CKD_EPI, id=as.factor(FID), data =plaque, family = binomial(link = "logit"), corstr="independence")
  a=summary(fit)        
  beta[i]=a$coefficients[2,1]
  OR[i]<-exp(beta[i])
  se[i]=a$coefficients[2,2]]
  lower_CI[i]=exp(beta[i] - 1.96 * se[i])
  upper_CI[i]=exp(beta[i] + 1.96 * se[i])
  pvals[i]=a$coefficients[2,4]
  print(i)
}
result=as.data.frame(cbind(colnames(lipids_select),beta,se,OR,lower_CI,upper_CI,pvals)) 
write.table(result,'output_p4.csv',sep=',',row.names = F, col.names=T)
