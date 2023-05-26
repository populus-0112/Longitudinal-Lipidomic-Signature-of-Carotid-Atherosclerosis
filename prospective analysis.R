plaque=read.table('plaque_progression.csv',sep=',',header=T)
# change column names
col_name <- scan('plaque_progression.csv', what = "", sep = ",", 
                 nlines = 1, quiet = TRUE, skip = 0, strip.white = TRUE)
colnames(plaque)=col_name

plaque=plaque[!is.na(plaque[,'S4BMI']),]
plaque=plaque[!is.na(plaque[,'S4LDL']),]
plaque=plaque[!is.na(plaque[,'S4GFR_CKD_EPI']),] 

# atherosclerosis progression 
table(plaque[,'EFFSEG4'],useNA='always')
table(plaque[,'AFFSEG5'],useNA='always')

plaque[,'progression']=plaque[,'AFFSEG5']-plaque[,'EFFSEG4']
table(plaque[,'progression'],useNA='always')

plaque[,'progression']=ifelse(plaque[,'progression']>0,1,0)
plaque=plaque[!is.na(plaque[,'progression']),] # delete subjects with missing progression
table(plaque[,'progression'],useNA='always')

plaque$SEX<-relevel(as.factor(plaque$SEX), ref = "M")
plaque$S4SMOKE <- gsub("N","2",plaque$S4SMOKE)
plaque$S4SMOKE <- gsub("E","2",plaque$S4SMOKE)
plaque$S4SMOKE <- gsub("Y","1",plaque$S4SMOKE)
plaque$S4SMOKE<-relevel(as.factor(plaque$S4SMOKE), ref = 2)
plaque$S4ADADM2020<-relevel(as.factor(plaque$S4ADADM2020), ref = "NFG")
table(plaque[,'S4HTN2'],useNA='always')
plaque=plaque[plaque[,'S4HTN2']!='',]
plaque$S4HTN2<-relevel(as.factor(plaque$S4HTN2), ref = "N") 


transform<-function(var){
  plaque[,var]<-(plaque[,var]-mean(plaque[,var],na.rm=T))/sd(plaque[,var],na.rm=T)
}

plaque[,'S4AGE']=transform('S4AGE')
plaque[,'S4BMI']=transform('S4BMI')
plaque[,'S4GFR_CKD_EPI']=transform('S4GFR_CKD_EPI')


plaque=plaque[plaque[,'CENTER']!='OK',] # AZ SD OK
table(plaque[,'progression'],useNA='always')

nm='AC(10:0)'
if (nm%in%colnames(plaque)){
  start_lipid=which(colnames(plaque)=='AC(10:0)')
} else {
  start_lipid=which(colnames(plaque)=='AC.10.0.')
}

selected_lipids_names=colnames(plaque)[start_lipid:(start_lipid+1541)]
col_name=col_name[start_lipid:(start_lipid+1541)]
n_lipids=1542
lipids_select<-plaque[,colnames(plaque)%in%selected_lipids_names]
n_lipids=ncol(lipids_select)
beta<-rep(0,n_lipids)
se<-rep(0,n_lipids)
OR<-rep(0,n_lipids)
upper_CI<-rep(0,n_lipids)
lower_CI<-rep(0,n_lipids)
pvals<-rep(0,n_lipids)
library("geepack")
for(i in 1:n_lipids){
  lipid=scale(lipids_select[,i],scale=TRUE,center=TRUE) 
  fit<-geeglm(progression ~ lipid+S4AGE+SEX+S4SMOKE+S4BMI+S4ADADM2020+S4HTN2+S4GFR_CKD_EPI, id=as.factor(FID), data =plaque, family = binomial(link = "logit"), corstr="independence")
  a=summary(fit)  
  beta[i]=a$coefficients[2,1]               
  OR[i]<-exp(beta[i])
  se[i]=a$coefficients[2,2]
  lower_CI[i]=exp(beta[i] - 1.96 * se[i])
  upper_CI[i]=exp(beta[i] + 1.96 * se[i])
  pvals[i]=a$coefficients[2,4]
  print(i)
}
selected_lipids_names=selected_lipids_names[pvals<0.05]
n_lipids=length(selected_lipids_names)
n_lipids

result=as.data.frame(cbind(colnames(lipids_select),col_name,beta,se,OR,lower_CI,upper_CI,pvals))
write.table(result,'output_gee.csv',sep=',',row.names = F, col.names=T)
