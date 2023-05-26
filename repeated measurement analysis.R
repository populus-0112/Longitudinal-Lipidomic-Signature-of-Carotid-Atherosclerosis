plaque=read.table('plaque_change.csv',sep=',',header=T)
col_name <- scan('plaque_change.csv', what = "", sep = ",", 
                 nlines = 1, quiet = TRUE, skip = 0, strip.white = TRUE)
colnames(plaque)=col_name

# delete prevalent CVD
require(haven)
cvd=read_sas('shs_fcvd2018.sas7bdat')
n_cvd=colnames(cvd)
cvd=cvd[,c("IDNO",'S4CVDFREE','S5CVDFREE','S5CHFFREE','S5MIFREE',"S5DSTKFREE",'ANYCVDS5')]
id_cvd=as.matrix(cvd[,1])
id_plaque=plaque[,1]
ID_set=id_plaque[id_plaque %in% id_cvd]
plaque=plaque[id_plaque %in% ID_set,]
colnames(cvd)%in%colnames(plaque)
plaque=merge.data.frame(plaque,cvd,by.x = "IDNO",by.y="IDNO")
plaque=plaque[plaque[,'S5CVDFREE']==1&plaque[,'S5CHFFREE']==1&plaque[,'S5MIFREE']==1&plaque[,'S5DSTKFREE']==1,]
# delete prevalent CVD


plaque[,'progression']=plaque[,'AFFSEG5']-plaque[,'EFFSEG4']

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

# delete missing values
plaque=plaque[!is.na(plaque[,'IR_c']),]
plaque=plaque[!is.na(plaque[,'BMI_c']),] 
plaque=plaque[!is.na(plaque[,'waist_c']),]
plaque=plaque[!is.na(plaque[,'eGFR_c']),] 
plaque=plaque[!is.na(plaque[,'G0_c']),]
plaque=plaque[!is.na(plaque[,'INSU_c']),]
plaque=plaque[!is.na(plaque[,'SBP_c']),]
plaque=plaque[!is.na(plaque[,'DBP_c']),]

# baseline lipids
c=read.table('plaque_progression.csv',sep=',',header=T)
id_4=as.character(c[,1])
id_o=as.character(plaque[,1])
ID_set=id_4[id_4 %in% id_o] 
length(ID_set)
c=c[id_4 %in% ID_set,]
plaque=plaque[id_o %in% ID_set,]
dim(c)
dim(plaque)
c=c[order(ID_set),]
plaque=plaque[order(ID_set),]
# check
sum(plaque[,1]==c[,1])

table(plaque[,'progression'],useNA='always')
# plaque score: p4: EFFSEG4   p5: AFFSEG5
plaque[,'score_c']=plaque[,'AFFSEG5']-plaque[,'EFFSEG4']

lipid_start0=which(colnames(c)=='AC.10.0.')
lipid_p4=c[,lipid_start0:(lipid_start0+1541)]
# baseline lipids

transform<-function(var){
  plaque[,var]<-(plaque[,var]-mean(plaque[,var],na.rm=T))/sd(plaque[,var],na.rm=T)
}

plaque[,'S4AGE']=transform('S4AGE')
plaque[,'BMI_c']=transform('BMI_c')
plaque[,'waist_c']=transform('waist_c')
plaque[,'G0_c']=transform('G0_c')
plaque[,'SBP_c']=transform('SBP_c')
plaque[,'DBP_c']=transform('DBP_c')
plaque[,'INSU_c']=transform('INSU_c')
plaque[,'IR_c']=transform('IR_c')
plaque[,'eGFR_c']=transform('eGFR_c') 
plaque[,'S4BMI']=transform('S4BMI')
plaque[,'S4INSU']=transform('S4INSU')
plaque[,'S4SBP']=transform('S4SBP')
plaque[,'S4DBP']=transform('S4DBP')
plaque[,'S4IR']=transform('S4IR')

nm='AC(10:0)'
if (nm%in%colnames(plaque)){
  start_lipid=which(colnames(plaque)=='AC(10:0)')
} else {
  start_lipid=which(colnames(plaque)=='AC.10.0.')
}

name=colnames(plaque)[start_lipid:(start_lipid+1541)]
n_lipids=1542

lipids_select<-plaque[,colnames(plaque)%in%name]
####################################
beta<-rep(0,n_lipids)
se<-rep(0,n_lipids)
upper_CI<-rep(0,n_lipids)
lower_CI<-rep(0,n_lipids)
pvals<-rep(0,n_lipids)
library(geepack)
for(i in 1:n_lipids){
  plaque[,'lipid']<-scale(lipids_select[,i],scale=TRUE,center=TRUE)
  plaque[,'lipid_4']<-scale(lipid_p4[,i],scale=TRUE,center=TRUE)
  fit<-geeglm(lipid ~ BMI_c+S4BMI+score_c+lipid_4+S4AGE+SEX+S4SMOKE+S4ADADM2020, id=as.factor(FID), data =plaque[!is.na(plaque[,'lipid']),], family = 'gaussian', corstr="independence")
  a=summary(fit)  
  beta[i]=a$coefficients[2,1]  
  se[i]=a$coefficients[2,2]
  lower_CI[i]=beta[i] - 1.96 * se[i]
  upper_CI[i]=beta[i] + 1.96 * se[i]
  pvals[i]=a$coefficients[2,4]
  print(i)
}
sum(pvals<0.05)
result=as.data.frame(cbind(colnames(lipids_select),beta,se,lower_CI,upper_CI,pvals))
write.table(result,'output_BMI.csv',sep=',',row.names = F, col.names=T)

