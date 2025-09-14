
# nlasso with 500 reps------------
# scaling before running the model
X<-data.matrix(scale_lassodata)
Y<-mydata_dummy_compl3$treatment
# create two columns data matrix for non-binary outcome if using binomial lasso model 
Y<-cbind(1-mydata$outcome,mydata$outcome)
data<-data.frame(matrix(0,ncol=500,nrow =82))
for (n in 1: 500){
  cv_output<-cv.glmnet(X,as.matrix(Y),family= "binomial",alpha=1,nfolds =4,
                       type.measure = "deviance")
  # set unpenalised factors
  #penalty.factor=c(0,1,0,1,0,rep(1,46),0,0,rep(1,13),0,1,rep(0,5),rep(1,8),0))
  beta_full<-t(as.matrix(coef(cv_output,s="lambda.min")))
  beta_full<-beta_full[, colnames(X)]
  selected <- ifelse(beta_full != 0, yes = 1, no = 0)
  print(selected)
  print(beta_full)
  data[,n]<-ifelse(beta_full != 0,beta_full,0)
  n=n+1
}

#write.csv(data,"500_lasso_CV.csv")
# compute selection proportions for each varb
data<-round(data,3)
data$names<-colnames(scale_lassodata)
data1<-melt(data,id.vars="names")
selected_per<-round(1-rowSums(data==0)/500,2)
selection<-as.data.frame(cbind(selected_per,data$names))
selection<-selection[order(selection$selected_per),]
selection<-selection%>%rename(varb=V2)
selection$oder<-c(1:ncol(X))

# CTMLE-----------------
# partial correlation 
#initial estimate
unadjusted_results<-mydata%>%group_by(treatment_id)%>%
  summarise_at(vars(Y11:Y31),mean,na.rm = TRUE)
baseline_results<-mydata%>%group_by(treatment_id)%>%
  summarise_at(vars(Y10:Y30),mean,na.rm = TRUE)
Ydata1<-mydata$Y11
Q<-cbind(rep(unadjusted_results[1,4],length(Ydata)),rep(unadjusted_results[2,4],length(Ydata)))
R.data<-data.frame(mydata$Y11,mydata$treatment_id)
R.data$Q0<-ifelse(mydata$treatment_id==1,unadjusted_results[2,4],unadjusted_results[1,4])
R.data$Q0<-as.numeric(R.data$Q0)
R.data$R<-R.data$mydata.Y11-R.data$Q0
# cbind with covariates set W
R.data2<-cbind(R.data,W)
with(mydata,by(Y11,treatment_id,summary))
W<-as.matrix(W)
A<-R.data$mydata.treatment_id
R<-R.data$R
partial_corr<-function(column){
  pcor.test(R,column,A,method = "pearson")$estimate
}

partial_cor_results<-apply(W,2,partial_corr)  
print(partial_cor_results)
partial_cor_results<-as.data.frame(partial_cor_results)
partial_cor_results$abs<-abs(partial_cor_results$partial_cor_results)
partial_cor_results$order<-rank(-partial_cor_results$abs)
order_cov<-partial_cor_results$order
ordered_partial<-partial_cor_results[order(-partial_cor_results[,2]),]
ordered_partial$order<-c(1:ncol(W))
# run CTMLE
Ydata<-mydata$Y11
time_scable<-system.time({
  ctmle_discrete_fit2 <- ctmleDiscrete(Y = Ydata, A = A, W = W,
                                       preOrder = TRUE, Qbounds=c(0,1),
                                       order=order_cov,detailed = TRUE)
})
summary(ctmle_discrete_fit)


# LOGISTIC strategy
logistic_ctmle_data<-W
logistic_ctmle_data$treatment_id<-mydata$treatment_id
logistic_ctmle_data$Y<-mydata$Y11
logistic_varb_list<-colnames(W)
mse_all<-numeric(length(logistic_varb_list))
for (i in 1: length( logistic_varb_list)){
  varb<- logistic_varb_list[i]
  formula_string<-paste("treatment_id","~",varb)
  formula<-as.formula(formula_string)
  model<-glm(formula,data=logistic_ctmle_data,family = "binomial")
  #print(summary(model)$deviance)
  predict.pscore<-predict(model,type = "response")
  logistic_ctmle_data$H.treated<-logistic_ctmle_data$treatment_id/predict.pscore
  logistic_ctmle_data$H.untreated<-(1-logistic_ctmle_data$treatment_id)/(1-predict.pscore)
  logistic_ctmle_data$clever<-logistic_ctmle_data$H.treated-logistic_ctmle_data$H.untreated
  formula_string2<-paste("Y","~treatment_id+",varb)
  formula2<-as.formula(formula_string2)
  Q.model<-glm(formula2,data=logistic_ctmle_data)
  init.pred<-predict(Q.model,type = "response")
  eps.mod<-lm(Y~-1+H.treated+H.untreated+offset(init.pred),data=logistic_ctmle_data)
  predict.Q<-predict(eps.mod)
  residuals<-logistic_ctmle_data$Y-predict.Q
  mse<-mean(residuals^2)
  mse_all[i]<-mse
  print(mse_all)
  
}

logistic_result<-data.frame(mse_all,logistic_varb_list)
logistic_result$order<-rank(logistic_result$mse_all)
order_cov_logistic<-logistic_result$order

time_scable<-system.time(
  ctmle_discrete_fit3 <- ctmleDiscrete(Y = Y, A = A, W = W, 
                                       preOrder = TRUE, Qbounds=c(0,1),
                                       order=order_cov_logistic,detailed = TRUE)
)
summary(ctmle_discrete_fit2)

# user-defined function for sharp package------------
BinomialLasso <- function(xdata,ydata,penalty.factor,Lambda,family=NULL,group_x=NULL){
  #print(ydata)
  #print(mean(ydata))
  ydata<-cbind(1-ydata, ydata)
  mymodel <- glmnet(xdata,as.matrix(ydata),family= "binomial",lambda = Lambda,penalty.factor=penalty.factor)
  #,penalty.factor=c(0,1,0,1,0,rep(1,46),0,0,rep(1,13),0,1,rep(0,5),rep(1,8),0)
  beta_full <- t(as.matrix(mymodel$beta))
  beta_full <- beta_full[, colnames(xdata)]
  #beta_full<-beta_full[which(names(beta_full)!="(Intercept)")]
  selected <- ifelse(beta_full != 0, yes = 1, no = 0)
  
  return(list(selected = selected, beta_full = beta_full))

}
# checking
X<-data.matrix(scale_lassodata)
Y<-mydata$Y11
mygrid <- LambdaSequence(lmax = 10, lmin = 1.0e-05, cardinal = 200)
myBinomialLasso<-BinomialLasso(xdata=X,ydata=Y,Lambda=mygrid)
# run stability-enhanced lasso model
MYstab<-VariableSelection(
  xdata=X,
  ydata =Y,
  Lambda = mygrid,
  #penalty.factor=c(0,1,0,1,0,rep(1,45),0,0,rep(1,8),0,rep(1,4),0,rep(1,18)),
  K = 500,
  tau = 0.5,
  seed = 1001,
  n_cat = 3,
  resampling = "subsampling",
  pi_list = seq(0.5, 0.99, by = 0.01),
  family = "gaussian",
  group_x = NULL,
  implementation = BinomialLasso,
  verbose = TRUE
)
pdf(file = "cab_BinomialLasso.sharp.pdf",width = 16, height = 20)
par(mar=c(8,8,6,12))
CalibrationPlot(MYstab)
dev.off()

pdf(file = "stab_BinomialLasso.sharp.pdf",width = 16, height = 20)
par(mar=c(8,8,6,12))
plot(MYstab)
dev.off()
