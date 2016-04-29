library(dataedat)
library(lme4)
library(synthpop)
library(R2jags)
library(SweaveLst)


predictor.matrix1=matrix(c(1,0,1,1,0,1,1,1,0,0,0,1,0,0,0,0),4,4)
predictor.matrix2=matrix(c(1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0),4,4)
predictor.matrix3=matrix(c(1,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0),4,4)

tabs<-lapply(list(predictor.matrix1,predictor.matrix2,predictor.matrix3),function(pred){
  syn(tab,predictor.matrix=pred,method=c("","","parametric","parametric"),drop.pred.only=FALSE,
          maxfaclevels=nlevels(tab$SCH_ID))$syn})


model=as.formula("BYTXMIRR~BYSEX+BYG10EP+(BYG10EP|SCH_ID)")



U0=lmer(model,data=tab)
save(U0,file="U0.rda")
load("U0.rda")

syn1=lapply(tabs,function(tab1){summary(lmer(model,data=tab1))})
save(syn1,file="syn1.rda")
load("syn1.rda")



nrep=1000
synrep<-lapply(list(predictor.matrix1,predictor.matrix2,predictor.matrix3),function(pred){
  syn(tab,predictor.matrix=pred,method=c("","","parametric","parametric"),drop.pred.only=FALSE,
          maxfaclevels=751,m=nrep)$syn})


tab2<-tab

U0s<-summary(U0)
Sigma<-U0s$varcor$SCH_ID[,]
#a=eigen(Sigma)
toto<-plyr::rlply(nrep,(function(){
  randome<-mvrnorm(n=nlevels(tab$SCH_ID),mu=rep(0,2),Sigma=Sigma)
#  b=eigen(var(randome))
#  randome<-a$vectors%*%diag(sqrt(a$values))%*%t(a$vectors)%*%b$vectors%*%diag(1/sqrt(b$values))%*%t(b$vectors)%*%randome
  randome<-randome-rep(apply(randome,2,mean),each=nrow(randome))
  
    tab2$BYTXMIRR<-model.matrix(~BYSEX+BYG10EP,data=tab)%*%U0s$coefficients[,1,drop=FALSE]+rnorm(nrow(tab),0,U0s$sigma)+
    model.matrix(~0+SCH_ID,data=tab)%*%randome[,1]+(model.matrix(~0+SCH_ID,data=tab)*tab$BYG10EP)%*%randome[,2]
  tab2})())
synrep<-c(synrep,list(toto))

save(synrep,file="synrep.rda")


estimates<-plyr::llply(synrep,function(x){plyr::llply(x,function(y){summary(lmer(model,data=y))})})
names(estimates)<-1:4
varestimates<-plyr::ldply(estimates,function(x){plyr::laply(x,function(y){unlist(y$varcor)})})
save(estimates,varestimates,file="estimates.rda")
load("estimates.rda")
names(varestimates)[1]<-c("Method")

library(ggplot2)
library(plyr)
plot1<-ggplot(data = varestimates,aes(x=SCH_ID1,color=Method))+geom_density()+geom_vline(xintercept = unlist(summary(U0)$varcor[1])[1])+xlab("$\\hat{\\sigma}^2_{\\gamma}$")+
  geom_vline(data=ddply(varestimates,~Method,function(x){mean(x$SCH_ID1)}), aes(xintercept=V1,  colour=Method))
x=graphtikzcode("print(plot1)",scale=.8,yxratio=.7,caption="Density plot of REML estimates of the variance of the intercept coefficients between schools")

y=stargazer2(U0,style="ajps",title="REML estimates")
y2=stargazer2(summary(U0)$varcor,style="ajps",title="REML estimates")

#save(U0,x,y,file="/home/daniel/Dropbox/Travail/Recherche/Projets/Synthetic_data/Cours_JPSMSurv662/x.rda")





