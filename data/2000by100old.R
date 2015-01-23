# Generate simulated data set 2000 patients x 100 variables (additional IV, risk factors, colliders, interrcorrel noise)
# and 2000 patients x 10 variables (Setoguchi 2008, Scenario E)
# 
# 01-Apr-14 Yen Low
###############################################################################

#############SET PARAMETERS######################
rm(list=ls(all=T))
setwd("/mnt/hgfs/projects/psm/data")
getwd()

source("../hdpsFunctions.R")


n=2000 #number of patients
p=100 #number of variables
#set logit(exposure) coefficients; includes b0
b=c(0,0.8,-0.25,0.6,-0.4,-0.8,-0.5,0.7,0,0,0,0.1,0.2,0.8,0.6,0,0,0,0,0,0,0,0,rep(0,100-22))
#set logit(outcome) coefficients; includes a0
a=c(-3.85,0.3,-0.36,-0.73,-0.2,0,0,0,0.71,-0.19,0.26,0,0,0,0,0.1,0.2,0.6,0.8,0,0,0,0,rep(0,100-22))
gamma=-0.4 #coeff for ps covariate adjustment

set.seed(123)
#generate 10 variables based on standard normal distribution
#introduce intercorrelation between 0.2 and 0.9
v1=rnorm(n,0,1)
v3=rnorm(n,0,1)

x1=as.integer(v1>=0)
x2=rnorm(n,0,1)
x3=as.integer(v3>=0)
x4=rnorm(n,0,1)

v5=sample(c(1,-1),1)*((1-sqrt(0.7))*v1+sqrt(0.7)*rnorm(n,0,1))
v6=sample(c(1,-1),1)*((1-sqrt(0.1))*x2+sqrt(0.1)*rnorm(n,0,1))
v8=sample(c(1,-1),1)*((1-sqrt(0.7))*v3+sqrt(0.7)*rnorm(n,0,1))
v9=sample(c(1,-1),1)*((1-sqrt(0.1))*x4+sqrt(0.1)*rnorm(n,0,1))

x5=as.integer(v5>=0)
x6=as.integer(v6>=0)
x7=rnorm(n,0,1) #strong IV
x8=as.integer(v8>=0)
x9=as.integer(v9>=0)
x10=rnorm(n,0,1) #weak risk factor

#additional variables
x11=rnorm(n,0,1) #weak IV
x12=rnorm(n,0,1) #weak IV
x13=rnorm(n,0,1) #strong IV
x14=rnorm(n,0,1) #strong IV
x15=rnorm(n,0,1) #weak risk factor
x16=rnorm(n,0,1) #weak risk factor
x17=rnorm(n,0,1) #strong risk factor
x18=rnorm(n,0,1) #strong risk factor

#create colliders with various correl to IV or risk factor
x19=sample(c(1,-1),1)*((1-sqrt(0.7))*x11+(1-sqrt(0.7))*x15+(1-2*(1-sqrt(0.7)))*rnorm(n,0,1))
x20=sample(c(1,-1),1)*((1-sqrt(0.7))*x12+(1-sqrt(0.1))*x17+(sqrt(0.7)+sqrt(0.1)-1)*rnorm(n,0,1))
x21=sample(c(1,-1),1)*((1-sqrt(0.1))*x13+(1-sqrt(0.7))*x18+(sqrt(0.7)+sqrt(0.1)-1)*rnorm(n,0,1))
x22=sample(c(1,-1),1)*((1-sqrt(0.1))*x14+(1-sqrt(0.1))*x16+(1-2*(1-sqrt(0.1)))*rnorm(n,0,1))

#binarize some variables
x11=as.integer(x11>=0.5) #weak IV
x13=as.integer(x13>=0.5)#strong IV
x16=as.integer(x16>=0.5) #weak risk factor
x17=as.integer(x17>=0.5) #strong risk factor


#create noise variables
xmat=cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22)
dim(xmat)
for(i in 23:50) xmat=cbind(xmat,rnorm(n,0,1)) #std normal noise
for(i in 51:p) xmat=cbind(xmat,round(runif(n,0,1))) #binary noise
colnames(xmat)=paste("x",1:p,sep="")

#set intercorrelations among noise
xmat[,"x23"]=(1-sqrt(0.7))*xmat[,"x50"]+sqrt(0.7)*rnorm(n,0,1)
xmat[,"x24"]=(1-sqrt(0.1))*xmat[,"x49"]+sqrt(0.1)*rnorm(n,0,1)
xmat[,"x51"]=as.integer(((1-sqrt(0.2))*xmat[,"x100"]+sqrt(0.2)*rnorm(n,0,1))>=0)
xmat[,"x52"]=as.integer(((1-sqrt(0.05))*xmat[,"x99"]+sqrt(0.05)*rnorm(n,0,1))>=0)

#view intercorrelations
#dim(xmat)
#cormat=cor(xmat)
#write.table(cormat,file="100x100cormat.txt",sep="\t",col.names=F, row.names=F)

#png("100x100cormat.png",width=7,height=7,res=100,units="in", bg="transparent")
#heatmap(cormat,Rowv=NA,symm=T,col=cm.colors(10),scale="none")
#dev.off()

#cormat[abs(cormat)<0.1]=NA #suppress for easy viewing

#save(xmat,file="xmat2000by100.RData")

#Simulate propensity scores based on scenario E based on eq in Appendices 1-2
#handle b0 (i.e. b[1]) differently
Pexp=c(1/(1+exp(-(	b[1] + b[-1] %*% t(xmat)
						+b[3]*x2*x2
						+0.5*b[2]*x1*x3
						+0.7*b[3]*x2*x4
						+0.5*b[5]*x4*x5
						+0.5*b[6]*x5*x6))))
#sum(mapply(function(x,y) x*y, b[-1], xmat[2,]))

#simulate true exposure
rdmNum=runif(n)
temp=as.data.frame(cbind(rdmNum,Pexp))
temp$exposure=0
temp$exposure[temp$rdmNum<temp$Pexp]=1

#Simulate estimated prob of outcome based on eq in Appendices 1-2
Poutcome=c(1/(1+exp(-(a[1]+a[-1] %*% t(xmat) + gamma*temp$exposure))))

#simulate true outcome
temp$rdmNum2=runif(n)
temp$Poutcome=Poutcome
temp$outcome=0
temp$outcome[temp$rdmNum2<temp$Poutcome]=1

########## put matrices together
ds=cbind(outcome=temp$outcome,exposure=temp$exposure,ps=temp$Pexp,Poutcome=temp$Poutcome,
				as.data.frame(xmat))
ds10=ds[,1:14]
		
###### save
#save(ds,temp,file="2000x100.RData")
#load(file="../../data/2000x100.RData")

#get true OR (0.67)
exp(gamma) #pre-defined OR (0.67)
mod=glm(paste("outcome ~ exposure +", paste("x",c(1:4,8:10,15:18),collapse=" + ",sep="")),
				data=ds,family="binomial")
coefficients(mod)["exposure"]
exp(coefficients(mod)["exposure"])   #actual re-derived OR from data = 0.45

#get crude OR (0.58)
tab=table(temp$outcome,temp$exposure)
#chisq.test(tab)
OR=tab[1,1]*tab[2,2]/tab[1,2]/tab[2,1]

#get true standardized mean difference (-0.09)
trueSmd=smd(ds10,exposed=exposed,outcome=outcome)

#check underlying distributions
#png("2000x100histograms.png",width=7,height=11,res=100,units="in", bg="transparent")
#par(mfrow=c(3,2))
#hist(Pexp)
#h=hist(logit(Pexp),prob=T)
#lines(density(logit(Pexp)))

#hist(Poutcome)
#hist(logit(Poutcome),prob=T)
#lines(density(logit(Poutcome)))

#hist(temp$exposure)
#hist(temp$outcome)
#dev.off()
