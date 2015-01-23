# To simulated data 2000x10, do expertPS, hdPS, lasso logit, RF,
# Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman correl
# Nsim=1000
# Collect pval, smd of baseline variables
#
# 15-Oct-14 use rfPS in reg mode (set min node size to 5% of sample size, i.e. 100, Malley 2012) 
# 11-Sep-14 adjust rfPS by OOB error by sourcing rfFunctions.R
# 26-Aug-14 Yen Low
##############################################################################
require(rJava)
require(Matching)
require(Epi) #clogistic
require(glmnet)
require(randomForest)
require(vegan)
require(FNN)
require(Matrix)

#parallel core computing for lassoMV bootstrapping
require(doParallel)
require(foreach)
#specify number of cores req
cl=makeCluster(10) 
registerDoParallel(cl)

#############SET PARAMETERS######################
#RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
#ProjectPath=Sys.getenv("ProjectPath") #Project is an environment variable, path="/mnt/hgfs/projects"
RScriptPath="~/scripts"
ProjectPath="~/projects"
setwd(paste(ProjectPath,"/psm/models/sim2000x10",sep=""))
getwd()

baseDir = getwd();
tempDir = paste(baseDir,"tmp",sep="/") #output directory
patientfile=paste(baseDir,"patients.txt",sep="")
dir.create(tempDir)

# initialize the Java subsystem.  include a jdbc object
.jinit(classpath="/home/yenlow/software/pharmacoepi.jar", force.init = TRUE,parameters="-Xmx2g");
.jclassPath()

#required functions
source(paste(RScriptPath,"/R/logit.R",sep=""))
source(paste(RScriptPath,"/R/patientChar.R",sep=""))
source("../../hdpsFunctions.R")
#source("../../rfFunctions.R")   #only req if using padjfromRF

id="id"
exposed="exposure" #name of exposure variable (must be string)
outcome="outcome"
expertvariables=c("x1","x2","x3","x4","x5","x6")  #exclude x7, IV
outcomemodelvariables_exclmatchingvar=c("x8","x9","x10") 
#empirical variables may be matched depending on HD-PS
empvariables=paste("x",1:10,sep="")
empvariables_cat=c("x1","x3","x5","x6","x8","x9")
empvariables_num=setdiff(empvariables,empvariables_cat)

#formula for expertPS model: exposure ~ var1 + var 2 + ...
fPS = as.formula(paste(exposed," ~ ",paste(expertvariables,collapse=" + "),sep=""))

#formula for outcome after expertPS
fmod=paste("outcome ~ exposure +", paste(outcomemodelvariables_exclmatchingvar,collapse=" + ",sep=""))

fPS
fmod

desiredOrder=c("n1","ORmatched","ORlow_matched","ORupp_matched","coeff_matched","se_matched","bias_matched","inCI_matched",
               "n0","ORadj","ORlow_adj","ORupp_adj",            "coeff_adj","se_adj",       "bias_adj","inCI_adj",
               "SMD","KS","KL","Cstat")

#to store consolidated results
startAll=proc.time()[1:3]

Nsim=1000
Ndraw=100
Nmethods=14 #14 methods (incl. 3 rdm sampling),
resultsArray=array(dim=c(Nmethods,20,Nsim))  #12methods (incl rdm downsampling), 20 metrics
#drawArray=array(dim=c(Nmethods,16,Ndraw))  #3methods (incl rdm downsampling), 20-4 metrics
pmatArray=array(dim=c(Nmethods+1,10,Nsim))  #12methods (incl rdm downsampling)+before, 10 variables
smdmatArray=pmatArray  #12methods (incl rdm downsampling), 10 variables
Noutcomesmat=matrix(nrow=Nsim,ncol=2)
time=array(dim=c(Nmethods,3,Nsim))
dimnames(time)[[1]]=c("expertPS","hdPS","lassoPS","rfPS","euclidean",
                      "jaccard","dice","cosine","pearson","spearman","lassoMV","rdmSamp1x","rdmBag","rdmJK")
dimnames(time)[[2]]=c("user","system","elapsed")
dimnames(resultsArray)[[1]]=dimnames(time)[[1]]
#dimnames(drawArray)[[1]]=dimnames(time)[[1]]
#dimnames(drawArray)[[2]]=c("n0","n1","KS","KL","Cstat","SMD",
#                           "ORadj","ORlow_adj","ORupp_adj","coeff_adj","se_adj",
#                           "ORmatched","ORlow_matched","ORupp_matched","coeff_matched","se_matched")

rfImpt=matrix(nrow=10,ncol=Nsim)  #10 rows for x1-x10
rownames(rfImpt)=paste("x",1:10,sep="")
lassoPSBeta=rfImpt
lassoMVvar=wantedVarList=matchedID=var_corOutcomeList=rdmBag=rdmJK=list()
Nexposed=trueORvec=trueSMDvec=trueCoeffvec=bestlambda=c()


#foreach(i=1:Nsim) %doPar%{
for(i in 1:Nsim){
	
cat("\n\n\n------------------------",i,"---------------------\n\n")

#generate simulated data
source("../../data/2000by10.R")
#load("../../data/2000x10.RData",verbose=T)
#rm(list=ls()[!(ls() %in% "ds10")]) #remove all objects except ds10
ds10$id=1:nrow(ds10)
nminor=min(table(ds10[,exposed]))

trueORvec[i]=trueOR
trueSMDvec[i]=trueSmd[2]
trueCoeffvec[i]=trueCoeff
Nexposed[i]=sum(ds10[,exposed])
Noutcomesmat[i,]=Noutcomes

###### expertPSM ######
#scale continuous variables (already scaled to standard normal)

#calculate estimated PS
print("############### expertPS ############")

start=proc.time()[1:3]

expertPsMod=glm(fPS, data=ds10[,c(id,exposed,expertvariables)], family="binomial")
expertResults=extractResults(ps=expertPsMod$fitted.values,exposurevec=expertPsMod$y,fmod=fmod,
							data=ds10,id=id,exposed=exposed, outcome=outcome,logitFlag=TRUE,
							outfile=NULL, verbose=FALSE)

end=proc.time()[1:3]
time["expertPS",,i]=end-start

####### HDPS ######################
print("############### hdPS ############")

start=proc.time()[1:3]

#prepare data into format required by hdps
#read data file as single string (to use addPatientsFromBuffer)
patientheader=paste(id,exposed,outcome,sep="\t")
patientstring=paste(paste(ds10[,id],ds10[,exposed],ds10[,outcome],sep="\t"),collapse="\n")
#for reading in whole file as single string
#patientstring=readChar(patientfile, file.info(patientfile)$size)
#patientstring=gsub("\"","",patientstring)
#patients=dataset[,c("pid_org",exposed,outcome[1])] #get patient id, exposure, outcome
#write.table(patients,file=patientfile,sep="\t",col.names=T,row.names=F,na="") #output data to patient.txt;
datainstring=paste(patientheader, patientstring, sep="\n")
variables=empvariables_cat
dimdata=ds10[,c(id,variables)]

### INVOKE pharmacoepi.jar ### (handles cat variables only)
hdpsobj=hdps(datainstring,dimdata,outDir=tempDir,Nmostfreq=10,k=10,stratifyDim=FALSE,
		outfile="output_cohort.txt",FullOutput=TRUE,verbose=T,ZeroCellCorrection=F)
hdpsobj$selectedvariables #(see Fig2, Schneeweiss 2008 for recoded definitions)
wantedvar=hdpsobj$selectedvariables[grep("1Once$",hdpsobj$selectedvariables)]

var_corExposed=empvariables_num[abs(cor(ds10[,exposed],ds10[,empvariables_num]))>0.05]
var_corOutcome=empvariables_num[abs(cor(ds10[,outcome],ds10[,empvariables_num]))>0.05] #include in PS
IV=setdiff(var_corExposed,var_corOutcome) #exclude from PS

# Estimate the PS (force numerical variables into PS model)
dataPS=cbind(ds10[,c(id,exposed,var_corOutcome)],hdpsobj$hdpsdata[,wantedvar])
hdPsMod=glm(paste(exposed,"~ . -",id), data=dataPS, family="binomial")
names(hdPsMod$fitted.values)=as.character(dataPS[,id])
summary(hdPsMod$fitted.values)

hdResults=extractResults(ps=hdPsMod$fitted,exposurevec=hdPsMod$y,fmod=NULL,
                      		data=ds10,id=id,exposed=exposed, outcome=outcome,logitFlag=TRUE,
                      		outfile=NULL,verbose=FALSE)

end=proc.time()[1:3]
time["hdPS",,i]=end-start


####### lasso ##################
#"exposure ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10"
#lassoPs=glmnet(paste(exposed, "~", paste("x",1:10,collapse=" + ",sep="")), data=ds10, family="binomial")

print("############### lassoPS ############")

start=proc.time()[1:3]

xmat=as.matrix(ds10[,paste("x",1:10,sep="")])
#penalty.factor=c(0,rep(1,ncol(xmat)-1)) #set 0 to force variable into model
#tune lambda for glmnet using 5-fold CV
lassoPsmod=cv.glmnet(xmat,ds10[,exposed],alpha=1,family="binomial",standardize=F,nfold=5)
bestlambda[i]=lassoPsmod$lambda.1se
lassoPSBeta[,i]=coeffAtlambda(lassoPsmod)[-1]  #exclude intercept

#get estimated ps
psl=unlogit(as.numeric(predict(lassoPsmod,xmat,s=lassoPsmod$lambda.1se))) #in logit form
names(psl)=ds10[,id]

lassoResults=extractResults(ps=psl,exposurevec=ds10[,exposed],fmod=NULL,
                        		data=ds10,id=id,exposed=exposed, outcome=outcome,logitFlag=TRUE,
                        		outfile=NULL,verbose=FALSE)

end=proc.time()[1:3]
time["lassoPS",,i]=end-start

##if verbose
#plot lassmod.cv
#plot(lassoPsmod)
#plot(lassoPsmod$lambda,lassoPsmod$glmnet.fit$dev.ratio)
#abline(v=lassoPsmod$lambda.1se)

#plot beta coeff
#coeff=na.omit(coeffAtlambda(lassoPsmod))
#coeff
#par(mar=c(1,15,1,1))
#boxplot(t(coeff),horizontal=T,las=2,cex.axis=0.7,pch=20)
#abline(v=0,lty="dotted")
#summary(ps)

####### RF ##################
print("############### rfPS ############")

start=proc.time()[1:3]

rfPsMod=randomForest(xmat,ds10[,exposed],ntree=100,importance=T,nodesize=100)  #reg mode for rel class freq
#ps=predict(rfPsMod, type="prob")[,2]
#ps=scales::rescale(ps,to=c(0.001,0.999))   #Method 1: rescale to avoid zeros and 1 which become Inf after logit
#ps=padjfromRF(rfPsMod)   #Method 2: oob-adjust prob from RF to avoid 100% prob
ps=rfPsMod$predicted
names(ps)=ds10[,id]

rfImpt[,i]=rfPsMod$importance[,1]

tryobj=try(extractResults(ps=ps,exposurevec=ds10[,exposed],
                         data=ds10,fmod=NULL,id=id,exposed=exposed, outcome=outcome,logitFlag=TRUE,
                         outfile=NULL,verbose=FALSE))
if(class(tryobj)!="try-error") rfResults=tryobj else rfResults=NA
end=proc.time()[1:3]
time["rfPS",,i]=end-start

##if verbose
#par(mfrow=c(2,2))
#hist(ps[ds10[,exposed]==0],xlim=0:1)
#hist(ps[ds10[,exposed]==1],xlim=0:1)
#summary(ps[ds10[,exposed]==0])
#summary(ps[ds10[,exposed]==1])

####### Similarity ##################
print("############### similarity ############")

xmat_ctrl=xmat[ds10[,exposed]==0,]
xmat_trted=xmat[ds10[,exposed]==1,]
rownames(xmat_ctrl)=ds10[ds10[,exposed]==0,id]
rownames(xmat_trted)=ds10[ds10[,exposed]==1,id]

runSimPS<-function(method="jaccard",caliper=0.7,nsd=3,algorithm="kd_tree"){
  matchedSim=matchByDist(xmat_ctrl,xmat_trted,method=method,k_neighbors=5,caliper=caliper,nsd=nsd,algorithm=algorithm)
  simResults=extractResults(ps=matchedSim,exposurevec=NULL,
                                data=ds10,fmod=NULL,id=id,exposed=exposed, outcome=outcome,logitFlag=F,
                                outfile=NULL,verbose=FALSE)
  return(simResults)
}

start=proc.time()[1:3]
euclideanResults=runSimPS(method="euclidean",nsd=3,algorithm="brute")
end=proc.time()[1:3]
time["euclidean",,i]=end-start

start=proc.time()[1:3]
jaccardResults=runSimPS(method="jaccard",caliper=0.1)
end=proc.time()[1:3]
time["jaccard",,i]=end-start

start=proc.time()[1:3]
diceResults=runSimPS(method="dice",caliper=0.1)
end=proc.time()[1:3]
time["dice",,i]=end-start

start=proc.time()[1:3]
cosineResults=runSimPS(method="cosine",caliper=0.1)
end=proc.time()[1:3]
time["cosine",,i]=end-start

start=proc.time()[1:3]
pearsonResults=runSimPS(method="pearson",caliper=0.1)
end=proc.time()[1:3]
time["pearson",,i]=end-start

start=proc.time()[1:3]
spearmanResults=runSimPS(method="spearman",caliper=0.1)
end=proc.time()[1:3]
time["spearman",,i]=end-start


####### lasso multivariate model (NOT PS) ##################
#"exposure ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10"
#lassoPs=glmnet(paste(exposed, "~", paste("x",1:10,collapse=" + ",sep="")), data=ds10, family="binomial")

print("############### lasso MV model ############")

start=proc.time()[1:3]

xmat=Matrix(as.matrix(ds10[,c(exposed,paste("x",1:10,sep=""))]),sparse=T)
penalty.factor=c(0,rep(1,ncol(xmat)-1)) #set 0 to force variable into model
#tune lambda for glmnet using 5-fold CV
lassoMvmod=cv.glmnet(xmat,ds10[,outcome],alpha=1,family="binomial",standardize=F,nfold=5,penalty.factor=penalty.factor)

coeff=na.omit(sort(coeffAtlambda(lassoMvmod)))

#bootstrap, repeat 100 times for 95%CI
#use global lambda (lambda.1se from lassomod.cv)
lassomodbootstrap=genCI(xmat,ds10[,outcome],ntrials=100,
                        lambda=lassoMvmod$lambda.1se,alpha=1,penalty.factor=penalty.factor)
lassomodCI=setCL(lassomodbootstrap)
lassoMVvar[[i]]=lassomodCI[[1]]  #get coeff

#jack-knife 80% sample, repeat 5*200 times
#use global lambda (lambda.1se from lassomod.cv)
#lassomodjacknife=genCI(xmat,ds10[,outcome],nfold=5,ntrials=200,replacement=FALSE,
#                       type="glmnet",lambda=lassoMvmod$lambda.1se,alpha=1,penalty.factor=penalty.factor)
#lassomodCI=setCL(lassomodjacknife)

#convert beta coeff to OR
#round(lassomodCI[[1]],3)
ORCInonZero=as.data.frame(exp(lassomodCI$betaCIlim[lassomodCI$beta_nonZero,]))

#png("lassoOR.png",width=7,height=11,units="in",res=300,bg="transparent")
#forestplot(ORCInonZero[(ORCInonZero$lowlim!=1 & ORCInonZero$median!=1) |(ORCInonZero$upplim!=1 & ORCInonZero$median!=1),],rownames(ORCInonZero))
#dev.off()
#write.table(ORCInonZero[nrow(ORCInonZero):1,],file="OR_lasso.txt",col.names=NA,sep="\t")

Smd=smd(ds10[,c(exposed,outcome)],exposed=exposed,variable=outcome,verbose=FALSE,categorical=TRUE)
se_adj=(lassomodCI$betaCIlim[exposed,3]-lassomodCI$betaCIlim[exposed,1])/2/1.96
results=c(n0=nrow(ds10),n1=NA,
          KS=NA,KL=NA,Cstat=NA,SMD=Smd,
          ORadj=ORCInonZero[exposed,2],
          ORlow_adj=ORCInonZero[exposed,1],
          ORupp_adj=ORCInonZero[exposed,3],
          coeff_adj=lassomodCI$betaCIlim[exposed,2],
          se_adj=se_adj,
          rep(NA,5))
lassoMVResults=list(results=results,matchedID=NULL)

end=proc.time()[1:3]
time["lassoMV",,i]=end-start


####### random sampling ##################

print("############### random downsampling 1x ############")

start=proc.time()[1:3]

matchedID_rdm=cbind(ds10[ds10[,exposed]==1,id],sample(ds10[ds10[,exposed]==0,id],nminor,replace=FALSE))
rdm1xResults=extractResultsRdm(matchedID_rdm,ds10,fmod=NULL,id=id,exposed=exposed,outcome=outcome,verbose=FALSE)

end=proc.time()[1:3]
time["rdmSamp1x",,i]=end-start


print("############### random bagging (repeatedly downsample major class Ndraw times WITH replacement) ############")

start=proc.time()[1:3]
for(j in 1:Ndraw){
  matchedID_bag=cbind(ds10[ds10[,exposed]==1,id],sample(ds10[ds10[,exposed]==0,id],nminor,replace=TRUE))
  rdmBag[[j]]=extractResultsRdm(matchedID_bag,ds10,fmod=NULL,id=id,exposed=exposed,outcome=outcome,verbose=FALSE)
}

temp=matrix(unlist(lapply(rdmBag,`[`,"results")),nrow=Ndraw,ncol=length(rdmBag[[1]][[1]]),byrow=T)
colnames(temp)=names(rdmBag[[1]][[1]])

rdmBagResults=list(results=colMeans(temp,na.rm=T),matchedID=NULL)

end=proc.time()[1:3]
time["rdmBag",,i]=end-start



print("############### random jacknife (repeatedly downsample major class Ndraw times WITHOUT replacement) ############")

start=proc.time()[1:3]
for(j in 1:Ndraw){
  #jackknife by shuffling PS scores and then resample
  matchedID_JK=cbind(ds10[ds10[,exposed]==1,id],sample(ds10[ds10[,exposed]==0,id],nminor))
  rdmJK[[j]]=extractResultsRdm(matchedID_JK,ds10,fmod=NULL,id=id,exposed=exposed,outcome=outcome,verbose=FALSE)
}

temp=matrix(unlist(lapply(rdmJK,`[`,"results")),nrow=Ndraw,ncol=length(rdmJK[[1]][[1]]),byrow=T)
colnames(temp)=names(rdmJK[[1]][[1]])

rdmJKResults=list(results=colMeans(temp,na.rm=T),matchedID=NULL)

end=proc.time()[1:3]
time["rdmJK",,i]=end-start


######consolidate matchedIDs
matchedID[[i]]=list(expertResults[[2]],hdResults[[2]],lassoResults[[2]],rfResults[[2]],
                    euclideanResults[[2]],
                    jaccardResults[[2]],diceResults[[2]],cosineResults[[2]],
                    pearsonResults[[2]],spearmanResults[[2]],
                    lassoMVResults[[2]],
                    rdm1xResults[[2]],rdmBagResults[[2]],rdmJKResults[[2]])
names(matchedID[[i]])=dimnames(time)[[1]]


#######check if baseline variables are balanced
beforeMatching=matchStats(numvar=empvariables_num,catvar=empvariables_cat,treatment=exposed,
                          data=ds10,
                          outXlsx=NULL,verbose=FALSE)
afterMatching=list()
for(j in 1:12){
  #lassoMV, rdmBag and rdmJK req special handling (shouldn't be matched)
  if(j==11){
    afterMatching[[j]]=afterMatching[[1]]
  }else{
    afterMatching[[j]]=matchStats(numvar=empvariables_num,catvar=empvariables_cat,treatment=exposed,
                                  data=ds10[matchedID[[i]][[j]],],
                                  outXlsx=NULL,verbose=FALSE)
  }
}
names(afterMatching)=names(matchedID[[i]])[-(13:14)]

#extract pval matrix
pmat=extractPval(beforeMatching) #initialize pmat vector with beforeMatching
pmat=cbind(pmat,sapply(afterMatching,extractPval)) #append pmat vector with afterMatching
colnames(pmat)[1]="before"
pmat=cbind(pmat,bootstrap=NA,jackknife=NA)
pmat[,"lassoMV"]=NA

#extract smd matrix
smdmat=extractSmd(beforeMatching) #initialize pmat vector with beforeMatching
smdmat=cbind(smdmat,sapply(afterMatching,extractSmd)) #append pmat vector with afterMatching
colnames(smdmat)[1]="before"
smdmat=cbind(smdmat,bootstrap=NA,jackknife=NA)
smdmat[,"lassoMV"]=NA


############consolidate results
resultsmat=as.data.frame(rbind(expertResults[[1]],hdResults[[1]],lassoResults[[1]],rfResults[[1]],
                               euclideanResults[[1]],
                               jaccardResults[[1]],diceResults[[1]],cosineResults[[1]],
                               pearsonResults[[1]],spearmanResults[[1]],
                               lassoMVResults[[1]],
                               rdm1xResults[[1]],rdmBagResults[[1]],rdmJKResults[[1]]))
rownames(resultsmat)=dimnames(time)[[1]]
resultsmat$bias_adj=resultsmat$coeff_adj-trueCoeff
resultsmat$bias_matched=resultsmat$coeff_matched-trueCoeff
resultsmat$inCI_adj=(resultsmat$ORlow_adj<=trueOR & resultsmat$ORupp_adj>=trueOR)
resultsmat$inCI_matched=(resultsmat$ORlow_matched<=trueOR & resultsmat$ORupp_matched>=trueOR)
resultsArray[,,i]=as.matrix(resultsmat[,desiredOrder])
wantedVarList[[i]]=wantedvar
var_corOutcomeList[[i]]=var_corOutcome
pmatArray[,,i]=t(pmat[paste("x",1:10,sep=""),])
smdmatArray[,,i]=t(smdmat[paste("x",1:10,sep=""),])


##########save periodidically
  if(i %in% seq(50,900,by=50)){
    save(resultsArray,time,wantedVarList,rfImpt,lassoPSBeta,lassoMVvar,bestlambda,var_corOutcomeList,matchedID,
         Noutcomesmat,Nexposed,trueORvec,trueSMDvec,trueCoeffvec,pmatArray,smdmatArray,file="sim2000x10_1000MC.RData")
  }
}   #end of Nsim Monte Carlo simulations

endAll=proc.time()[1:3]
endAll-startAll

#remember to terminate the cores when done
stopCluster(cl)

dimnames(resultsArray)[1:2]=list(rownames(resultsmat),desiredOrder)
dimnames(pmatArray)[1:2]=list(colnames(pmat),paste("x",1:10,sep=""))
dimnames(smdmatArray)[1:2]=list(colnames(smdmat),paste("x",1:10,sep=""))
save(resultsArray,time,wantedVarList,rfImpt,lassoPSBeta,lassoMVvar,bestlambda,var_corOutcomeList,matchedID,
     Noutcomesmat,Nexposed,trueORvec,trueSMDvec,trueCoeffvec,pmatArray,smdmatArray,file="sim2000x10_1000MC.RData")

round(apply(resultsArray,c(1,2),mean,na.rm=T),3)
apply(time[,"elapsed",],1,summary)
apply(rfImpt,1,summary)
apply(lassoPSBeta,1,summary)

rowSums(resultsArray[,"inCI_matched",],na.rm=T)  #excludes lassoMV
rowSums(resultsArray[,"inCI_adj",],na.rm=T)  #for lassoMV
#round(rowMeans(abs(resultsArray["expertPS",c("n1","KS","KL","Cstat","se_adj","bias_adj","se_matched","bias_matched"),]),na.rm=T),3)
rowMeans(rbind(trueORvec,resultsArray[,"ORadj",]),na.rm=T)
rowMeans(rbind(trueORvec,resultsArray[,"ORmatched",]),na.rm=T)





