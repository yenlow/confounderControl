# To simulated data 2000x100, perform similarity search 
# by Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman correl
# Tune caliper (Euc: mean+3sd, Jac: 0.6, others: 0.7 or 0.8)
# 
# 03-Apr-14 Yen Low
##############################################################################
#
#
# Read in the original patients file 
require(rJava)
require(Matching)
require(Epi) #clogistic
require(glmnet)
require(randomForest)
require(vegan)
require(FNN)

#############SET PARAMETERS######################
#RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
#ProjectPath=Sys.getenv("ProjectPath") #Project is an environment variable, path="/mnt/hgfs/projects"
RScriptPath="~"
ProjectPath="~/projects"
setwd(paste(ProjectPath,"/psm/models/sim2000x100",sep=""))
getwd()

#required functions
source("../../hdpsFunctions.R")
source(paste(RScriptPath,"/scripts/R/logit.R",sep=""))

id="id"
exposed="exposure" #name of exposure variable (must be string)
outcome="outcome"

############tune nsd
desiredOrder=c("n1","ORmatched","ORlow_matched","ORupp_matched","coeff_matched","se_matched","bias_matched","inCI_matched",
               "n0","ORadj","ORlow_adj","ORupp_adj",            "coeff_adj","se_adj",       "bias_adj","inCI_adj",
               "SMD","KS","KL","Cstat")
nsds=rev(c(0,1,1.5,2,2.5,3,10))

Nsim=100
resultsArrayEuc=array(dim=c(length(nsds),10,Nsim))
Nexposed=Noutcomesmat=matrix(nrow=Nsim,ncol=2)
trueORvec=trueSMDvec=trueCoeffvec=c()

for(i in 1:Nsim){
  
  print(i)
  #generate simulated data
  source("../../data/2000by100.R")
  ds$id=1:nrow(ds)
  
  trueORvec[i]=trueOR
  trueSMDvec[i]=trueSmd
  trueCoeffvec[i]=trueCoeff
  Noutcomesmat[i,]=Noutcomes
  Nexposed[i]=sum(ds[,exposed])
  
  xmat_ctrl=xmat[ds[,exposed]==0,]
  xmat_trted=xmat[ds[,exposed]==1,]
  rownames(xmat_ctrl)=ds[ds[,exposed]==0,id]
  rownames(xmat_trted)=ds[ds[,exposed]==1,id]

  resultsmat=c()
  for(q in 1:length(nsds)){
    matchedEuclidean=matchByDist(xmat_ctrl,xmat_trted,method="euclidean",k_neighbors=5,caliper=0.8,nsd=nsds[q])
    euclideanResults=extractResults(ps=matchedEuclidean,exposurevec=NULL,
                                    data=ds,fmod=NULL,id=id,exposed=exposed, outcome=outcome,logitFlag=F,
                                    outfile=NULL,verbose=FALSE)
    
    ############consolidate results
    resultsmat=rbind(resultsmat,euclideanResults)
  }
    colnames(resultsmat)=names(euclideanResults)
    resultsmat=as.data.frame(resultsmat)
    resultsmat$bias_matched=resultsmat$coeff_matched-trueCoeff
    resultsmat$inCI_matched=(resultsmat$ORlow_matched<=trueOR & resultsmat$ORupp_matched>=trueOR)
    resultsArrayEuc[,,i]=as.matrix(resultsmat[,desiredOrder[c(1:9,17)]])  
}   # end of Nsim loop for Euc

dimnames(resultsArrayEuc)[1:2]=list(nsds,desiredOrder[c(1:9,17)])
save(resultsArrayEuc,trueORvec,trueSMDvec,trueCoeffvec,file="tuneNsdEuc.RData")

#get final results matrix
reportmat=
rbind(
  c(n, rowMeans(resultsArrayEuc[,"n1",])),  #1-3 sd (1sd is equiv to 1316 patients total)
  rowMeans(rbind(trueORvec,resultsArrayEuc[,"ORmatched",])),   #1-3 sd
  rowMeans(rbind(trueSMDvec,resultsArrayEuc[,"SMD",])),    #any
  c(0,rowMeans(abs(resultsArrayEuc[,"bias_matched",]))),  #2-3 sd
  c(0,rowMeans(resultsArrayEuc[,"se_matched",])),  #1.5-3 sd
  c(Nsim,rowSums(resultsArrayEuc[,"inCI_matched",]))/Nsim #2-3 sd
)
rownames(reportmat)=c("n1","ORmatched","SMD","bias_matched","se_matched","inCI_matched")
colnames(reportmat)[1]="true_value"
write.table(reportmat,file="tuneNsdEuc.txt",sep="\t",col.names=NA)


#plot final results matrix
png(file="tuneNsdEuc.png",width=7,height=11,units="in",res=300,bg="transparent")
par(mfrow=c(3,2),mar=c(2,3,3,1))
barplot(reportmat[1,-1],space=0.1,ylim=c(0,2000),main="Number of patients",cex.main=1)
abline(h=reportmat[1,1],lwd=5)
barplot(reportmat[2,-1],space=0.1,main="Odds ratio",cex.main=1)
abline(h=reportmat[2,1],lwd=5)
barplot(reportmat[3,-1],space=0.1,main="Effect size (standardized mean difference)",cex.main=1)
abline(h=reportmat[3,1],lwd=5)
barplot(reportmat[4,-1],space=0.1,main="Bias",cex.main=1)
abline(h=0,lwd=5)
barplot(reportmat[5,-1],space=0.1,ylim=c(0,0.4),main="Std error of beta",cex.main=1)
abline(h=0,lwd=5)
barplot(reportmat[6,-1],space=0.1,ylim=c(0,1),main="Coverage of true OR within 95% CI",cex.main=1)
abline(h=1,lwd=5)
dev.off()


############tune caliper for similarity methods
calipers=c(0,0.1,0.2,0.3,0.4,0.5,0.6)
#calipers=c(0,0.6)


Nsim=100
resultsArrayCal=array(dim=c(5,10,length(calipers),Nsim))
Noutcomesmat=matrix(nrow=Nsim,ncol=2)
trueORvec=trueSMDvec=trueCoeffvec=c()

for(i in 1:Nsim){
  
  print(i)
  #generate simulated data
  source("../../data/2000by100.R")
  ds$id=1:nrow(ds)
  
  trueORvec[i]=trueOR
  trueSMDvec[i]=trueSmd
  trueCoeffvec[i]=trueCoeff
  Noutcomesmat[i,]=Noutcomes
  
  xmat_ctrl=xmat[ds[,exposed]==0,]
  xmat_trted=xmat[ds[,exposed]==1,]
  rownames(xmat_ctrl)=ds[ds[,exposed]==0,id]
  rownames(xmat_trted)=ds[ds[,exposed]==1,id]
  
  for(p in 1:length(calipers)){
    matchedJaccard=matchByDist(xmat_ctrl,xmat_trted,method="jaccard",k_neighbors=5,caliper=calipers[p])
    matchedDice=matchByDist(xmat_ctrl,xmat_trted,method="dice",k_neighbors=5,caliper=calipers[p])
    matchedCosine=matchByDist(xmat_ctrl,xmat_trted,method="cosine",k_neighbors=5,caliper=calipers[p])
    matchedPearson=matchByDist(xmat_ctrl,xmat_trted,method="pearson",k_neighbors=5,caliper=calipers[p])
    matchedSpearman=matchByDist(xmat_ctrl,xmat_trted,method="spearman",k_neighbors=5,caliper=calipers[p])
    
    jaccardResults=extractResults(ps=matchedJaccard,exposurevec=NULL,
                                  data=ds,fmod=NULL,id=id,exposed=exposed, outcome=outcome,logitFlag=F,
                                  outfile=NULL,verbose=FALSE)
    diceResults=extractResults(ps=matchedDice,exposurevec=NULL,
                               data=ds,fmod=NULL,id=id,exposed=exposed, outcome=outcome,logitFlag=F,
                               outfile=NULL,verbose=FALSE)
    cosineResults=extractResults(ps=matchedCosine,exposurevec=NULL,
                                 data=ds,fmod=NULL,id=id,exposed=exposed, outcome=outcome,logitFlag=F,
                                 outfile=NULL,verbose=FALSE)
    pearsonResults=extractResults(ps=matchedPearson,exposurevec=NULL,
                                  data=ds,fmod=NULL,id=id,exposed=exposed, outcome=outcome,logitFlag=F,
                                  outfile=NULL,verbose=FALSE)
    spearmanResults=extractResults(ps=matchedSpearman,exposurevec=NULL,
                                   data=ds,fmod=NULL,id=id,exposed=exposed, outcome=outcome,logitFlag=F,
                                   outfile=NULL,verbose=FALSE)
    
    ############consolidate results
    resultsmat=as.data.frame(rbind(jaccardResults,diceResults,cosineResults,pearsonResults,spearmanResults))
    rownames(resultsmat)=c("jaccard","dice","cosine","pearson","spearman")
    resultsmat$bias_matched=sum(resultsmat$coeff_matched,-trueCoeff)
    resultsmat$inCI_matched=(resultsmat$ORlow_matched<=trueOR & resultsmat$ORupp_matched>=trueOR)
    resultsArrayCal[,,p,i]=as.matrix(resultsmat[,desiredOrder[c(1:9,17)]])
    
  }
  for(i in seq(10,90,by=10)) save(resultsArrayCal,trueORvec,trueSMDvec,trueCoeffvec,file="tuneCalSim.RData")
}  #end of Nsim loop for similarities
dimnames(resultsArrayCal)[1:3]=list(rownames(resultsmat),desiredOrder[c(1:9,17)],calipers)
save(resultsArrayCal,trueORvec,trueSMDvec,trueCoeffvec,Noutcomesmat,Nexposed,file="tuneCalSim.RData")

dim(resultsArrayCal)
apply(resultsArrayCal,c(1,2,3),mean,na.rm=T)

sink("tuneCaliperSim.txt")
print("jaccard")
round(apply(resultsArrayCal["jaccard",,,],c(1,2),mean,na.rm=T),3)  #0.6
print("dice")
round(apply(resultsArrayCal["dice",,,],c(1,2),mean,na.rm=T),3)     #0.7 (or 0.8)
print("cosine")
round(apply(resultsArrayCal["cosine",,,],c(1,2),mean,na.rm=T),3)   #0.8 (or 0.7)
print("pearson")
round(apply(resultsArrayCal["pearson",,,],c(1,2),mean,na.rm=T),3)  #0.7 (or 0.8)
print("spearman")
round(apply(resultsArrayCal["spearman",,,],c(1,2),mean,na.rm=T),3) #0.7
sink()

