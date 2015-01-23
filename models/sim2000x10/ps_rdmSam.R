# To simulated data 2000x10, do random sampling repeatedly
# 
# 29-Apr-14 Yen Low
##############################################################################


#############SET PARAMETERS######################
#RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
#ProjectPath=Sys.getenv("ProjectPath") #Project is an environment variable, path="/mnt/hgfs/projects"
RScriptPath="~"
ProjectPath="~/projects"
setwd(paste(ProjectPath,"/psm/models/sim2000x10",sep=""))
getwd()

#required functions
source("../../hdpsFunctions.R")
source(paste(RScriptPath,"/scripts/R/logit.R",sep=""))

id="id"
exposed="exposure" #name of exposure variable (must be string)
outcome="outcome"
desiredOrder=c("n1","ORmatched","ORlow_matched","ORupp_matched","coeff_matched","se_matched","bias_matched","inCI_matched",
               "n0","ORadj","ORlow_adj","ORupp_adj",            "coeff_adj","se_adj",       "bias_adj","inCI_adj",
               "SMD","KS","KL","Cstat")

#to store consolidated results
startAll=proc.time()[1:3]

Nsim=1000
Ndraw=100
Nmethods=3
resultsArray=array(dim=c(Nmethods,20,Nsim))  #3methods (incl rdm downsampling), 20 metrics
drawArray=array(dim=c(Nmethods,16,Ndraw))  #3methods (incl rdm downsampling), 20-4 metrics
Noutcomesmat=matrix(nrow=Nsim,ncol=2)
time=array(dim=c(Nmethods,3,Nsim))
dimnames(time)[[1]]=c("rdm1x","rdmBag","rdmJK")
dimnames(time)[[2]]=c("user","system","elapsed")
dimnames(resultsArray)[[1]]=dimnames(time)[[1]]
dimnames(drawArray)[[1]]=dimnames(time)[[1]]
dimnames(drawArray)[[2]]=c("n0","n1","KS","KL","Cstat","SMD",
                            "ORadj","ORlow_adj","ORupp_adj","coeff_adj","se_adj",
                            "ORmatched","ORlow_matched","ORupp_matched","coeff_matched","se_matched")
Nexposed=trueORvec=trueSMDvec=trueCoeffvec=bestlambda=c()

for(i in 1:Nsim){
	
cat("\n\n\n###############",i,"###########\n\n")

#generate simulated data
source("../../data/2000by10.R")
#load("../../data/2000x10.RData",verbose=T)
#rm(list=ls()[!(ls() %in% "ds10")]) #remove all objects except ds10
ds10$id=1:nrow(ds10)

trueORvec[i]=trueOR
trueSMDvec[i]=trueSmd[2]
trueCoeffvec[i]=trueCoeff
Nexposed[i]=sum(ds10[,exposed])
Noutcomesmat[i,]=Noutcomes


nminor=min(table(ds10[,exposed]))

print("############### random downsampling 1x ############")

start=proc.time()[1:3]

matchedID=cbind(ds10[ds10[,exposed]==1,id],sample(ds10[ds10[,exposed]==0,id],nminor,replace=FALSE))
rdm1xResults=extractResultsRdm(matchedID,ds10,fmod=NULL,id=id,exposed=exposed,outcome=outcome,verbose=FALSE)
drawArray[1,,1]=rdm1xResults[[1]]

end=proc.time()[1:3]
time[1,,i]=end-start



print("############### random bagging (repeatedly downsample major class Ndraw times WITH replacement) ############")

start=proc.time()[1:3]
for(j in 1:Ndraw){
  matchedID=cbind(ds10[ds10[,exposed]==1,id],sample(ds10[ds10[,exposed]==0,id],nminor,replace=TRUE))
  rdmBagResults=extractResultsRdm(matchedID,ds10,fmod=NULL,id=id,exposed=exposed,outcome=outcome,verbose=FALSE)
  drawArray[2,,j]=rdmBagResults[[1]]
}

end=proc.time()[1:3]
time[2,,i]=end-start



print("############### random jacknife (repeatedly downsample major class Ndraw times WITHOUT replacement) ############")

start=proc.time()[1:3]
for(j in 1:Ndraw){
#jackknife by shuffling PS scores and then resample
  matchedID=cbind(ds10[ds10[,exposed]==1,id],sample(ds10[ds10[,exposed]==0,id],nminor))
  rdmJKResults=extractResultsRdm(matchedID,ds10,fmod=NULL,id=id,exposed=exposed,outcome=outcome,verbose=FALSE)
  drawArray[3,,j]=rdmJKResults[[1]]
}

end=proc.time()[1:3]
time[3,,i]=end-start



######can't consolidate matchedIDs
############consolidate results
resultsmat=as.data.frame(apply(drawArray,c(1,2),mean,na.rm=T))
resultsmat$bias_adj=resultsmat$coeff_adj-trueCoeff
resultsmat$bias_matched=resultsmat$coeff_matched-trueCoeff
resultsmat$inCI_adj=(resultsmat$ORlow_adj<=trueOR & resultsmat$ORupp_adj>=trueOR)
resultsmat$inCI_matched=(resultsmat$ORlow_matched<=trueOR & resultsmat$ORupp_matched>=trueOR)
resultsArray[,,i]=as.matrix(resultsmat[,desiredOrder])

if(i %in% seq(50,900,by=50)) save(resultsArray,time,Noutcomesmat,Nexposed,trueORvec,trueSMDvec,trueCoeffvec,file="sim2000x10_1000MC_rdm.RData")

}   #end of Nsim Monte Carlo simulations

endAll=proc.time()[1:3]
endAll-startAll


dimnames(resultsArray)[1:2]=list(rownames(resultsmat),desiredOrder)
save(resultsArray,time,Noutcomesmat,Nexposed,
     trueORvec,trueSMDvec,trueCoeffvec,bestlambda,file="sim2000x10_1000MC_rdm.RData")

round(apply(resultsArray,c(1,2),mean,na.rm=T),3)
apply(time[,"elapsed",],1,summary)

rowSums(resultsArray[,"inCI_matched",],na.rm=T)
rowMeans(rbind(trueORvec,resultsArray[,"ORmatched",]),na.rm=T)





