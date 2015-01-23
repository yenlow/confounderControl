# To sim data 2000 x 10 variables, do expertPS, hdPS, lasso logit, RF,
# Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman correl
# 
# 20-Apr-14 Yen Low
##############################################################################
#
#
#############SET PARAMETERS######################
RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
ProjectPath=Sys.getenv("ProjectPath") #Project is an environment variable, path="/mnt/hgfs/projects"
#RScriptPath="~"
#ProjectPath="~/projects"
setwd(paste(ProjectPath,"/psm/models/sim2000x10_compare",sep=""))
getwd()


load(file="../sim2000x10/tuneCalSim.RData")

apply(resultsArrayCal[,"ORmatched",,],c(1,2),mean,na.rm=T)


data=as.data.frame(resultsArray[,grepl("^OR",dimnames(resultsArray)[[2]]),1])
plotOR(data,filename="PAD_10sparseCid_forest.png",title="411 cid features (> 10% sparsity)",width=5,height=5,bg="white")



plotOR(data)

pdf(file="PAD_10sparseCid_forest.pdf",width=8.5,height=11)
par(mfrow=c(3,2))
for(i in 1:10){
  data=as.data.frame(resultsArray[,grepl("^OR",dimnames(resultsArray)[[2]]),i])
  plotOR(data,title=dimnames(resultsArray)[[3]][i],legend=FALSE)
}
dev.off()
