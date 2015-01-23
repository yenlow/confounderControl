# To PAD data 5757 x 447 vars (dem + 10% sparsity, 412 cids), 
# do expertPS, hdPS, lasso logit, RF,
# Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman correl
# Stitch the results from pad_10sparseCid, pad_10sparseCid2, pad_10sparseCid3
# 
# 04-Sep-14 Yen Low
##############################################################################
#
#

#############SET PARAMETERS######################
#RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
#ProjectPath=Sys.getenv("ProjectPath") #Project is an environment variable, path="/mnt/hgfs/projects"
RScriptPath="~/scripts"
ProjectPath="~/projects"
setwd(paste(ProjectPath,"/psm/models/pad_10sparseCid",sep=""))

#load data
load(file="PAD10sparse_models.RData",verbose=T)
load(file="PAD10sparse_results.RData",verbose=T)
resultsArray_stitched=resultsArray
matchedID_stitched=matchedID

load(file="../pad_10sparseCid2/PAD10sparse_results.RData",verbose=T)
resultsArray_stitched[,,10:17]=resultsArray[,,10:17]
for(i in 10:17) matchedID_stitched=c(matchedID_stitched,list(matchedID[[i]]))

load(file="../pad_10sparseCid3/PAD10sparse_results.RData",verbose=T)
resultsArray_stitched[,,18:20]=resultsArray[,,18:20]
for(i in 18:20) matchedID_stitched=c(matchedID_stitched,list(matchedID[[i]]))

resultsArray=resultsArray_stitched
matchedID=matchedID_stitched
save(resultsArray,matchedID,file="PAD10sparse_results.RData")

