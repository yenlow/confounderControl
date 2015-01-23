# To PAD data 5757 x 447 vars (dem + 10% sparsity, 412 cids), 
# do expertPS, hdPS, lasso logit, RF,
# Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman correl
# Stitch the noagesex enforced hdpsobj from pad_10sparseCid, pad_10sparseCid2, pad_10sparseCid3
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
load(file="hdpsobj.RData",verbose=T)
hdpsobj_stitch=hdpsobj

load(file="../pad_10sparseCid2/hdpsobj.RData",verbose=T)
for(i in 10:17) hdpsobj_stitch=c(hdpsobj_stitch,list(hdpsobj[[i]]))

load(file="../pad_10sparseCid3/hdpsobj.RData",verbose=T)
for(i in 18:20) hdpsobj_stitch=c(hdpsobj_stitch,list(hdpsobj[[i]]))

hdpsobj=hdpsobj_stitch
save(hdpsobj_stitch,file="hdpsobj_stitched.RData")

