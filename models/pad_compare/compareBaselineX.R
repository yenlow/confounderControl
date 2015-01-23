# To PAD data of various sparsities, compare how well baseline characteristics match
# after expertPS, hdPS, lasso logit, RF,
# Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman correl
# flip colors of bubble plot (blue=protective, red=risky)
# extract pval and smd
# 
# 11-Sep-14 moved rdmSamp1x to last row, lassoPS has poor results
# 27-Aug-14 Yen Low
##############################################################################
#
#
#############SET PARAMETERS######################
#RScriptPath=Sys.getenv("RScriptPath") #Rscriptpath is an environment variable, path="/mnt/hgfs/Dropbox"
#ProjectPath=Sys.getenv("ProjectPath") #Project is an environment variable, path="/mnt/hgfs/projects"
RScriptPath="~/scripts"
ProjectPath="~/projects"
setwd(paste(ProjectPath,"/psm/models/pad_compare",sep=""))
getwd()

require(ggplot2)
require(reshape2)
require("grid")

source("../../charts.R")
source("../../hdpsFunctions.R")
source(paste(RScriptPath,"/R/utils.R",sep=""))
source(paste(RScriptPath,"/R/patientChar.R",sep=""))
source(paste(RScriptPath,"/R/logit.R",sep=""))

modelOrder=c("expertPS","hdPS","lassoPS","rfPS",
             "euclidean","jaccard","dice","cosine","pearson","spearman",
             "lassoMV","rdmBag","rdmJK","rdmSamp1x")

#load results from running various methods on PAD_10sparseCid
load(file="../pad_10sparseCid/PAD10sparse_results.RData",verbose=TRUE)

### bubble plot for OR ###
#scale_color_gradient2(limits=c(-2.5,2.5),low="red",mid="white",high="blue",midpoint=0)
#colorRampPalette(c(red,white))
colscale=colorpanel(11, "blue", "white", "red")

#prepare data for bubble plot
data=melt(resultsArray[,c("coeff_matched","se_matched","coeff_adj","se_adj"),])
data=dcast(data,X1+X3~X2)
colnames(data)[1:2]=c("method","outcome")
data$outcomelevel=factor(gsub(".after","",data$outcome),gsub(".after","",dimnames(resultsArray)[[3]]))
data$methodlevel=factor(data$method,rev(modelOrder))
#replace quasi-sep model (rdmmatched)
data[which(data$se_matched>3),c("coeff_matched","se_matched")]=NA

#overlay tigher adj on top of matched values
png("PAD_10sparseCid_bubble.png",width=11,height=7,units="in",res=300,bg="transparent")
p = ggplot(data,aes(x=outcomelevel,y=methodlevel))
p = p + theme_bw()
p = p + theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1))
#p = p + scale_color_gradientn(values=seq(-2.5,2.5,length=11),colours=colscale,rescaler = function(x, ...) x) + scale_size_area(max_size=15)
p = p + scale_color_gradient2(limits=c(-2,2),low="blue",mid="white",high="red",midpoint=0) + scale_size_area(max_size=15)
p = p + geom_point(aes(colour=coeff_matched,size=se_matched),alpha=1,na.rm=T) 
p = p + geom_point(aes(colour=coeff_adj,size=se_adj),alpha=1,na.rm=T)
plot(p)
dev.off()


##### forest plot
pdf(file="PAD_10sparseCid_forest.pdf",width=11,height=8.5)
par(mfrow=c(2,3))
for(i in 1:20){
#for(i in c(1,8,10)){
  data=as.data.frame(resultsArray[,grepl("^OR",dimnames(resultsArray)[[2]]),i])
  data[is.infinite(as.matrix(data))]=NA
  plotOR(data,title=dimnames(resultsArray)[[3]][i],legend=FALSE)
}
dev.off()


#### baseline characteristics
#load data
load("../../data/pad/5757x412_10sparseCid.RData")
load("../../data/pad/anna/PAD_PLOSOnePatientData.RData")

id="pid_org"
exposed="cilostazol" #name of exposure variable (must be string)
outcome=c(  "mace.after","cardiac.arrest.after","death","ec.after","mi.after","stroke.after","sudden.cardiac.death.after",             
            "male.after","amputation.after","angioplasty.after","bypass.after","revasc.after",
            "arrythmias.after","af.after","bradycardia.after","tachycardia.after","vf.after","vt.after",
            "dizziness.after","palpitations.after")
demvariables=c("min_pad_age_scaled","gender.male","gender.unknown","race.white","race.asian","race.black","race.unknown")
expertvariables=c(	"chf.before", "dyslipidemias.before", "renal.failure.before",
                   "aspirin.before", "clopidogrel.before", "warfarin.before","antiarrhythmics.before",
                   "mace.before", "male.before",
                   "hypertension.before","statins.before", "ace.inhibitors.before","arrythmias.before")
empvariables=colnames(rawCidDatamat_rmCol)[-1] #except pid_org
empvariables_cat=empvariables[-1]   #except min_pad_age_scaled
empvariables_num="min_pad_age_scaled"

#formula for expertPS model: exposure ~ var1 + var 2 + ...
fPS = as.formula(paste(exposed," ~ ",paste(unique(c(demvariables,expertvariables)),collapse=" + "),sep=""))
fPS

#remove bias and inCI (not applicable in this case since no gold std)
desiredOrder=c("n1","ORmatched","ORlow_matched","ORupp_matched","coeff_matched","se_matched",
               "n0","ORadj","ORlow_adj","ORupp_adj",            "coeff_adj","se_adj",
               "SMD","KS","KL","Cstat")

#put data together
#ds=merge(variables.data.final[,c(id,exposed,outcome,"min_pad_age",expertvariables)],rawCidDatamat_rmCol,by=id,all.y=T,sort=T)

#remove unnecesssary objects to save memory
rm(rawCidDatamat_rmCol,cil.control.m5)
lsos()
gci()
gc()

#expertvariables
priorDrugs=c("statins.before","ace.inhibitors.before","aspirin.before","clopidogrel.before","warfarin.before","antiarrhythmics.before")
priorDiseases=c("mace.before","chf.before","arrythmias.before","dyslipidemias.before","hypertension.before","renal.failure.before","male.before")


#generate spreadsheet tables of baseline char (hdps depends on MACE outcome; the other matchingID are outcome-independent)
beforeMatching=matchStats(numvar="min_pad_age",catvar=c("gender","race",priorDrugs,priorDiseases),treatment=exposed,
                          data=variables.data.final,
                          outXlsx="beforeMatching.xlsx",verbose=FALSE)
methodNames=names(matchedID[[20]])[-(13:14)]  #can't match rdmBag and rdmJK
afterMatching=list()
for(i in setdiff(c(1:12),11)){  #no matched metrics for lassoMV
  outXlsx=paste("after",methodNames[i],".xlsx",sep="")
  afterMatching[[i]]=matchStats(numvar="min_pad_age",catvar=c("gender","race",priorDrugs,priorDiseases),treatment=exposed,
                                              data=variables.data.final[variables.data.final[,id] %in% matchedID[[20]][[i]],],
                                              outXlsx=outXlsx,verbose=FALSE)
}
afterMatching[[11]]=afterMatching[[1]]   #set empty lassoMV to expertMV to avoid NA (set NA later)
names(afterMatching)=methodNames


#extract pval matrix
pmat=extractPval(beforeMatching) #initialize pmat vector with beforeMatching
pmat=cbind(pmat,sapply(afterMatching,extractPval)) #append pmat vector with afterMatching
colnames(pmat)[1]="before"
pmat=cbind(pmat,rdmBag=NA,rdmJK=NA)
pmat[,"lassoMV"]=NA
pmat=pmat[,c("before",modelOrder)]

#extract pval matrix
smdmat=extractSmd(beforeMatching) #initialize pmat vector with beforeMatching
smdmat=cbind(smdmat,sapply(afterMatching,extractSmd)) #append pmat vector with afterMatching
colnames(smdmat)[1]="before"
smdmat=cbind(smdmat,rdmBag=NA,rdmJK=NA)
smdmat[,"lassoMV"]=NA
smdmat=smdmat[,c("before",modelOrder)]

## spaghetti plot
colors()
col=c("black","gray","red","orange","magenta",
      "lightgreen","blue","olivedrab","cyan","green",
      "darkgreen","maroon","pink")
png("PAD_10sparseCid_spaghetti.png",width=11,height=7,units="in",res=300,bg="white")
par(mar=c(5,6,1,2),xpd=F)
plot(pmat[,1],type="l",ylim=c(10^-10,1),log="y",
     col=col[1],lwd=3,xaxt="n",yaxt="n",xlab="",ylab="")
axis(1,1:nrow(pmat),rownames(pmat),las=2)
yticks=c(0.05,10^c(-10,-5,-3,-1,0))
axis(4,yticks,yticks,las=3,cex=0.5)
abline(h=yticks,col="gray",lwd=0.5)
abline(h=.05,lwd=5,lty="dotted")
#lines(pmat[,2],col=col[2],lwd=3)
#for(i in 3:5) lines(pmat[,i],col=col[i],lwd=3)
#for(i in 6:11) lines(pmat[,i],col=col[i],lwd=3)
#lines(pmat[,13],col=col[13],lwd=3)
for(i in ncol(pmat):2) lines(pmat[,i],col=col[i],lwd=3)
par(xpd=T)
legend("bottomleft",colnames(pmat)[-12],col=col[-12],lwd=3,
       cex=0.8,inset=c(-0.13,0.1),horiz=F)
dev.off()


#heatmaps
colscale=colorpanel(10, "white", "black")

png("pad_heatmap_pval.png",width=6.5,height=7,units="in",res=150,bg="transparent")
heatmap.2(t(-log10(pmat)),Rowv=F,Colv=F,col=colscale,scale="none",
          trace="none",keysize=1,density.info="none",na.rm=T,na.col="transparent")
#heatmap.2(t(1-pmat),Rowv=F,Colv=F,col=colscale,scale="none",
#          trace="none",keysize=1,density.info="none")
dev.off()

png("pad_heatmap_smd.png",width=6.5,height=7,units="in",res=150,bg="transparent")
heatmap.2(t(abs(smdmat)),Rowv=F,Colv=F,col=colscale,scale="none",
          trace="none",keysize=1,density.info="none")
dev.off()

save(pmat,smdmat,afterMatching,beforeMatching,file="baselineMatching.RData")


###### check models rfPsMod, lassoPsmod, lassoMVmod, lassoMVvar
load(file="../pad_10sparseCid/PAD10sparse_models.RData",verbose=T)

demvariables=c("min_pad_age_scaled","gender.male","race.white","race.asian","race.black","race.unknown")
expertvariables=c(  "chf.before", "dyslipidemias.before", "renal.failure.before",
                    "aspirin.before", "clopidogrel.before", "warfarin.before","antiarrhythmics.before",
                    "mace.before", "male.before",
                    "hypertension.before","statins.before", "ace.inhibitors.before","arrythmias.before")
empvariables=names(rfImpt)
empvariables_cat=setdiff(empvariables,c("min_pad_age_scaled"))   #except min_pad_age_scaled

#check rfPS
rfPScoeffInterested=c(rfImpt[rfImpt>0],rfImpt[c(demvariables,expertvariables)])
rfPScoeffInterested=rfPScoeffInterested[!duplicated(names(rfPScoeffInterested))]

#check lassoPS
lassoPScoeff=as.vector(coef(lassoPsmod))
names(lassoPScoeff)=rownames(coef(lassoPsmod))
lassoPScoeffInterested=c(lassoPScoeff[lassoPScoeff!=0],lassoPScoeff[c(demvariables,expertvariables)])
lassoPScoeffInterested=lassoPScoeffInterested[!duplicated(names(lassoPScoeffInterested))]

#check lassoMVmod
lassoMVcoeff=as.vector(coef(lassoMvmod))
names(lassoMVcoeff)=rownames(coef(lassoMvmod))
lassoMVcoeffInterested=c(lassoMVcoeff[lassoMVcoeff!=0],lassoMVcoeff[c(demvariables,expertvariables)])
lassoMVcoeffInterested=lassoMVcoeffInterested[!duplicated(names(lassoMVcoeffInterested))]

#number of times variables were used in hdPS model
hdVarID=gsub("D|V001Once","",wantedvar)
mode(hdVarID)="numeric"
timesSel_hdps=sort(c("min_pad_age_scaled",empvariables_cat[unlist(hdVarID)]),dec=T)

sink("hdPSvarForPalpitations.txt")
timesSel_hdps
sink()

pdf(file="PAD_10sparseCid_varImpt.pdf",width=11,height=8.5)
par(mfrow=c(2,2),mar=c(3,5,2,1))
plot(0:1,0:1,type="n",axes=F,main="hdPS",xlab="",ylab="")
legend("topleft",timesSel_hdps[1:34],pch="",bty="n",cex=0.4)
legend("top",timesSel_hdps[35:68],pch="",bty="n",cex=0.4)
legend("topright",timesSel_hdps[69:101],pch="",bty="n",cex=0.4)
dotchart(sort(lassoPScoeffInterested[!(names(lassoPScoeffInterested) %in% "(Intercept)")]),
         pch=20,main="lassoPS",xlab="beta coefficient")
abline(v=0)
dotchart(sort(rfPScoeffInterested),main="rfPS",xlab="variable importance",pch=20)
abline(v=0)
dotchart(sort(lassoMVcoeffInterested[!(names(lassoMVcoeffInterested) %in% "(Intercept)")]),
        pch=20,main="lassoMV",xlab="beta coefficient")
abline(v=0)
dev.off()


#### output results to suppl table
desiredCol=c("n1","coeff_matched","se_matched","coeff_adj","se_adj")
sink("PAD_10sparseCid_betaSE.txt")
resultsArray[,desiredCol,]
sink()
