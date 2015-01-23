# To PAD data 5757 x 447 vars (dem + 10% sparsity, 412 cids), 
# do expertPS, hdPS, lasso logit, RF,
# Euclidean dist, Jaccard, Dice, Cosine similarities, Pearson and Spearman correl
# force age, sex into model where possible (expertPS, hdPS, lassoPS, lassoMV)
# 
# 04-Sep-14 Yen Low
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
require(Matrix)

#parallel core computing for lasso
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
setwd(paste(ProjectPath,"/psm/models/pad_10sparseCid",sep=""))
getwd()

baseDir = getwd();

# initialize the Java subsystem.  include a jdbc object
#set JVM max heap size to to 12GB (set 12GB to get 60% of it, i.e. ~7GB)
#use pharmacoepi_juan.jar for more than 100 variables
#.jinit(classpath="/home/yenlow/software/pharmacoepi_juan.jar", force.init = TRUE, parameters="-Xmx12g")
.jinit(classpath="/home/yenlow/software/pharmacoepi_juan.jar", force.init = TRUE, parameters="-Xmx20g")
.jclassPath()

#required functions
source("../../hdpsFunctions.R")
source(paste(RScriptPath,"/R/logit.R",sep=""))

#load data
load("../../data/pad/5757x412_10sparseCid.RData",verbose=T)
load("../../data/pad/anna/PAD_PLOSOnePatientData.RData",verbose=T)

id="pid_org"
exposed="cilostazol" #name of exposure variable (must be string)
outcome=c(	"mace.after","cardiac.arrest.after","death","ec.after","mi.after","stroke.after","sudden.cardiac.death.after",             
		"male.after","amputation.after","angioplasty.after","bypass.after","revasc.after",
		"arrythmias.after","af.after","bradycardia.after","tachycardia.after","vf.after","vt.after",
		"dizziness.after","palpitations.after")
#remove "gender.unknown" as there are only 16 such patients
demvariables=c("min_pad_age_scaled","gender.male","race.white","race.asian","race.black","race.unknown")
expertvariables=c(	"chf.before", "dyslipidemias.before", "renal.failure.before",
		"aspirin.before", "clopidogrel.before", "warfarin.before","antiarrhythmics.before",
		"mace.before", "male.before",
		"hypertension.before","statins.before", "ace.inhibitors.before","arrythmias.before")
#outcomemodelvariables_exclmatchingvar=c("x8","x9","x10") 
#447 empirical variables may be matched depending on HD-PS
empvariables=c(demvariables,
               colnames(variables.data.final)[grep(".before$",colnames(variables.data.final))],
               colnames(rawCidDatamat_rmCol)[grep("^[0-9]+$",colnames(rawCidDatamat_rmCol))])
empvariables_cat=setdiff(empvariables,c("min_pad_age_scaled"))   #except min_pad_age_scaled
empvariables_num="min_pad_age_scaled"

#formula for expertPS model: exposure ~ var1 + var 2 + ...
fPS = as.formula(paste(exposed," ~ ",paste(unique(c(demvariables,expertvariables)),collapse=" + "),sep=""))

#formula for outcome after expertPS
#fmod = paste("outcome ~ exposure +", paste("x",c(8:10),collapse=" + ",sep=""))

fPS
#fmod

#remove bias and inCI (not applicable in this case since no gold std)
desiredOrder=c("n1","ORmatched","ORlow_matched","ORupp_matched","coeff_matched","se_matched",
               "n0","ORadj","ORlow_adj","ORupp_adj",            "coeff_adj","se_adj",
               "SMD","KS","KL","Cstat")

#put data together
ds=merge(variables.data.final,
         rawCidDatamat_rmCol[,c(id,setdiff(colnames(rawCidDatamat_rmCol),colnames(variables.data.final)))],by=id,all.y=T,sort=T)
ds=ds[,c(id,exposed,outcome,empvariables)]
dim(ds)   #5757 469
colnames(ds)

#remove unnecesssary objects to save memory
rm(rawCidDatamat_rmCol,cil.control.m5)
lsos()

Noutcome=length(outcome)
#Noutcome=2
Ndraw=100
Nmethods=14 #14 methods (incl. 3 rdm sampling),
nminor=min(table(ds[,exposed]))
resultsArray=array(dim=c(Nmethods,16,Noutcome))   #20-4 metrics
#Noutcomesmat=matrix(nrow=Noutcome,ncol=2)
Nexposed=trueORvec=trueSMDvec=trueCoeffvec=NA
lassoMVvar=matchedID=hdpsobj=rdmBag=rdmJK=list()

#### get PS ###########
###get expertPS values previously calculated by Anna
ps_exp=variables.data.final$psm.score[order(variables.data.final$pid_org)]
names(ps_exp)=variables.data.final$pid_org[order(variables.data.final$pid_org)]


### hdPS values depend on outcome; do within outcomes loop
#generate binary x matrix for hdps
dimdata=ds[,c(id,empvariables_cat)]
dim(dimdata)   #5757  447

###get lassoPS values
xmat=as.matrix(ds[,empvariables])
#set penalty to 0 to force age, sex into model; 1 otherwise
penalty.factor=rep(1,ncol(xmat))
#penalty.factor[colnames(xmat) %in% c("min_pad_age_scaled","gender.male")]=0

#tune lambda for glmnet using 5-fold CV
lassoPsmod=cv.glmnet(xmat,ds[,exposed],alpha=1,family="binomial",standardize=F,nfold=5,penalty.factor=penalty.factor,parallel=TRUE)
plot(lassoPsmod)
psl=as.vector(predict(lassoPsmod,xmat,s=lassoPsmod$lambda.1se,type="response"))
names(psl)=ds[,id]


###get rfPS values
#if regression
#rfPsMod=randomForest(xmat,ds[,exposed],ntree=100,importance=F)
#ps_rf=predict(rfPsMod)
#hist(ps_rf)

#if clasification
rfPsMod=randomForest(xmat,as.factor(ds[,exposed]),ntree=100,importance=T)
ps_rf=predict(rfPsMod, type="prob")[,2]
range(ps_rf,na.rm=T)
#scaling is done in extractResults
#ps_rf=scales::rescale(ps_rf,to=c(max(0.001,min(ps_rf,na.rm=T)),min(0.999,max(ps_rf,na.rm=T))))   #rescale to avoid zeros and 1 which become Inf after logit
hist(ps_rf)
names(ps_rf)=ds[,id]
rfImpt=rfPsMod$importance[,1]

### get random downsampling
matchedID_rdm1x=cbind(ds[ds[,exposed]==1,id],sample(ds[ds[,exposed]==0,id],nminor,replace=FALSE))

###get similarity scores
xmat_ctrl=xmat[ds[,exposed]==0,]
xmat_trted=xmat[ds[,exposed]==1,]
rownames(xmat_ctrl)=ds[ds[,exposed]==0,id]
rownames(xmat_trted)=ds[ds[,exposed]==1,id]

runSimPS<-function(method="jaccard",caliper=0.7,nsd=3,algorithm="kd_tree",i){
  matchedSim=matchByDist(xmat_ctrl,xmat_trted,method=method,k_neighbors=5,caliper=caliper,nsd=nsd,algorithm=algorithm)
  simResults=extractResults(ps=matchedSim,exposurevec=NULL,
                            data=ds,fmod=NULL,id=id,exposed=exposed, outcome=outcome[i],logitFlag=F,
                            outfile=NULL,verbose=FALSE)
  return(simResults)
}

#create sparse matrix for lassoMV
smat=Matrix(cbind(ds[,exposed],xmat),sparse=T)
colnames(smat)[1]=exposed
#set penalty.factor=0 to force cilostazol, age, sex into model; 1 otherwise
penalty.factorMV=rep(1,ncol(smat))
penalty.factorMV[colnames(smat) %in% c("cilostazol")]=0

#save models
save(rfPsMod,rfImpt,lassoPsmod,
     ps_exp,psl,ps_rf,
     file="PAD10sparse_models.RData")
#load(file="PAD10sparse_models.RData")

#### Loop through 20 outcomes

#for(i in 1:Noutcome){
for(i in c(1:9)){

cat("-------", i, outcome[i],"-------\n")

if(i>1) file.remove(paste(tempDir,"output_cohort.txt",sep="/"))
###### expertPSM ######
#scale continuous variables (already scaled to standard normal)

#calculate estimated PS
print("############### expertPS ############")

#start=proc.time()[1:3]

#expertPsMod=glm(fPS, data=ds[,c(id,exposed,demvariables,expertvariables)], family="binomial")
#names(expertPsMod$fitted.values)=ds[,id]
#expertResults=extractResults(ps=expertPsMod$fitted.values,exposurevec=expertPsMod$y,fmod=NULL,
#							data=ds,id=id,exposed=exposed, outcome=outcome[2],logitFlag=TRUE,
#							outfile="PAD10sparse_expertPS.pdf", verbose=TRUE)
#round(expertResults,3)

	expertResults=extractResults(ps=ps_exp,exposurevec=ds[,exposed],fmod=NULL,
								data=ds,id=id,exposed=exposed, outcome=outcome[i],logitFlag=TRUE,
								outfile="PAD10sparse_expertPS_anna.pdf", verbose=TRUE ,printvalues=FALSE,ylim=NULL)

#end=proc.time()[1:3]
#time["expertPS",,i]=end-start

####### HDPS ######################
print("############### hdPS ############")

tempDir = paste(baseDir,"tmp",i,sep="/") #output directory
patientfile=paste(tempDir,"patients.txt",sep="")
dir.create(tempDir)

#prepare data into format required by hdps
#read data file as single string (to use addPatientsFromBuffer)
patientheader=paste(id,exposed,outcome[i],sep="\t")
patientstring=paste(paste(ds[,id],ds[,exposed],ds[,outcome[i]],sep="\t"),collapse="\n")
datainstring=paste(patientheader, patientstring, sep="\n")

### INVOKE pharmacoepi.jar ### (handles cat variables only)
hdpsobj[[i]]=hdps(datainstring,dimdata,outDir=tempDir,Nmostfreq=200,k=200,stratifyDim=FALSE,
             outfile="output_cohort.txt",FullOutput=FALSE,verbose=T,ZeroCellCorrection=F)
#hdpsobj$selectedvariables #(see Fig2, Schneeweiss 2008 for recoded definitions)
save(hdpsobj,file="hdpsobj.RData")
#load(file="hdpsobj.RData",verbose=T)
wantedvar=hdpsobj[[i]]$selectedvariables[grep("1Once$",hdpsobj[[i]]$selectedvariables)]

var_corExposed=empvariables_num[abs(cor(ds[,exposed],ds[,empvariables_num]))>0.05]
var_corOutcome=empvariables_num[abs(cor(ds[,outcome[i]],ds[,empvariables_num]))>0.05] #include in PS
IV=setdiff(var_corExposed,var_corOutcome) #exclude from PS

# Estimate the PS (force numerical variables into PS model)
dataPS=cbind(ds[,c(id,exposed,var_corOutcome)],hdpsobj[[i]]$hdpsdata[,wantedvar])
hdPsMod=glm(paste(exposed,"~ . -",id), data=dataPS, family="binomial")
names(hdPsMod$fitted.values)=as.character(dataPS[,id])
summary(hdPsMod$fitted.values)

hdResults=extractResults(ps=hdPsMod$fitted,exposurevec=hdPsMod$y,fmod=NULL,
                      		data=ds,id=id,exposed=exposed, outcome=outcome[i],logitFlag=TRUE,
                      		outfile="PAD10sparse_hdPS.pdf",verbose=TRUE,printvalues=FALSE,ylim=NULL)


####### lasso ##################

print("############### lassoPS ############")

lassoResults=extractResults(ps=psl,exposurevec=ds[,exposed],fmod=NULL,
                        		data=ds,id=id,exposed=exposed, outcome=outcome[i],logitFlag=TRUE,
                        		outfile="PAD10sparse_lassoPS.pdf",verbose=T,printvalues=FALSE,ylim=90)

####### RF ##################
print("############### rfPS ############")

tryobj=try(extractResults(ps=ps_rf,exposurevec=ds[,exposed],
                         data=ds,fmod=NULL,id=id,exposed=exposed, outcome=outcome[i],logitFlag=TRUE,
                         outfile="PAD10sparse_rfPS.pdf",verbose=TRUE,printvalues=FALSE,ylim=5))
if(class(tryobj)!="try-error") rfResults=tryobj else rfResults=list(NULL,NULL)

####### Similarity ##################
print("############### similarity ############")

euclideanResults=runSimPS(method="euclidean",nsd=3,algorithm="brute",i=i)
jaccardResults=runSimPS(method="jaccard",caliper=0.1,i=i)
diceResults=runSimPS(method="dice",caliper=0.1,i=i)
cosineResults=runSimPS(method="cosine",caliper=0.1,i=i)
pearsonResults=runSimPS(method="pearson",caliper=0.1,i=i)
spearmanResults=runSimPS(method="spearman",caliper=0.1,i=i)

####### lasso multivariate model (NOT PS) ##################
print("############### lasso MV model ############")

#tune lambda for glmnet using 5-fold CV
lassoMvmod=cv.glmnet(smat,ds[,outcome[i]],alpha=1,family="binomial",standardize=F,nfold=5,parallel=T,penalty.factor=penalty.factorMV)
plot(lassoMvmod)
coeff=na.omit(sort(coeffAtlambda(lassoMvmod)))

#bootstrap, repeat 1000 times for 95%CI
#use global lambda (lambda.1se from lassomod.cv)
lassomodbootstrap=genCI(smat,ds[,outcome[i]],ntrials=100,lambda=lassoMvmod$lambda.1se,alpha=1,penalty.factor=penalty.factorMV)
lassomodCI=setCL(lassomodbootstrap)

#convert beta coeff to OR
lassoMVvar[[i]]=lassomodCI[[1]]
ORCInonZero=as.data.frame(exp(lassomodCI$betaCIlim[lassomodCI$beta_nonZero,]))
#png("lassoOR.png",width=7,height=11,units="in",res=300,bg="transparent")
#forestplot(ORCInonZero[(ORCInonZero$lowlim!=1 & ORCInonZero$median!=1) |(ORCInonZero$upplim!=1 & ORCInonZero$median!=1),],rownames(ORCInonZero))
#dev.off()
#write.table(ORCInonZero[nrow(ORCInonZero):1,],file="OR_lasso.txt",col.names=NA,sep="\t")

Smd=smd(ds[,c(exposed,outcome[i])],exposed=exposed,variable=outcome[i],verbose=FALSE,categorical=TRUE)
se_adj=(lassomodCI$betaCIlim[exposed,3]-lassomodCI$betaCIlim[exposed,1])/2/1.96
results=c(n0=nrow(smat),n1=NA,
          KS=NA,KL=NA,Cstat=NA,SMD=Smd,
          ORadj=ORCInonZero[exposed,2],
          ORlow_adj=ORCInonZero[exposed,1],
          ORupp_adj=ORCInonZero[exposed,3],
          coeff_adj=lassomodCI$betaCIlim[exposed,2],
          se_adj=se_adj,
          rep(NA,5))
lassoMVResults=list(results=results,matchedID=NULL)


####### random sampling ensemble ##################
print("############### random downsampling 1x ############")
rdm1xResults=extractResultsRdm(matchedID_rdm1x,ds,fmod=NULL,id=id,exposed=exposed,outcome=outcome[i],verbose=FALSE)

for(j in 1:Ndraw){
  ############### random bagging (repeatedly downsample major class Ndraw times WITH replacement) ############
  matchedID_bag=cbind(ds[ds[,exposed]==1,id],sample(ds[ds[,exposed]==0,id],nminor,replace=TRUE))
  rdmBag[[j]]=extractResultsRdm(matchedID_bag,ds,fmod=NULL,id=id,exposed=exposed,outcome=outcome[i],verbose=FALSE)
  
  ############### random jackknife (repeatedly downsample major class Ndraw times WITHOUT replacement) ############
  matchedID_JK=cbind(ds[ds[,exposed]==1,id],sample(ds[ds[,exposed]==0,id],nminor))
  rdmJK[[j]]=extractResultsRdm(matchedID_JK,ds,fmod=NULL,id=id,exposed=exposed,outcome=outcome[i],verbose=FALSE)
}

temp=matrix(unlist(lapply(rdmBag,`[`,"results")),nrow=Ndraw,ncol=length(rdmBag[[1]][[1]]),byrow=T)
colnames(temp)=names(rdmBag[[1]][[1]])
rdmBagResults=list(results=colMeans(temp,na.rm=T),matchedID=NULL)

temp=matrix(unlist(lapply(rdmJK,`[`,"results")),nrow=Ndraw,ncol=length(rdmJK[[1]][[1]]),byrow=T)
colnames(temp)=names(rdmJK[[1]][[1]])
rdmJKResults=list(results=colMeans(temp,na.rm=T),matchedID=NULL)

#consolidate matchedIDs
matchedID[[i]]=list(expertResults[[2]],hdResults[[2]],lassoResults[[2]],rfResults[[2]],
                    euclideanResults[[2]],
                    jaccardResults[[2]],diceResults[[2]],cosineResults[[2]],
                    pearsonResults[[2]],spearmanResults[[2]],
                    lassoMVResults[[2]],
                    rdm1xResults[[2]],rdmBagResults[[2]],rdmJKResults[[2]])
names(matchedID[[i]])=c("expertPS","hdPS","lassoPS","rfPS",
                        "euclidean","jaccard","dice","cosine","pearson","spearman",
                        "lassoMV","rdmSamp1x","rdmBag","rdmJK")

############consolidate results
resultsmat=as.data.frame(rbind(expertResults[[1]],hdResults[[1]],lassoResults[[1]],rfResults[[1]],
                               euclideanResults[[1]],
                               jaccardResults[[1]],diceResults[[1]],cosineResults[[1]],
                               pearsonResults[[1]],spearmanResults[[1]],
                               lassoMVResults[[1]],
                               rdm1xResults[[1]],rdmBagResults[[1]],rdmJKResults[[1]]))
rownames(resultsmat)=names(matchedID[[i]])
resultsArray[,,i]=as.matrix(resultsmat[,desiredOrder])

}   #end of Noutcomes


dimnames(resultsArray)=list(rownames(resultsmat),desiredOrder,outcome)
names(matchedID)=outcome[1:9]

save(rfPsMod,rfImpt,lassoPsmod,lassoMVvar,lassomodbootstrap,lassoMvmod,
     ps_exp,psl,ps_rf,wantedvar,hdpsobj,
     file="PAD10sparse_models.RData")

save(resultsArray,matchedID, file="PAD10sparse_results.RData")

#remember to terminate the cores when done
stopCluster(cl)

#load(file="PAD10sparse_methods.RData")

# resultsArray
# 
# round(apply(resultsArray,c(1,2),mean,na.rm=T),3)
# apply(time[,"elapsed",],1,summary)
# apply(rfImpt,1,summary)
# 
# rowSums(resultsArray[,"inCI_matched",],na.rm=T)  #excludes lassoMV
# rowSums(resultsArray[,"inCI_adj",],na.rm=T)  #for lassoMV
# round(rowMeans(abs(resultsArray["expertPS",c("n1","KS","KL","Cstat","se_adj","bias_adj","se_matched","bias_matched"),]),na.rm=T),3)
# round(rowMeans(abs(resultsArray["hdPS",c("n1","KS","KL","Cstat","se_adj","bias_adj","se_matched","bias_matched"),]),na.rm=T),3)
# round(rowMeans(abs(resultsArray["lassoPS",c("n1","KS","KL","Cstat","se_adj","bias_adj","se_matched","bias_matched"),]),na.rm=T),3)
# round(rowMeans(abs(resultsArray["rfPS",c("n1","KS","KL","Cstat","se_adj","bias_adj","se_matched","bias_matched"),]),na.rm=T),3)

#rowMeans(rbind(trueORvec,resultsArray[,"ORadj",]),na.rm=T)
