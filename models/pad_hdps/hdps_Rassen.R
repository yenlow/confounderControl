# Adapted from HDPS_PROTOTYPE.R by Jeremy A. Rassen <jrassen@post.harvard.edu>
# 
# 24-Mar-14 Yen Low
##############################################################################
#

library(rJava) #hdps
#library(rmysql)
library(flexmix) #KLdiv
library(ROCR) #Cstat
library(Matching)
library(survival) #clogit
library(Epi) #clogistic

#############SET PARAMETERS######################
setwd("/mnt/hgfs/projects/psm/models/pad_hdps/")
getwd()

source("../../hdpsFunctions.R")

baseDir = "/mnt/hgfs/projects/psm/models/pad_hdps/";
tempDir = paste(baseDir,"tmp",sep="") #output directory
patientfile=paste(baseDir,"patients.txt",sep="")
dir.create(tempDir)

# initialize the Java subsystem.  include a jdbc object
.jinit(classpath="/home/yenlow/software/mysql-connector-java-5.1.26/mysql-connector-java-5.1.26-bin.jar:/home/yenlow/software/pharmacoepi.jar", parameters="-Xmx1g -Xms1g", force.init = TRUE);
.jclassPath()

#########INPUT#######################
# Read in the original patients file 
load("../../data/pad/PAD_PLOSOnePatientData.RData",verbose=T)
#load(paste(baseDir,"PAD_PLOSOnePSMGLM1.RData",sep="")) #matching variables
dataset=variables.data.final
rm(variables.data.final) #to save memory

id="pid_org"
exposed="cilostazol" #name of exposure variable (must be string)
outcome=c(	"mace.after","cardiac.arrest.after","death","ec.after","mi.after","stroke.after","sudden.cardiac.death.after",             
			"male.after","amputation.after","angioplasty.after","bypass.after","revasc.after",
			"arrythmias.after","af.after","bradycardia.after","tachycardia.after","vf.after","vt.after",
			"dizziness.after","palpitations.after")
#demographic and expert selected variables must be matched
demvariables=c("min_pad_age_scaled", "gender", "race")
#demvariables=c("gender.male","gender.female","gender.unknown" ,
#				"race.white","race.asian","race.black",
#				"race.nativeamerican","race.other" ,"race.unknown",
#				"age.pediatric","age.youngadults","age.middleaged","age.geriatric")
expertvariables=c(	"chf.before", "dyslipidemias.before", "renal.failure.before",
					"aspirin.before", "clopidogrel.before", "warfarin.before","antiarrhythmics.before",
					"mace.before", "male.before",
					"hypertension.before","statins.before", "ace.inhibitors.before","arrythmias.before")
#expertvariables unused in Anna's code:
#"hypertension.before","statins.before", "ace.inhibitors.before","arrythmias.before","revasc.before","bypass.before","angioplasty.before"
#other variables will NOT be matched
othervariables=c("min_pad_timeoffset", "followup", "min_cil_age", "min_cil_timeoffset", "timepoint","note.count","psm.score")
#empirical variables may be matched depending on HD-PS
empvariables=colnames(dataset)[grep("before",colnames(dataset))]
#		"hypertension.before","statins.before", "ace.inhibitors.before","arrythmias.before",
#colnames(dataset)[grepl("after$",colnames(dataset))]
#setdiff(colnames(dataset),c(id,exposed,outcome,expertvariables,othervariables))
#wantedvariables=c(exposed,outcome,expertvariables)
#wantedvariables %in% colnames(dataset) #check variables are in main data table


###### expertPSM ######
#scale continuous variables
dataset[,"min_pad_age_scaled"]=scale(dataset[,"min_pad_age"])
par(mfrow=c(1,2))
hist(dataset[,"min_pad_age"])
hist(dataset[,"min_pad_age_scaled"])

# Create formula: exposure ~ var1 + var 2 + ...
fwithdemvariables=paste(demvariables,sep="",collapse=" + ")
fwithexpertvariables=paste(expertvariables,sep="",collapse=" + ")
fPS = as.formula(paste(exposed," ~ ",paste(fwithdemvariables,fwithexpertvariables,sep=" + "),sep=""))
fPS

xmat=as.matrix(dataset[,c(demvariables,expertvariables)])
xmat=cbind(	model.matrix(~gender+race-1,dataset)[,c("genderMale","raceAsian","raceBlack","raceNative American","raceWhite")],
			xmat[,-grep("gender$|race$",colnames(xmat))])
mode(xmat)="numeric"
head(xmat,3)

#calculate PS
expertPsMod=glm(fPS, data=dataset[,c(id,exposed,demvariables,expertvariables)], family="binomial")
ps=expertPsMod$fitted.values
names(ps)=dataset[,id]

#Define optimal caliper
#using Austin 2011 recommendation of 0.2*sd(logit(PS)) - assumes normality after transformation
#Rosenbaum recommends 0.25sd(PS) in his book
#the two converges if ps is small
optcaliper=0.2*sd(logit(ps)) #logit=predict(expertPsMod)

# Downsample by matching
#using hdPS (calls Java org.drugepi.match.Match)
expertPsM=hdpsMatch(dataset[,c(id,exposed)],ps,mode="GREEDY_CALIPER",k=1,caliper=optcaliper,output_matchfile="matches_expertPS.txt")

#diagnostics for matching
expertPSresults=psDiag(expertPsM,glm=expertPsMod,"matchedplots_expertPS.pdf")

#using Matching package (uses number of sd of PS, not logit(PS)
expertPsM1=Match(Tr=expertPsMod$y, X=expertPsMod$fitted, replace=F, M=1, caliper=optcaliper/sd(ps), ties=F)
expertPsM1$ecaliper #check actual caliper used is same as optcaliper
expertPSresults_Matching=psDiag(expertPsM1,glm=expertPsMod,"matchedplots_expertPS_Matching.pdf")
#mb1=MatchBalance(fPS, data=dataset[,c(id,exposed,demvariables,expertvariables)],
#		match.out=expertPsM1,nboots=0)

save(expertPSresults,expertPsMod,expertPsM,file="psm_pad.RData")

######################
# get OR #
i=1
logitdata=merge(dataset[,c(id,exposed,"psm.score",outcome[i])],expertPsM$matches[,c("pat_id","ps","set_num")],by.x=id,by.y="pat_id",all.y=T,sort=F)
logitdata=logitdata[order(logitdata$set_num),]
logitdata$set_num=as.factor(logitdata$set_num)
modExpertPSadj=clogistic(paste(outcome[i],"~",exposed,"+ ps"),strata=set_num, data=logitdata)

expertOR=getOR(modExpertPSadj)
expertSmd=smd(logitdata,exposed=exposed,outcome=outcome[i])




####### HDPS ######################

### INVOKE pharmacoepi.jar ###
#read data file as single string (to use addPatientsFromBuffer)
for(i in 1:length(outcome)){
	patientheader=paste(id,exposed,outcome[i],sep="\t")
	patientstring=paste(paste(dataset[,id],dataset[,exposed],dataset[,outcome[i]],sep="\t"),collapse="\n")
#for reading in whole file as single string
#patientstring=readChar(patientfile, file.info(patientfile)$size)
#patientstring=gsub("\"","",patientstring)
#patients=dataset[,c("pid_org",exposed,outcome[1])] #get patient id, exposure, outcome
#write.table(patients,file=patientfile,sep="\t",col.names=T,row.names=F,na="") #output data to patient.txt;
}


datainstring=paste(patientheader, patientstring, sep="\n")
variables=unique(c(expertvariables,empvariables))
dimdata=dataset[,c(id,variables)]

hdpsobj=hdps(datainstring,dimdata,outDir=tempDir,Nmostfreq=30,k=20,outfile="output_cohort.txt",FullOutput=T,verbose=T,ZeroCellCorrection=F)
selectedvariables=colnames(hdpsobj$hdpsdata)[-1]


### Calculate propensity scores using R ###
# Estimate the PS
dataPS=merge(dataset[,c(id,exposed,demvariables)],hdpsobj$hdpsdata,by.x=id,by.y="patient_id",all.x=T,sort=F)
hdPsMod=glm(paste(exposed,"~ . -",id), data=dataPS, family="binomial")
names(hdPsMod$fitted.values)=as.character(dataPS$pid_org)

#create cohort data matrix: should contain 3 columns (id, exposed, ps)
ps=hdPsMod$fitted.values
names(ps)=dataPS[,id]
optcaliper_logit=0.2*sd(logit(ps))
optcaliper=exp(optcaliper_logit)/(1+exp(optcaliper_logit))

### Downsample by matching (calls Java org.drugepi.match.Match)####
hdpsM=hdpsMatch(dataPS[,c(id,exposed)],ps,mode="GREEDY_CALIPER",k=1,caliper=optcaliper,output_matchfile="matches_hdPS.txt")

#diagnostics for matching
hdPSresults=psDiag(hdpsM,glm=hdPsMod,outfile="matchedplots_hdPS.pdf")


### get OR ###
logitdata=merge(dataset[,c(id,exposed,othervariables,outcome[i])],hdpsM$matches[,c("pat_id","ps")],by.x=id,by.y="pat_id",all.y=T,sort=F)
modHdPSadj=glm(paste(outcome[i],"~ cilostazol + ps"), data=logitdata, family="binomial")
#modHdPSadj=glm(paste(outcome[i],"~ cilostazol + psm.score"), data=logitdata, family="binomial")

hdpsOR=getOR(modHdPSadj)
hdpsSmd=smd(logitdata,exposed=exposed,outcome=outcome[i])

plot(hdPsMod$fitted.values,expertPsMod$fitted.values)

save(	hdpsobj,hdpsM,hdPsMod,
		expertPsM,expertPsMod,expertPSresults,
		file="hdps.RData")
