
# Parse arguments
args = commandArgs(TRUE)
action=args[1] # should be either "construct" or "predict"
INT = args[2]
methods= args[3]; methods = unlist(strsplit(methods, ",")) # comma-separated string of used methods
lmPriority= args[4]
lmThreshold= as.numeric(args[5])
dir.create(paste0(INT,"/LM"), recursive = T, showWarnings = F)

# Install required packages
list.of.packages = c("ROCR","boot","caret")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(ROCR); library(boot); library(caret)

# Load all genes
genes = read.table(paste0(INT,"/SimulatedGenes/GenesinBams.txt")); colnames(genes)[1]="gene"
hits = read.table(paste0(INT,"/TrueEventlist.txt")); colnames(hits)=c("gene", "asType")
genes[,"hits"]=FALSE; genes[genes[,1]%in%hits[,1],"hits"]=TRUE
genes[,"found"]=FALSE
for (mname in methods){
  genes[,mname] = "none"
  for (type in c("A3","A5","RI","SE")){
    met = read.table(paste0(INT,"/GSEA/",mname,"/LE",type,"genes.txt"))
    genes[genes[,1]%in%met[,1],mname] = type
    genes[genes[,1]%in%met[,1],"found"] = TRUE
  }
  genes[,mname] = as.factor(genes[,mname])
}
genes=genes[genes$found==TRUE,]
row.names(genes)=genes[,1];genes = genes[,-c(1,3)]

######MODEL BUILDING######
if (action=="construct"){
k=10
iFolds = caret::createFolds(1:nrow(genes), k = k, list = T)
rocVector = rep(NA,k)

# Prepare drawing roc curve
cl = rainbow(10)
# svg(paste0(INT,"/LM/LMRoc.svg"))

for(i in 1:k){
  # Hold out subset
  genesTraining = genes[-unlist(iFolds[i]),]
  genesTest = genes[unlist(iFolds[i]),]
  
  # Fit linear model
  mod = glm(paste0("hits~",paste(methods,collapse="+")),family=binomial(link="logit"),data=genesTraining)
  
  # Perform ROC analysis
  predTest = predict(mod,genesTest)
  predTest = boot::inv.logit(predTest)
  compTest = data.frame(prediction=predTest,label=genesTest[,"hits"],predBinom=F)
  mod.pred = prediction(compTest[,1],compTest[,2])
  perf = performance(mod.pred, measure="tpr", x.measure="tnr")
  rocData = data.frame(TNR=unlist(perf@x.values), TPR=unlist((perf@y.values)),alpha=unlist(perf@alpha.values))
  
  # Select threshold, store in glm object
  if(lmPriority == "TPR")
    rocThreshold = max(rocData[rocData[,"TPR"]>lmThreshold,3])
  if(lmPriority == "TNR")
    rocThreshold = min(rocData[rocData[,"TNR"]>lmThreshold,3])
  
  # # Draw ROC-curve
  # rocCurve = performance(mod.pred, measure="tpr", x.measure="fpr")
  # # auc = performance(mod.pred, measure = "auc")
  # if(i==1)
  #   plot(rocCurve, col=cl[i], main="ASsembly linear model performance (10-fold)", xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
  # else
  #   plot(rocCurve, col=cl[i], add=T)
  # text(1,0.1,pos=2,paste("AUC = ",round(auc@y.values[[1]],3), sep=""))
  rocVector[i] = rocThreshold
}
abline(a=0,b=1,lty=2)

# Add point for plain union of methods
genes[,"union"] = F
genes[genes$MATS!="none" | genes$DEXSeq!="none" | genes$ASpli!="none" | genes$Cufflinks!="none","union"] = T
genes[,"union2"] = F
for (i in 1:nrow(genes))
  genes[i,"union2"] = any(genes[i,methods]!="none")
unionTpr = sum(genes$union2 & genes$hits) / sum (genes$hits)
unionFpr = sum(genes$union2 & !genes$hits) / sum (!genes$hits)
points(unionFpr, unionTpr, pch=19)

# Build final model on all data
mod = glm(hits~.,family=binomial(link="logit"),data=genes)

# Calculate roc threshold as average of values obtained through k-fold cv
mod$rocThreshold = mean(rocVector)

save(mod,file=paste0(INT,"/LM/ASsemblyLM.mod"))


} else if(action=="predict"){
######PREDICTION USING EXISTING MODEL######
  load(paste0(INT,"/LM/ASsemblyLM.mod"))
  predTest = predict(mod,genes)
  predTest=boot::inv.logit(predTest)
  compTest=data.frame(prediction=predTest,predBinom=F)
  compTest$predBinom[compTest$prediction>mod$rocThreshold]=T
  write.table(compTest,paste0(INT,"/LM/ASsemblyLMresults.txt"))
}