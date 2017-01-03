
# Parse arguments
args = commandArgs(TRUE)
action=args[1] # should be either "construct" or "predict"
INT = args[2]
methods= args[3]; methods = unlist(strsplit(methods, ",")) # comma-separated string of used methods
lmPriority= args[4]
lmThreshold= as.numeric(args[5])
dir.create(paste0(INT,"/LM"), recursive = T, showWarnings = F)

# Install required packages
list.of.packages = c("ROCR","boot","caret","MuMIn")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(ROCR); library(boot); library(caret); library(MuMIn)

# Load all genes
genes = read.table(paste0(INT,"/SimulatedGenes/GenesinBams.txt")); colnames(genes)[1]="gene"
hits = read.table(paste0(INT,"/TrueEventlist.txt")); colnames(hits)=c("gene", "asType")
genes[,"hits"]=FALSE; genes[genes[,1]%in%hits[,1],"hits"]=TRUE
genes[,"found"]=FALSE
for (mname in methods){
  genes[,mname] = "none"
  for (type in c("A3","A5","RI","SE")){
    met = read.table(paste0(INT,"/GSEA/",mname,"/LeadingEdge_genes/LE",type,"genes.txt"))
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
options(na.action = "na.fail")
iFolds = caret::createFolds(1:nrow(genes), k = k, list = T)
rocVector = rep(NA,k)
cvct = data.frame(accuracy=rep(NA,k),TPR=NA,TNR=NA)

# Prepare drawing roc curve
cl = rainbow(10)
svg(paste0(INT,"/LM/LMRoc.svg"))

for(i in 1:k){
  # Hold out subset
  genesTraining = genes[-unlist(iFolds[i]),]
  genesTest = genes[unlist(iFolds[i]),]
  
  # Construct and fit linear model
  mod = glm(paste0("hits~",paste(methods,collapse="+")),family=binomial(link="logit"),data=genesTraining)
  dredgeObj = MuMIn::dredge(mod, rank=AIC)[1]
  for (n in methods){
    levels(dredgeObj[,n]) = c("+","-")
    dredgeObj[which(is.na(dredgeObj[,n])),n] = "-"
  }
  bestModVars = names(dredgeObj[,c(1,which(dredgeObj=="+"))])[-1]
  mod = glm(paste0("hits~",paste(bestModVars,collapse="+")),family=binomial(link="logit"),data=genes)
  
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
  if(lmPriority == "TNR"){
    rocThreshold = min(rocData[rocData[,"TNR"]>lmThreshold,3])
  }
  
  # Fill accuracy table
  cvct[i,c("TPR","TNR")] = rocData[rocData[,"alpha"]==rocThreshold,c("TPR","TNR")]
  compTest[compTest$prediction>rocThreshold,"predBinom"] = T
  ct = table(compTest[,"label"],compTest[,"predBinom"])
  cvct[i,"accuracy"] = sum(diag(ct))/sum(ct)
  
  # Draw ROC-curve
  rocCurve = performance(mod.pred, measure="tpr", x.measure="fpr")

  if(i==1){
    plot(rocCurve, col=cl[i], main="ASsembly linear model performance (10-fold)", xlim=c(0,1), ylim=c(0,1), xaxs="i",yaxs="i")
    dredgedMods = dredgeObj
  }else{
    plot(rocCurve, col=cl[i], add=T)
    dredgedMods[i,] = dredgeObj
  }
  rocVector[i] = rocThreshold
  
}
abline(a=0,b=1,lty=2)

# Add point for plain union of methods
genes[,"union"] = F
for (i in 1:nrow(genes))
  genes[i,"union"] = all(genes[i,methods]!="none")
unionTpr = sum(genes$union & genes$hits) / sum (genes$hits)
unionFpr = sum(genes$union & !genes$hits) / sum (!genes$hits)
points(unionFpr, unionTpr, pch=19)
dev.off()

# Build final model on all data
for (n in methods){
  levels(dredgedMods[,n]) = c("+","-")
  dredgedMods[which(is.na(dredgedMods[,n])),n] = "-"
  # if (sumdredgedMods[,n])
}
bestMod = data.frame(table(dredgedMods[,2:5]))
bestMod = bestMod[bestMod$Freq == max(bestMod$Freq), ]
bestModVars = names(bestMod[which(bestMod=="+")])
mod = glm(paste0("hits~",paste(bestModVars,collapse="+")),family=binomial(link="logit"),data=genes)

# Calculate roc threshold as average of values obtained through k-fold cv
bestModIndex = apply(dredgedMods, 1, function(r) all(r[2:(length(methods)+1)] == bestMod[1:length(methods)]))
mod$rocThreshold = mean(rocVector[bestModIndex])
save(mod,file=paste0(INT,"/LM/ASsemblyLM.mod"))
write.csv(cvct, file =paste0(INT,"/LM/accuracyTable.csv") )

} else if(action=="predict"){
######PREDICTION USING EXISTING MODEL######
  load(paste0(INT,"/LM/ASsemblyLM.mod"))
  predTest = predict(mod,genes)
  predTest=boot::inv.logit(predTest)
  compTest=data.frame(prediction=predTest,predBinom=F)
  compTest$predBinom[compTest$prediction>mod$rocThreshold]=T
  write.table(compTest,paste0(INT,"/LM/ASsemblyLMresults.txt"))
}