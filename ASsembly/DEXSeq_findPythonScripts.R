# Install packages if necessary (suppose this doesn't apply once we go run it for real)

list.of.packages = c("rPython")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!any(installed.packages()[,"Package"]=="BiocInstaller") | !any(installed.packages()[,"Package"]=="BiocParallel"))
  source("https://bioconductor.org/biocLite.R")

list.of.biocPackages = c("DEXSeq")
new.packages = list.of.packages[!(list.of.biocPackages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocInstaller::biocLite(new.packages)


pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
if(any(list.files(pythonScriptsDir)== "dexseq_count.py")){pythonScriptsDir
  }else{stop("Cannot find DEXSeq python script", call.=FALSE)}
