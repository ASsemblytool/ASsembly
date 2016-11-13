# ASsembly

*(Alpha phase)*

ASsembly is a tool for the assembly of the results of several methods
for detecting differential splicing. These methods tend to find different differential
splicing events, with little overlap within their results. Some methods work well for 
finding certain types of events but not others, and not all of them are good when novel 
isoforms appear. ASsembly aims to unify these results and find high 
confidence differential splicing events among them.

ASimulator is a tool for simulating differential splicing between two given conditions. 
These simulations can be used for testing the accuracy of differential splicing methods
for finding certain alternative splicing event types (simple and complex) and novel isoforms;
and checking their robustness against diffenrent isoform ratios between the samples, 
quality of sequencing and number of replicates.

Given the simulated files, ASsembly uses a GSEA approach to find high confidence thresholds
(e.g. for specific event types or novel isoform discovery) for each method. The simulated 
gene lists are used as reference gene sets against the ranked differentially spliced genes
lists found by the methods. Methods are then run on the original data and these thresholds
are used to identify how many of the differential spliced genes found can be trusted. 

More information can be found in the tools documentation.

**Authors: Carmen Bravo, Carlos de Lannoy, Emma Houben**
**Email: assemblytoolkuleuven@gmail.com**
**Github: https://github.com/ASsemblytool/ASsembly**
