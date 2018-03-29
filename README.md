To run the simulation:

1, create a folder in your computer, say, Folder1

2, Download all the files including the .Rdata into Folder1

3, In Rstudio, run the code "script_singlepath_DNSMI_MC_Accuracy.R" by typing: source("./script_singlepath_DNSMI_MC_Accuracy.R")

4, After the running is completed, the outcome is stored in the matrix "montecarloout" with each column representing the result from a Monte Carlo run and each row being a measurement as indicated by the row names.



Meaning of some of the parameters that can be set by user:

1, currentfolder

This is the scenario setting parameter with the possible setting values being one of "nw1, nw10, nw15, nw1010, nw1015, nw10s15, nw10s25"

nw1: Low signal to nonsignal ratio
nw10:
