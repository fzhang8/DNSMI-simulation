To run the simulation:

1, create a folder in your computer, say, Folder1

2, Download all the files including the .Rdata into Folder1

3, In Rstudio, run the code "script_singlepath_DNSMI_MC_Accuracy.R" by typing: source("./script_singlepath_DNSMI_MC_Accuracy.R")

4, After the running is completed, the outcome is stored in the matrix "montecarloout" with each column representing the result from a Monte Carlo run and each row being a measurement as indicated by the row names.



Meaning of some of the parameters that can be set by user:

1, currentfolder

This is the scenario setting parameter with the possible setting values being one of "nw1, nw10, nw15, nw1010, nw1015, nw10s15, nw10s25".

nw1: Low signal to nonsignal ratio

nw10: Middle signal to nonsignal ratio / Low sample size / High sparsity

nw15: High signal to nonsignal ratio

nw1010: Middle sample size

nw1015: High sample size

nw10s15: Middle sparsity

nw10s25: Low sparsity


2, MC

This is the Monte Carlo number setting. 


3, gridc1 

This is the marginal number for search gird on c1 dimension.


4, gridc2 

This is the marginal number of search gird on c2 dimension.

The number of search grid points will be gridc1*gridc2.


5, repeatsN 

This is the number of subsamples for DNSMI.
