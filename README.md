# Replication files for "Mixing in Early Childhood:  Evidence from SyrianRefugees in Turkey"

## Vincent Boucher, Semih Tumen, Michael Vlassopoulos, Jackline Wahba, and Yves Zenou

+ The file *simgmm.R* is a R script that has to be run FIRST. It produces the estimates for Table 9 (Panel A) and Table A6. Results are saved and used by the other scripts.
+ The file *postestimation.R* is a R script that has to be run AFTER *simgmm.R*. It produces the estimates for Table 9 (Panel B) and Table A7. The code also produces Figures 3--18.
+ The file *model_fit.R* is a R script that has to be run AFTER *simgmm.R*. It produces the estimates for Table 10 and Table 11.
+ The file *netgraph.R* is a R script that has to be run AFTER *simgmm.R*. It produces Figures 1--2. 
+ The file *fcts_gmm.R* contains all the functions used by the other scripts.
