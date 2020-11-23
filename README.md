# Replication files for "Mixing in Early Childhood:  Evidence from SyrianRefugees in Turkey"

## Vincent Boucher, Semih Tumen, Michael Vlassopoulos, Jackline Wahba, and Yves Zenou

### Stata codes
+ file1
+ file2

### R codes
+ The file *simgmm.R* is a R script that has to be run FIRST. It produces the estimates for Table 8 (Panel A) and Table A6. Results are saved and used by the other scripts.
+ The file *postestimation.R* is a R script that has to be run AFTER *simgmm.R*. It produces the estimates for Table 8 (Panel B) and Table A7. The code also produces Figures 1--16 and A4--A9
+ The file *model_fit.R* is a R script that has to be run AFTER *simgmm.R*. It produces the estimates for Table 9 and Table 10.
+ The file *fcts_gmm.R* contains all the functions used by the other scripts.
