# Replication files for "Mixing in Early Childhood:  Evidence from SyrianRefugees in Turkey"

## Vincent Boucher (Laval University), Semih Tumen (TED University, IZA, and ERF), Michael Vlassopoulos (University of Southampton and IZA), Jackline Wahba (University of Southampton and IZA) and Yves Zenou (Monash University and CEPR)

+ The file *fcts_gmm.R* includes all functions for the other codes.
+ The file *simgmm.R* has to be run BEFORE *postestimation.R*. It produces the estimates for Table 9 (Socialization) and Table A6. Results are saved and used by *postestimation.R*
+ The file *postestimation.R* has to be run AFTER *simgmm.R*. It produces the estimates for Table 9 (Network) and Table A7. The code also produces Figures 1--8.
+ The file *Additional_regressions.R* produces the results in Table A8.
