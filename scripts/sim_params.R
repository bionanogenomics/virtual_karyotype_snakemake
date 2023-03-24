#coverage
Coverages <- 45 
min.len <- 150e3
max.len <- 10e6
avg.len <- 275e3

#FP and FN of labels
FNs <- c('1'= 0.095)
FPs <- c('1'=  1.2)
#For simulating from real data where the molecule already have some base FP/FN rates
FNs.base <- c('1'=0.0)
FPs.base <- c('1'=0.0)

#sizing error parameters
#global sizing error
meanStretch <- 1.035
sdStretch <- 0.0056


#per-interval sizing error
#sd <- -0.10
sf <- 0.112
sr <- 0.014
sd <- -0.037#-1*(4*sf^2*sr^2)^0.25  #we automatically calculate a valid sd base on sf and sr


#Knots, Fold and Chimera
fold.prob <- 0.0
knot.prob <- 0.2     #noted knot probability are now length dependent, it is specified as number of knots/100kb of sequence
Breakages.prob <- 0.0 #breakage probability are deprecated
chimera.prob <- 0.2

#label resolution parameters
res <- 2.4
resSD <- 0.28

ref.aligner.path <- '/home/jestabrook/repositories/Solve3.7_20221104_163244_27/RefAligner/12432.12642rel/RefAligner'
