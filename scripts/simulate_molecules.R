#configuration parameters for simulations



args<-commandArgs()
#file.arg.name <- "--file="
#script.name <- sub(file.arg.name, "", args[grep(file.arg.name, args)])
#script.basename <- dirname(script.name)

#Input file and output folder can optionally parsed from
#command-line argument

cmapFile <- ""
outputFolder <- ""
genome_name <- "Genome_Sim"
params_file <- NULL

for(i in 1:length(args)){
  if(args[[i]] == '-r'){
    print(paste0("qry genome: ", args[i+1]));
    cmapFile=args[[i+1]]
  }
  if(args[[i]] == '-o'){
    print(paste0("output: ", args[i+1]));
    outputFolder=args[[i+1]]
  }
  if(args[[i]] == '-n'){
    print(paste0("name: ", args[i+1]));
    genome_name=args[[i+1]]
  }
  if(args[[i]] == '-n'){
    print(paste0("name: ", args[i+1]));
    genome_name=args[[i+1]]
  }
  if(args[[i]] == '-p'){
    print(paste0("name: ", args[i+1]));
    params_file=args[[i+1]]
  }
}

#Input cmap file
begin.time <- Sys.time()
print(paste0("Begin timestamp ", begin.time))

if(!file.exists(cmapFile) || !dir.exists(outputFolder)){
  stop("please supply a valid genome cmap file and output folder")
}

mol_file <- dir(outputFolder, paste0(genome_name, ".*mols.RData"), full.names=T)

if(length(mol_file) > 0){
  stop("There is molecule file in the output directory with the same genome name, please use a different name of different output directory")
}
library(BionanoR)
library(MolSim)
library(data.table)

#coverage
#BR: change here for requested coverage
Coverages <- 20 
min.len <- 150e3
max.len <- 1000e3
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

ref.aligner.path <- '/home/users/jwang/bin/RefAligner'


#load user parameters if specified
if(!is.null(params_file)){
  source(params_file)
}


ref <- as.data.table(read_bng_maps(cmapFile))

mol_err_free <- generate_molecules(ref, minLen=min.len, maxLen = max.len, avglen = avg.len, coverage = Coverages, outputPath = outputFolder, genome.name = genome_name)

mol_file <- dir(outputFolder, paste0(genome_name, ".*mols.RData"), full.names=T)


print(paste0("Error free molecule file: ", mol_file))

add_molecular_noises(mol_file, outputFolder,
                      FPs, FNs, FPs.base, FNs.base,
                      coverage=Coverages, sd = sd, sf = sf, sr = sr, mean.stretch = meanStretch, sd.stretch=sdStretch,
                      fold.prob=fold.prob, knot.prob=knot.prob, chimera.prob=chimera.prob, breakage.prob=Breakages.prob, res=res, resSD=resSD,
                       ref.aligner.path=ref.aligner.path)

end.time <- Sys.time()
print(paste0("End time ", end.time))
print(paste0("Time used ", end.time - begin.time))

