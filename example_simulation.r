#### Description ####
# A script like this used to generate the simulations for the paper: X Chromosome Association Testing in Genome Wide Association Studies. Peter Hickey (peter.hickey@gmail.com) 23/12/2010
# This script is an example of how to use the simulation code.

#### Load necessary functions ####
source("/Users/hickey/Desktop/chrX_GWAS_simulation_R_code_example/functions.R") # This is the path to the functions.R file in XGWAS. The user will need to specify this path.
testNames <- c("Z_A.sq", "Z_C.sq", "Z_mfA.sq", "Z_mfG.sq", "S_1", "S_2", "plink.trend.f", "plink.geno.f", "chi4df") # This shouldn't be altered by the user. Names of tests (in order).
results <- c("MAF", "Model", "t", "case", "control", testNames) # This shouldn't be altered by the user. Header line for results table. 

#### Specify simulation parameters ####
reps <- 10000 # reps is the number of replicates per parameter combination
alpha <- 0.0001 # alpha is the nominal significance level
MAF <- 0.2 # MAF is the minor allele frequency, which we assume to be the risk allele frequency. Can be a vector.
CASE <- 0.5 # The proportion of cases that are female. Can be a vector.
CONTROL <- 0.5 # The proportion of control that are female. Can be a vector (if a vector, should be the same length as the CASE variable).
model <- "ADD" # The genetic model: "ADD" (additive constraint), "DOM" (dominant constraint), "REC" (recessive constraint) or "Null" (null model; used for computing size of tests)
n_cases <- 6000 # The number of cases
n_controls <- 6000 # The number of contols
k <- 0.02 # The population prevalence of the phenotype
t <- 1.25 # The risk in the homozygous BB genotype (see paper for details)

#### Loop over all combinations of parameters, simulate data and compute value of each test statistic on each simulated dataset. WARNING: This isn't the most efficient way of doing things in R but it is perhaps the simplest way to convey what operations are being performed. ####
for(maf in 1:length(MAF)){
	for(cohort in 1:length(CASE)){
		test_results <- XSim(n_reps = reps, N_cases = n_cases, N_controls = n_controls, case.sex.ratio = CASE[cohort], control.sex.ratio = CONTROL[cohort], MAF = MAF[maf], MODEL = model, t = t, K = k, results = NULL, alpha = alpha) # test_results contains the value of each test statistic under each simulation scenario
		params <- c(MAF[maf], model, t, CASE[cohort], CONTROL[cohort]) # params is the parameters used by each simulation
		paramsPlusResults <- append(params, makeTable(test_results, reps, alpha)) # paramsPlusResults stores the parameters plus the test_results for each combination of parameters
		results <- rbind(results, paramsPlusResults) # results simply tallies up all the paramsPlusResults instances.
	}
}

#### Print output files ####
write.table(results, quote = F, sep = "\t", row.names = F, col.names = F, file = "/Users/hickey/Desktop/chrX_GWAS_simulation_R_code_example/example_results.txt") # The user should specify this path. The resulting file is a tab-delimited version of 'results'
parameters <- paste("reps = ", reps, ", alpha = ", alpha, ", t = ", t, sep = "")
write(parameters, file = "/Users/hickey/Desktop/chrX_GWAS_simulation_R_code_example/example_parameters.txt") # The user should specifiy this path. The resulting file is a comma-separated file of the parameters used in each simulation.

#### FINISHED ####
