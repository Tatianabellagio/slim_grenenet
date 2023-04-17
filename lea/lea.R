setwd("/Users/tbellagio/Documents/grenephase1/scratch_tati/simulations/lea/slim_grenenet/")
library(LEA)
## well their functions for data format conversion are a piece of shit, so i had to do it myself in python 
# load simulated data


# run of pca
# Available options, K (the number of PCs),
#                    center and scale.
# Create files: genotypes.eigenvalues - eigenvalues,
#               genotypes.eigenvectors - eigenvectors, #               genotypes.sdev - standard deviations, #               genotypes.projections - projections,
# Create a pcaProject object: pc.
pc = pca("allele_freq.lfmm", scale = TRUE)
tw = tracy.widom(pc)

# plot the percentage of variance explained by each component 
plot(tw$percentage, pch = 19, col = "darkblue", cex = .8)


### ECOLOGICAL ASSOCIATIONS 
## lfmm() is Bayesian method that uses a Monte-Carlo Markov Chain algoritm, 
## and lfmm2() is  a  frequentist  approach  that  uses  least-squares  estimates
# (e.g. more than 1,000-10,000 genetic loci), the best is to use lfmm2().


# main options:
# K = the number of latent factors
## BECAUSE WE SAW FROM THE PCA that there seem to be 6 clusters..

# Runs with K = 6 using 5 repetitions.
project = NULL
project = lfmm("genotypes.lfmm",
               "gradients.env",
               K = 6,
               repetitions = 5,
               project = "new")
## lfmm uses a very naive imputation method which has low power
when genotypes are missing:  See impute() for a better imputation method.
## Note that lfmm has an improved estimation algorithm implemented in lfmm2, which should be the prefered option.





data <- read.csv("allele_freq.csv", header=FALSE)
env = read.csv("env.env", header=FALSE)


#data("offset_example")
# 200 diploid individuals genotyped at 510 SNP 
#Y <- offset_example$geno
# 4 environmental variables
#X <- offset_example$env

mod.lfmm2 <- lfmm2(input = data, env = env, K = 5)

#Simulate non-null effect sizes for 10 target loci #individuals
n = 100
#loci
L = 1000
# Environmental variable
X = as.matrix(rnorm(n))
# effect sizes
B = rep(0, L)
target = sample(1:L, 10)

# GEA significance test
# showing the K = 2 estimated factors
plot(mod.lfmm2@U, col = "grey", pch = 19, xlab = "Factor 1",
     ylab = "Factor 2")

B[target] = runif(10, -10, 10)

Create 3 hidden factors and their loadings
U = t(tcrossprod(as.matrix(c(-1,0.5,1.5)), X)) + matrix(rnorm(3*n), ncol = 3)
V <- matrix(rnorm(3*L), ncol = 3)


pv <- lfmm2.test(object = mod.lfmm2,
                 input = data,
                 env = env,
                 full = TRUE)
plot(-log10(pv$pvalues), col = 'grey', cex = .5, pch = 19)
abline(h = -log10(0.1/16), lty = 2, col = "orange")

-log10(pv$pvalues)
#sorted_pvalues = sort(pv$pvalues, decreasing = TRUE)
write.csv(pv$pvalues, "p_values_snps")
write.csv(-log10(pv$pvalues), "log10p_values_snps")


# Simulate a matrix containing haploid genotypes Y <-  tcrossprod(as.matrix(X), B) +
tcrossprod(U, V) +
  matrix(rnorm(n*L, sd = .5), nrow = n)
Y <- matrix(as.numeric(Y > 0), ncol = L)







