# Check if the package is already installed
if (!requireNamespace("lfmm", quietly = TRUE)) {
  # If the package is not installed, install it
  library(devtools)
  devtools::install_github("bcm-uga/lfmm")
}

library(lfmm)

geno_file = snakemake@input[['geno_file']]
env_file = snakemake@input[['env_file']]
num_components_file = snakemake@input[['num_components']]

pvalues_file = snakemake@output[['p_values_lfmm']]
qqplot_file = snakemake@output[['qq_plot']]

print(num_components_file)
# Read the num_components
num_components <- scan(num_components_file, what = integer(), n = 1)
#num_components <- readLines(num_components_file, warn = FALSE) 
#num_components <- as.integer(num_components[1])

# Read the geno and env tables 
geno <- read.csv(geno_file, sep = ',', header = TRUE)
geno <- geno[, !colnames(geno) %in% "chrom_pos"]
geno <- sapply(geno, as.numeric)
Y <- as.matrix(geno)
X <- as.matrix(read.csv(env_file, sep = '', header = FALSE))

# Check the number of columns in Y
if (ncol(Y) <= 3) {
  # Create the p-values file with 0 inside
  write.csv(0, file = pvalues_file, row.names = FALSE)
  
  # Create the empty PNG file
  file.create(qqplot_file)
} else {
  # Proceed with LFMM analysis
  
  ## Transpose Y because of the format needed for lfmm 
  Y_t <- t(Y)
  
  ## Fit an LFMM, i.e., compute B, U, V estimates
  ## using 5 because that is what I estimated based on the screeplot 
  mod.lfmm <- lfmm_ridge(Y = Y_t, 
                         X = X, 
                         K = num_components)
  
  ## Perform association testing using the fitted model:
  pv <- lfmm_test(Y = Y_t, 
                  X = X, 
                  lfmm = mod.lfmm, 
                  calibrate = "gif")
  
  pvalues <- pv$calibrated.pvalue
  
  ## Write p-values
  write.csv(pvalues, file = pvalues_file, row.names = FALSE)
  
  # Set the file path and name for the PNG file
  png_file <- qqplot_file
  
  # Open the PNG file device
  png(filename = png_file)
  
  # Generate the plot
  qqplot(rexp(length(pvalues), rate = log(10)),
         -log10(pvalues), xlab = "Expected quantile",
         pch = 19, cex = 0.4)
  abline(0, 1)
  
  # Close the device and save the file
  dev.off()
}