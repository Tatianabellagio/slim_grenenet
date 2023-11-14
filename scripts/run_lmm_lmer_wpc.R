library(dplyr)
library(data.table)
library(doParallel)
library(stats)
library(lmerTest)
library(lme4)
## 20 threads
cl <- makeCluster(20)
registerDoParallel(cl)

## load data
env_sites = snakemake@input[["env_sites"]]
p_norm = snakemake@input[["p_norm"]]
pop_structure = snakemake@input[["pop_structure"]]
lmm_results = snakemake@output[["lmm_results"]]

env_sites = read.csv(env_sites, row.names = 1, )
print(dim(env_sites))
p_norm = read.csv(p_norm, nrows=10)
print(dim(p_norm))
deltap = subset(p_norm, select = -chrom_pos)
print(dim(deltap))
pop_strc = read.csv(pop_structure,row.names = 1)
print(dim(pop_strc))

# assemble lmm objects
prep_lmm <- function(yy, env_sites, envvar, pop_strc) {
  mydata = cbind(yy, env_sites[c('sites', 'env')], pop_strc)
  return(mydata)
}

# get lmm model results
format_lmer <- function(mymodel) {
  l_ratio = drop1(mymodel,test="Chisq") #test="Chisq"
  outdt = c(fixef(model)['env'],
            summary(model)$coefficients[, "Std. Error"]['env'],
            anova(model)$"Pr(>F)",
            fixef(model)['(Intercept)'],
            summary(model)$coefficients[, "Std. Error"]['(Intercept)'],
            BIC(model),
            p_value = l_ratio$'Pr(>Chi)'[2]
  )
  return(outdt)
}

# test
print(dim(deltap))
yy = as.numeric(unlist(deltap[1,]))
print(yy)
#.GlobalEnv$myfm <- myfm # fix a global env bug
mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 
print(mydata)

# Combine all parts of the formula
formula_str <- paste("yy ~ env +", paste0("PC", 1:10, collapse = " + "), "+ (1|sites)")
# Print the formula string
print(formula_str)


model <- lmer(formula_str, data=mydata,  REML = FALSE)
l_ratio = drop1(model,test="Chisq") #test="Chisq"

outdt = c(fixef(model)['env'], '1',
            summary(model)$coefficients[, "Std. Error"]['env'],  '2',
            anova(model)['env', 'Pr(>F)'],  '3',
            fixef(model)['(Intercept)'], '4',
            summary(model)$coefficients[, "Std. Error"]['(Intercept)'], '5',
            BIC(model), '6',
            p_value = l_ratio$'Pr(>Chi)'[2], '7'
  )

print(outdt)

#output = format_lmer(model) # output model results
#print(output)

#nrow(deltap)
functions_to_export <- c("lmer", "fixef")
lmeres = foreach(ii = 1:nrow(deltap), .combine = rbind, .export = functions_to_export ) %dopar% {
  yy = as.numeric(unlist(deltap[ii,]))
  #.GlobalEnv$myfm <- myfm # fix a global env bug
  mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 
  model <- lmer(formula_str, data=mydata,  REML = FALSE)
  l_ratio = drop1(model,test="Chisq") #test="Chisq"
  outdt = c(fixef(model)['env'],
            summary(model)$coefficients[, "Std. Error"]['env'],
            anova(model)['env', 'Pr(>F)'],
            fixef(model)['(Intercept)'],
            summary(model)$coefficients[, "Std. Error"]['(Intercept)'],
            BIC(model),
            p_value = l_ratio['env', 'Pr(>F)'])
}


print(lmeres)

dimnames(lmeres)[[2]] = c('env_value', 'env_eror', 'p_value_env','intercept_value', 'intercept_eror', 'bic', 'lrt')
#colnames(final_results_lmer) <- c('env_value', 'env_eror', 'p_value_env','intercept_value', 'intercept_eror', 'bic')

print(lmeres)

write.csv(lmeres, file = lmm_results)
