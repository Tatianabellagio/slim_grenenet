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
p_norm = read.csv(p_norm)
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
  l_ratio = drop1(model,test="Chisq") #test="Chisq"
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
model <- lmer('yy ~ env + PC1 + PC2 + PC3 + (1|sites)', data=mydata,  REML = FALSE)

format_lmm(model, envvar) # output model results




output = format_lmm(model, envvar) # output model results
print(output)

lmeres = foreach(ii = 1:nrow(deltap), .combine = 'rbind', .errorhandling = 'remove') %dopar% {
  yy = as.numeric(unlist(deltap[ii,]))
  .GlobalEnv$myfm <- myfm # fix a global env bug
  mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 

  model <- lmer('yy ~ env + PC1 + PC2 + PC3 + (1|sites)', data=mydata,  REML = FALSE)

  format_lmm(model, envvar) # output model results
}

dimnames(lmeres)[[2]] = c('env_value', 'env_eror', 'p_value_env','intercept_value', 'intercept_eror', 'bic', 'lrt')

print(lmeres)
print(dim(deltap))
print(dim(lmeres))

print(dim(deltap))
print(dim(lmeres))

print(lmm_results)
write.csv(lmeres, file = lmm_results)
