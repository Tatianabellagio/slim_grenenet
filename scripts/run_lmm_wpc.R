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
format_lmm <- function(mymodel, envvar) {
  lmesum = summary(mymodel)
  lmer2 = MuMIn::r.squaredGLMM(mymodel) # contains two rows
  outdt = c(lmer2,
            lmesum$tTable[envvar, "Value"],
            lmesum$tTable[envvar, "p-value"],
            lmesum$BIC,
            )
  return(outdt)
}

## set up the model
envvar = 'env'
myfm = as.formula(paste0('yy ~ ', paste(c(paste0('PC',1:3), envvar), collapse = ' + ')))
print(myfm)

print(dim(deltap))
yy = as.numeric(unlist(deltap[1,]))
print(yy)
#.GlobalEnv$myfm <- myfm # fix a global env bug
mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 
print(mydata)
model = nlme::lme(fixed = myfm, random = ~ 1|sites, data = mydata, method = 'ML') # no popstr PCs

myfm_reduced = as.formula(paste0('yy ~ ', paste(c(paste0('PC',1:3)), collapse = ' + ')))
reduced_model = nlme::lme(fixed = myfm_reduced, random = ~ 1|sites, data = mydata, method = 'ML') # no popstr PCs

lrt_result <- anova(model, reduced_model, test = "Chisq")




output = format_lmm(model, envvar) # output model results
print(output)

lmeres = foreach(ii = 1:nrow(deltap), .combine = 'rbind', .errorhandling = 'remove') %dopar% {
  yy = as.numeric(unlist(deltap[ii,]))
  .GlobalEnv$myfm <- myfm # fix a global env bug
  mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 
  #print(mydata)
  model = nlme::lme(fixed = myfm, random = ~ 1|sites, data = mydata, method = 'ML') # no popstr PCs
  #summary(mymodel)
  format_lmm(model, envvar) # output model results
}
print(lmeres)
print(dim(deltap))
print(dim(lmeres))
dimnames(lmeres)[[2]] = c('R2m', 'R2c', 'beta', 'beta_p', 'BIC')

print(dim(deltap))
print(dim(lmeres))

print(lmm_results)
write.csv(lmeres, file = lmm_results)
