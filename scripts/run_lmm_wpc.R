library(dplyr)
library(data.table)
library(doParallel)
## 20 threads
cl <- makeCluster(20)
registerDoParallel(cl)

## load data
env_sites = snakemake@input[["env_sites"]]
p_norm = snakemake@input[["p_norm"]]
pop_structure = snakemake@input[["pop_structure"]]
lmm_results = snakemake@input[["lmm_results"]]

env_sites = read.csv(env_sites, row.names = 1, )
p_norm = read.csv(p_norm)
deltap = subset(p_norm, select = -chrom_pos)
pop_strc = read.csv(pop_structure,row.names = 1)

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
            lmesum$BIC)
  return(outdt)
}

## set up the model
envvar = 'env'
myfm = as.formula(paste0('yy ~ ', paste(c(paste0('PC',1:3), envvar), collapse = ' + ')))
print(myfm)

yy = as.numeric(unlist(deltap[1,]))
print(yy)
#.GlobalEnv$myfm <- myfm # fix a global env bug
mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 
print(mydata)
model = nlme::lme(fixed = myfm, random = ~ 1|sites, data = mydata) # no popstr PCs
format_lmm(model, envvar) # output model result


lmeres = foreach(ii = 1:nrow(deltap), .combine = 'rbind', .errorhandling = 'remove') %dopar% {
  yy = as.numeric(unlist(deltap[ii,]))
  print(yy)
  .GlobalEnv$myfm <- myfm # fix a global env bug
  mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 
  #print(mydata)
  model = nlme::lme(fixed = myfm, random = ~ 1|sites, data = mydata) # no popstr PCs
  #summary(mymodel)
  format_lmm(model, envvar) # output model results
}
print(lmeres)
print(dim(deltap))
print(dim(lmeres))
dimnames(lmeres)[[2]] = c('R2m', 'R2c', 'beta', 'beta_p', 'BIC')

print(dim(deltap))
print(dim(lmeres))

write.csv(lmeres, file = lmm_results)
