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
p_norm = read.csv(p_norm)
deltap = subset(p_norm, select = -chrom_pos)
pop_strc = read.csv(pop_structure,row.names = 1)


# assemble lmm objects
prep_lmm <- function(yy, env_sites, envvar) {
  mydata = cbind(yy, env_sites[c('sites', 'env')])
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

## set up the model
myfm = as.formula(paste0('yy ~ env'))

functions_to_export <- c("lmer", "fixef")  
#nrow(deltap)
lmeres = foreach(ii = 1:3, .combine = 'rbind', .errorhandling = 'remove', .export = functions_to_export ) %dopar% {
  yy = as.numeric(unlist(deltap[ii,]))
  .GlobalEnv$myfm <- myfm # fix a global env bug
  mydata = prep_lmm(yy, env_sites, envvar, pop_strc) 

  model <- lmer('yy ~ env + (1|sites)', data=mydata,  REML = FALSE)
  l_ratio = drop1(model,test="Chisq") #test="Chisq"
  outdt = c(fixef(model)['env'],
            summary(model)$coefficients[, "Std. Error"]['env'],
            anova(model)['env', 'Pr(>F)'],
            fixef(model)['(Intercept)'],
            summary(model)$coefficients[, "Std. Error"]['(Intercept)'],
            BIC(model),
            p_value = l_ratio['env', 'Pr(>F)'])
}

dimnames(lmeres)[[2]] = c('env_value', 'env_eror', 'p_value_env','intercept_value', 'intercept_eror', 'bic', 'lrt')
write.csv(lmeres, file = lmm_results)
