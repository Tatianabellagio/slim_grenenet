# simulate power calculation

n = 500 # individuals
p = 5000 # SNPs for both null and alternative
f = 0.5 # MAF
b.alt = 0.2 # effect size under the alternative hypothesis

## simulate genotypes at 1 snp for n individuals
# n independent trials, 2 possible outcomes, maf probability of sucess 
x = rbinom(n, 2, f) # genotypes at 1 SNP for n ind 

## simulate random phenotype with mean 0 and sd 1
y = scale( rnorm(n) ) # random phenotype normalized to have sample sd=1
hist(y)
# regressing y in x 
summary( lm( y ~ x ) )

fit =lm( y ~ x ) 
plot(x = x, y = y)
abline(fit, col = "red")

## run a linear regression between the genotpyes and the phenotpyes
summary( lm( y ~ x ) )$coeff
## save the standar error of the slope
## kind fo a meassure of the average distance between the 
## observed values and the regression line 
se = summary(lm( y ~ x ) )$coeff[2,2] # pick SE, and assume it stays constant and independent of beta

# simulate the 2 betas 
## beta under null, so mean effect is 0 
b.hat.null = rnorm(p, 0, se) # estimates under null
hist(b.hat.null)

## beta under alt hypo
b.hat.alt = rnorm(p, b.alt, se)  # estimates under alternative
hist(b.hat.alt)


par(mfrow=c(1,2))
# Plot observed densities of z-scores 
plot(NULL, xlim = c(-3,6), ylim = c(0,0.5), xlab = "z", 
     ylab = "density", col = "white") # empty panel for plotting

## the walt stats is basically the beta/se 
lines(density( (b.hat.null/se) ), col = "black", lwd = 2) # Wald statistic for null variants
lines(density( (b.hat.alt/se) ), col = "red", lwd = 2) # Wald statistic for alternative variants
# add theoretical densities for z-scores
x.seq = seq(-3, 6, 0.01)
x.seq
lines(x.seq, dnorm(x.seq, 0, 1), col = "blue", lty = 2) # for null
lines(x.seq, dnorm(x.seq, b.alt/se, 1), col = "orange", lty = 2) # for alternative

# Plot observed densities of z^2 
plot(NULL, xlim = c(0,35), ylim = c(0,1), xlab = expression(z^2), 
     ylab = "density", col = "white") # empty panel for plotting

## this just calculate the chi squared statsitci
# (beta/se)**2
lines(density( (b.hat.null/se)^2 ), col = "black", lwd = 2) # chi-square stat for null variants
lines(density( (b.hat.alt/se)^2 ), col = "red", lwd = 2) # chi-square stat for alternative variants
# Let's add theoretical densities of the chi-square distributions
x.seq = seq(0, 35, 0.01)

## now theoretical 
## dchisq(x.seq, df = 1, ncp = (b.alt/se)^2)
## this will create a chisq distribution for the values in xseq with i debree of 
## freedom and ncp equal = 0 for the null and equals (b.alt/se)^2 for the alt 
lines(x.seq, dchisq(x.seq, df = 1, ncp = 0), col = "blue", lty = 2) # ncp=0 for null
lines(x.seq, dchisq(x.seq, df = 1, ncp = (b.alt/se)^2), col = "orange", lty = 2) # ncp = (beta/se)^2 for alternative
legend("topright", leg = c("NULL obs'd","ALT obs'd","NULL theor","ALT theor"),
       col = c("black","red","blue","orange"), 
       lty = c(1,1,2,2), lwd = c(2,2,1,1) )
# Let's add significance thresholds corresponding to 0.05 and 5e-8
# By definition, the thresholds are always computed under the null.

## this si very important 
q.thresh = qchisq( c(0.05, 5e-8), df = 1, ncp = 0, lower = FALSE)
abline(v = q.thresh, col = c("darkgreen", "springgreen"), lty = 3)

## then we plot the 2 significance levels: the classical 5% and the adjusted 
# 
text( q.thresh+2, c(0.4,0.4), c("P<0.05","P<5e-8") )


## ok so basically in this whole code wha t we see is that, the walt statistic
## is useful to compare the null and alternative hypothesis 
## where the null is that the beta is 0 and the alternative
## is that the beta is 0.2 
## we can see how given a z2 scores obtaines we can reject or not the nullh


q.thresh = qchisq(c(0.05,5e-8), df = 1, ncp = 0, lower = FALSE) # signif. thresholds in chi-square units
pchisq(q.thresh, df = 1, ncp = (b.alt/se)^2, lower = FALSE) # corresponding right tail probabilities

# so based on this we have a 90% prob to detect a variant under a signf level of 0.05 (maf=0.2,b=500. maf=0.5)
# but for a significance level of 5e-8 we onlt have 1.3%power





