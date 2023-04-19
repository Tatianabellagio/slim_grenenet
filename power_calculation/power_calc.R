n = 500 # individuals
p = 5000 # SNPs for both null and alternative
f = 0.5 # MAF
b.alt = 0.2 # effect size under the alternative hypothesis
x = rbinom(n, 2, f) # genotypes at 1 SNP for n ind 
y = scale( rnorm(n) ) # random phenotype normalized to have sample sd=1
se = summary( lm( y ~ x ) )$coeff[2,2] # pick SE, and assume it stays constant and independent of beta
b.hat.null = rnorm(p, 0, se) # estimates under null
b.hat.alt = rnorm(p, b.alt, se)  # estimates under alternative

par(mfrow=c(1,2))
# Plot observed densities of z-scores 
plot(NULL, xlim = c(-3,6), ylim = c(0,0.5), xlab = "z", 
     ylab = "density", col = "white") # empty panel for plotting
lines(density( (b.hat.null/se) ), col = "black", lwd = 2) # Wald statistic for null variants
lines(density( (b.hat.alt/se) ), col = "red", lwd = 2) # Wald statistic for alternative variants
# add theoretical densities for z-scores
x.seq = seq(-3, 6, 0.01)
lines(x.seq, dnorm(x.seq, 0, 1), col = "blue", lty = 2) # for null
lines(x.seq, dnorm(x.seq, b.alt/se, 1), col = "orange", lty = 2) # for alternative

# Plot observed densities of z^2 
plot(NULL, xlim = c(0,35), ylim = c(0,1), xlab = expression(z^2), 
     ylab = "density", col = "white") # empty panel for plotting
lines(density( (b.hat.null/se)^2 ), col = "black", lwd = 2) # chi-square stat for null variants
lines(density( (b.hat.alt/se)^2 ), col = "red", lwd = 2) # chi-square stat for alternative variants
# Let's add theoretical densities of the chi-square distributions
x.seq = seq(0, 35, 0.01)
lines(x.seq, dchisq(x.seq, df = 1, ncp = 0), col = "blue", lty = 2) # ncp=0 for null
lines(x.seq, dchisq(x.seq, df = 1, ncp = (b.alt/se)^2), col = "orange", lty = 2) # ncp = (beta/se)^2 for alternative
legend("topright", leg = c("NULL obs'd","ALT obs'd","NULL theor","ALT theor"),
       col = c("black","red","blue","orange"), 
       lty = c(1,1,2,2), lwd = c(2,2,1,1) )
# Let's add significance thresholds corresponding to 0.05 and 5e-8
# By definition, the thresholds are always computed under the null.
q.thresh = qchisq( c(0.05, 5e-8), df = 1, ncp = 0, lower = FALSE)
abline(v = q.thresh, col = c("darkgreen", "springgreen"), lty = 3)
text( q.thresh+2, c(0.4,0.4), c("P<0.05","P<5e-8") )
