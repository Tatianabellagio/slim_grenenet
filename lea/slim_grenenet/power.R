## power calc for my simulations 

## minor allele freq 
f = 0.5
## beta under alt hypothesis 
b.alt = 0.5
# candidate values for n
ns = seq(10, 1000, 10) 


sigma = sqrt(1 - 2*f*(1-f)*b.alt^2) # error sd after SNP effect is accounted for (see next part for explanation)

ses = sigma/sqrt(ns*2*f*(1-f)) # SEs corresponding to each candidate n
q.thresh = qchisq(5e-8, df = 1, ncp = 0, lower = F) # chi-sqr threshold corresponding to alpha = 5e-8
pwr = pchisq(q.thresh, df = 1, ncp=(b.alt/ses)^2, lower=F) # power at alpha = 5e-8 for VECTOR of SE values
plot(ns, pwr, col = "darkgreen", xlab = "n", ylab = "power", 
     main = paste0("QT sd=1; MAF=",f,"; beta=",b.alt), t = "l", lwd = 1.5)
abline(h = 0.9, lty = 2)
