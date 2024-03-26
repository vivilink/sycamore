setwd("~/Desktop")
total_pop <- 40000
# sample_size <- 20000
# sample_size <- 10450 + 1954 + 1120 + 3135
sample_size <- 1000
# pop_prev <- 0.0003
pop_prev <- 0.01
# sample_prev <- 0.2
# sample_prev <- (1930 + 835 + 605 ) /sample_size
sample_prev <- 0.3
num_diseased_sample <- sample_size * sample_prev
num_healthy_sample <- sample_size - num_diseased_sample
cutoff <- qnorm(1 - pop_prev)
liability <- rnorm(total_pop,0,1)

print(sum(liability >= cutoff)) # this needs to be equal or larger to num_diseased_sample


pdf("binary_sampling.pdf", width=8, height=4)
par(mfrow=c(1,2))
hist(liability, freq=FALSE, breaks=60, xlab="Liability", main="Random sampling", las=2)
abline(v=cutoff)
hist(c(sample(liability[liability >= cutoff], num_diseased_sample), sample(liability[liability < cutoff], num_healthy_sample)), freq=FALSE, breaks=60, xlab="Liability", main="Ascertained sampling")
dev.off()