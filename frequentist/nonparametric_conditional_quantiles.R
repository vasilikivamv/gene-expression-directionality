library(tidyverse)
library(gcmr)
library(VineCopula)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(kdecopula)
library(copula)

net1<- read.delim(file=".../.../Network1/input data/net1_expression_data.tsv",
                  header = TRUE, sep = "\t")
gold <- read.delim(file=".../.../Network1/gold standard/DREAM5_NetworkInference_GoldStandard_Network1.tsv",
                  header = FALSE, sep = "\t")
interactions <- gold %>% subset(V3==1) %>% select(V1,V2)

# give a value at j to choose a pair of genes
genes <- data.frame(net1 %>% select_if(used$V1[j]==colnames(net1)),
                    net1 %>% select_if(used$V2[j]==colnames(net1)))
# inspect the data
plot(genes)

# Create the copula 

set.seed(1234)
udat = pobs(genes)

U <- udat[,1]        # marginal U ~ Uniform (0,1)  
V <- udat[,2]        # marginal V ~ Uniform (0,1)

# fits a nonparametric kernel density estimate for a copula
fit <- kdecop(udat)
summary(fit)

# plot the dependence structure
plot(fit)
contour(fit, col = "blue")

# Set the quantile values
p = c(0.01,0.05,0.1,0.9,0.95,0.99)



# Conditional copula V|U=u
# C2|1 (v|u) ---> V|U=u
u2 <- seq(0,1, length.out = length(U))
h21 <- matrix(nrow = length(u2), ncol = length(p))
for (j in 1:length(p)) {
  for (i in 1:length(u2)) {
    h21[i,j] <- hkdecop(c(p[j],u2[i]),fit,cond.var = 1)
  }
  colnames(h21) <- c("1%","5%","10%","90%","95%","99%")
}
# plot
matplot(h21, type = "l", lty = 1, lwd = 2,
        col = (cols <- seq_len(ncol(h21))), ylab = "",
        xlab = substitute(C["2|1"](v~"|"~u)~"as a function of"~
                            v~"for a"~{C} ~"copula"))
legend("bottomright", bty = "n", lwd = 2, col = cols,
       legend = as.expression(lapply(seq_along(p), function(j)
         substitute(p[1] == u2, list(u2 = p[j])))))



# Conditional copula U|V=v
# C1|2 (u|v) ---> U|V=v
u1 <- seq(0,1, length.out = length(U))
# Step 3: Fitting the quantile regression model:
h12 <- matrix(nrow = length(u1), ncol = length(p))
for (j in 1:length(p)) {
  for (i in 1:length(u1)) {
    h12[i,j] <- hkdecop(c(u1[i],p[j]),fit, cond.var = 2)
  }
  colnames(h12) <- c("1%","5%","10%","90%","95%","99%")
}

# plot
matplot(h12, type = "l", lty = 1, lwd = 2,
        col = (cols <- seq_len(ncol(h12))), ylab = "",
        xlab = substitute(C["1|2"](u~"|"~v)~"as a function of"~
                            u~"for a"~{C} ~"copula"))
legend("bottomright", bty = "n", lwd = 2, col = cols,
       legend = as.expression(lapply(seq_along(p), function(j)
         substitute(p[1] == u1, list(u1 = p[j])))))





# Inverse C2|1 (v|u) ---> V|U=u
u2 <- seq(0,1, length.out = length(U))
h21_inv <- matrix(nrow = length(u2), ncol = length(p))
for (j in 1:length(p)) {
  for (i in 1:length(u2)) {
    h21_inv[i,j] <- hkdecop(c(p[j],u2[i]),fit,cond.var = 1, inverse = TRUE)
  }
  colnames(h21_inv) <- c("1%","5%","10%","90%","95%","99%")
}
matplot(h21_inv, type = "l", lty = 1, lwd = 2,
        col = (cols <- seq_len(ncol(h21_inv))), ylab = "",
        xlab = substitute(C["2|1"]^(-1)~(v~"|"~u)~"as a function of"~
                            v~"for a"~{C} ~"copula"))
legend("bottomright", bty = "n", lwd = 2, col = cols,
       legend = as.expression(lapply(seq_along(p), function(j)
         substitute(p[1] == u2, list(u2 = p[j])))))






# # Inverse of C1|2 (u|v) ---> U|V=v
u1 <- seq(0,1, length.out = length(U))
# Step 3: Fitting the quantile regression model:
h12_inv<- matrix(nrow = length(u1), ncol = length(p))
for (j in 1:length(p)) {
  for (i in 1:length(u1)) {
    h12_inv[i,j] <- hkdecop(c(u1[i],p[j]),fit, cond.var = 2, inverse = TRUE)
  }
  colnames(h12_inv) <- c("1%","5%","10%","90%","95%","99%")
}
matplot(h12_inv, type = "l", lty = 1, lwd = 2,
        col = (cols <- seq_len(ncol(h12_inv))), ylab = "",
        xlab = substitute(C["1|2"]^(-1)~(u~"|"~v)~"as a function of"~
                            u~"for a"~{C} ~"copula"))
legend("bottomright", bty = "n", lwd = 2, col = cols,
       legend = as.expression(lapply(seq_along(p), function(j)
         substitute(p[1] == u1, list(u1 = p[j])))))


c(mean(h21_inv[,1]),mean(h21_inv[,2]),mean(h21_inv[,3]),mean(h21_inv[,4]),mean(h21_inv[,5]),mean(h21_inv[,6]))


c(mean(h12_inv[,1]),mean(h12_inv[,2]),mean(h12_inv[,3]),mean(h12_inv[,4]),mean(h12_inv[,5]),mean(h12_inv[,6]))

