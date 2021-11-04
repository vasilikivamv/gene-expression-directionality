# Libraries-----------------------------------------------------------------------------------------
library(tidyverse)
library(VineCopula)
library(ggpubr)
library(copula)
library(VC2copula)

# Data----------------------------------------------------------------------------------------------
net1<- read.delim(file="D:/DREAM5_data/Network1/input data/net1_expression_data.tsv",
                  header = TRUE, sep = "\t")
gold<- read.delim(file="D:/DREAM5_data/Network1/gold standard/DREAM5_NetworkInference_GoldStandard_Network1.tsv",
                  header = FALSE, sep = "\t")
used <- gold %>% subset(V3==1) %>% select(V1,V2)

# selected pair of genes (need to inpout a value for j)
genes <- data.frame(net1 %>% select_if(used$V1[j]==colnames(net1)),
                    net1 %>% select_if(used$V2[j]==colnames(net1)))

# Create the copula-------------------------------------------------------------------------------- 

set.seed(1234)
udat = pobs(genes)

U <- udat[,1]        # marginal U ~ Uniform (0,1)  
V <- udat[,2]        # marginal V ~ Uniform (0,1)


cop <- BiCopSelect(udat[,1], udat[,2], familyset = NA, rotations = TRUE)
summary(cop)
c = BiCop(family = cop$family,par = cop$par,par2 = cop$par2)

# check the surface and contour plot of the copula
plot(c)
contour(c, col = "blue")



# Calculate conditional quantiles-----------------------------------------------------------------


# set the quantile values
p = c(0.01,0.05,0.25,0.5,0.75,0.90,0.95)


# Conditional C2|1 (v|u) ---> V|U=u
u2 <- seq(0,1, length.out = nrow(udat))
h21 <- matrix(nrow = length(u2), ncol = 7)
for (j in 1:length(p)) {
  for (i in 1:length(u2)) {
    h21[i,j] <- BiCopHfunc1(p[j],u2[i],obj = c)
  }
  colnames(h21) <- c("1%","5%","25%","50%","75%","90%","95%")
}

# plot
matplot(h21, type = "l", lty = 1, lwd = 2,
        col = (cols <- seq_len(ncol(h21))), ylab = "",
        xlab = substitute(C["2|1"](v~"|"~u)~"as a function of"~
                            v~"for a"~{C} ~"copula"))
legend("bottomright", bty = "n", lwd = 2, col = cols,
       legend = as.expression(lapply(seq_along(p), function(j)
         substitute(p[1] == u2, list(u2 = p[j])))))




# Conditional C1|2 (u|v) ---> U|V=v
u1 <- seq(0,1, length.out = nrow(udat))
h12 <- matrix(nrow = length(u1), ncol = 7)
for (j in 1:length(p)) {
  for (i in 1:length(u1)) {
    h12[i,j] <- BiCopHfunc2(u1[i],p[j],obj = c)
  }
  colnames(h12) <- c("1%","5%","25%","50%","75%","90%","95%")
}

# plot
matplot(h12, type = "l", lty = 1, lwd = 2,
        col = (cols <- seq_len(ncol(h12))), ylab = "",
        xlab = substitute(C["1|2"](u~"|"~v)~"as a function of"~
                            u~"for a"~{C} ~"copula"))
legend("bottomright", bty = "n", lwd = 2, col = cols,
       legend = as.expression(lapply(seq_along(p), function(j)
         substitute(p[1] == u1, list(u1 = p[j])))))





# Conditional inverse C2|1 (v|u) ---> V|U=u
u2 <- seq(0,1, length.out = nrow(udat))
h21_inv <- matrix(nrow = length(u2), ncol = 7)
for (j in 1:length(p)) {
  for (i in 1:length(u2)) {
    h21_inv[i,j] <- BiCopHinv1(p[j],u2[i],obj = c)
  }
  colnames(h21_inv) <- c("1%","5%","25%","50%","75%","90%","95%")
}

# plot
matplot(h21_inv, type = "l", lty = 1, lwd = 2,
        col = (cols <- seq_len(ncol(h21_inv))), ylab = "",
        xlab = substitute(C["2|1"]^(-1)~(v~"|"~u)~"as a function of"~
                            v~"for a"~{C} ~"copula"))
legend("bottomright", bty = "n", lwd = 2, col = cols,
       legend = as.expression(lapply(seq_along(p), function(j)
         substitute(p[1] == u2, list(u2 = p[j])))))






# Conditional inverse of C1|2 (u|v) ---> U|V=v
u1 <- seq(0,1, length.out = nrow(udat))
h12_inv<- matrix(nrow = length(u1), ncol = 7)
for (j in 1:length(p)) {
  for (i in 1:length(u1)) {
    h12_inv[i,j] <- BiCopHinv2(u1[i],p[j],obj = c)
  }
  colnames(h12_inv) <- c("1%","5%","25%","50%","75%","90%","95%")
}

# plot
matplot(h12_inv, type = "l", lty = 1, lwd = 2,
        col = (cols <- seq_len(ncol(h12_inv))), ylab = "",
        xlab = substitute(C["1|2"]^(-1)~(u~"|"~v)~"as a function of"~
                            u~"for a"~{C} ~"copula"))
legend("bottomright", bty = "n", lwd = 2, col = cols,
       legend = as.expression(lapply(seq_along(p), function(j)
         substitute(p[1] == u1, list(u1 = p[j])))))



