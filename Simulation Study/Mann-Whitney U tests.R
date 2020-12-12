setwd("/Users/jarrettphillips/desktop/new_script/data-selection")

# pop size = 1000

a <- read.table(file.choose(), sep = ",") # 90
b <- read.table(file.choose(), sep = ",") # 4545 
c <- read.table(file.choose(), sep = ",") # 303030

# pop size = 10000

e <- read.table(file.choose(), sep = ",") # 90
f <- read.table(file.choose(), sep = ",") # 4545 
g <- read.table(file.choose(), sep = ",") # 303030

# pop size = 100000

i <- read.table(file.choose(), sep = ",") # 90
j <- read.table(file.choose(), sep = ",") # 4545 
k <- read.table(file.choose(), sep = ",") # 303030


# Mann-Whitney U tests - two-tailed at 5% significance level

wilcox.test(a$V2, e$V2) # W = 1305.5, p-value = 0.01797, significant
wilcox.test(a$V2, i$V2) # W = 1394.5, p-value = 0.002038, significant
wilcox.test(e$V2, i$V2) # W = 1130.5, p-value = 0.3424, not significant

wilcox.test(b$V2, f$V2) # W = 215, p-value = 1.386e-07, significant
wilcox.test(b$V2, j$V2) # W = 303, p-value = 1.347e-05, significant
wilcox.test(f$V2, j$V2) # W = 832.5, p-value = 0.2525, not significant

wilcox.test(c$V2, g$V2) # W = 540.5, p-value = 0.06236, not significant
wilcox.test(c$V2, k$V2) # W = 628, p-value = 0.001252, significant
wilcox.test(g$V2, k$V2) # W = 538, p-value = 0.06774, not significant
