# Code for comparing the BRT fits to select the set of parameters 
# (tree complexity, learning rate, and number of trees) that minimized 
# cross validation predictive deviance

library(plyr)
library(dplyr)
library(magrittr)

# One by one, read in the BRT fits, and extract residual deviance, number of trees, tree complexity, and learning rate
brt_dev <- adply(list.files("./output/BRT_fits", pattern = ".rds", full.names = T), 1, function(path_name){
  brt_fit <- readRDS(path_name)
  return(c(dev = brt_fit$cv.statistics$deviance.mean,  
           ntrees = brt_fit$n.trees,
           tc = brt_fit$gbm.call$tree.complexity,
           lr = brt_fit$gbm.call$learning.rate))
}, .id = NULL)

# Save the results to a csv to be included as a table in the supplement
write.csv(brt_dev, "output/results/brt_comparison.csv")

# Also create a plot of residual deviance over different tree complexities and learning rates for the supplementary materials
brt_dev = brt_dev[order(brt_dev$tc), ]
pdf("./output/results/brt_comparison.pdf", width = 6, height = 4)
par(mar = c(5, 6, 4, 2) + 0.1)
plot(dev ~ tc, subset(brt_dev, lr == 0.001), col = "red", ylim = c(0.002, 0.0026),
     type = "b", pch = 17, xlab = "Tree complexity", ylab = "",
     main = "Residual deviance by tree complexity and learning rate",
     las = 1, bty = "l")
points(dev ~ tc, subset(brt_dev, lr == 0.005), col = "blue",
       type = "b", pch = 16)
title(ylab = "Residual deviance", line=5)
text(x=8, y=0.0025, "Learning rate = 0.005",font=2, cex=0.9)
text(x=8, y=0.00205, "Learning rate = 0.001",font=2, cex=0.9)
dev.off()

