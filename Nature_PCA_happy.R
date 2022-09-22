### R script for Nature Primer Review on PCA (Greenacre et al., 2022) 
### Required packages (easyCODA needs installing)
### Some additional analyses in this script are not reported in the Review

require(easyCODA)

#######################################
### PCA of data from happiness report #
#######################################

### read 2021 happiness report data
HAPPY <- read.csv("world-happiness-report-2021.csv", check.names=FALSE)

### regions
table(happy[,2])
#        Central and Eastern Europe Commonwealth of Independent States 
#                                17                                 12 
#                         East Asia        Latin America and Caribbean 
#                                 6                                 20 
#      Middle East and North Africa              North America and ANZ 
#                                17                                  4 
#                        South Asia                     Southeast Asia 
#                                 7                                  9 
#                Sub-Saharan Africa                     Western Europe 
#                                36                                 21                               21 

### PCA on 5 variables, standardized (PCA function from easyCODA package)
happy <- HAPPY[,8:12]
happy.st <- scale(happy)
rownames(happy.st) <- HAPPY[,1]
happy.pca <- PCA(happy.st, weight=FALSE)

### Percentages of variance
round(100*happy.pca$sv^2 / sum(happy.pca$sv^2),1)
# [1] 47.0 24.5 14.1  9.6  4.9

### Cumulative variance explained for each variable
var.expl <- matrix(0,5,2)
rownames(var.expl) <- colnames(happy)
colnames(var.expl) <- paste("Dim", 1:2, sep="")
for(j in 1:5) {
  foo.lm <- lm(happy.st[,j] ~ happy.pca$rowpcoord[,1])
  var.expl[j,1] <- cor(predict(foo.lm),happy.st[,j])^2
  foo.lm <- lm(happy.st[,j] ~ happy.pca$rowpcoord[,1] + happy.pca$rowpcoord[,2])
  var.expl[j,2] <- cor(predict(foo.lm),happy.st[,j])^2
}
round(var.expl,3)
#                               Dim1  Dim2
# Social support               0.680 0.767
# Healthy life expectancy      0.744 0.816
# Freedom to make life choices 0.583 0.664
# Generosity                   0.000 0.782
# Perceptions of corruption    0.341 0.544

round(colMeans(var.expl),3)
#  Dim1  Dim2 
# 0.470 0.715 

### reverse all coordinates and rescale row standard coordinates to get contribution coordinates
happy.rpc <- -happy.pca$rowpcoord                           # row principal
happy.rcc <- -happy.pca$rowcoord * sqrt(happy.pca$rowmass)  # row contribution
happy.csc <- -happy.pca$colcoord                            # column standard
### higher than average row contributions
happy.ctr <- ( happy.rcc[,1]^2 > 1/nrow(happy) ) | ( happy.rcc[,2]^2 > 1/nrow(happy) )
sum(happy.ctr)
# [1] 82

### Cantrill ladder happiness score and log(GDP)
ladder <- HAPPY[,3]
logGDP <- HAPPY[,7]

### Correlation matrix of all variables
round(cor(cbind(ladder,logGDP,happy)),3)
                             ladder logGDP Social support Healthy life expectancy Freedom to make life choices Generosity Perceptions of corruption
ladder                        1.000  0.790          0.757                   0.768                        0.608     -0.018                    -0.421
logGDP                        0.790  1.000          0.785                   0.859                        0.432     -0.199                    -0.342
Social support                0.757  0.785          1.000                   0.723                        0.483     -0.115                    -0.203
Healthy life expectancy       0.768  0.859          0.723                   1.000                        0.461     -0.162                    -0.364
Freedom to make life choices  0.608  0.432          0.483                   0.461                        1.000      0.169                    -0.401
Generosity                   -0.018 -0.199         -0.115                  -0.162                        0.169      1.000                    -0.164
Perceptions of corruption    -0.421 -0.342         -0.203                  -0.364                       -0.401     -0.164                     1.000

### Correlations of two supplementary variables with PC1 and PC2
round(cor(cbind(ladder, logGDP, happy.rpc[,1:2])), 3)
       ladder logGDP             
ladder  1.000  0.790 0.850 -0.067
logGDP  0.790  1.000 0.818 -0.295
        0.850  0.818 1.000  0.000
       -0.067 -0.295 0.000  1.000

### Regression of happiness on the five indicators and square root of R^2
summary(lm(ladder ~ ., data=happy))
# Coefficients:
#                                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    -2.02782    0.64733  -3.133 0.002102 ** 
# `Social support`                3.49701    0.60713   5.760 4.96e-08 ***
# `Healthy life expectancy`       0.05741    0.01067   5.380 2.98e-07 ***
# `Freedom to make life choices`  1.93916    0.51023   3.801 0.000213 ***
# Generosity                      0.20131    0.32741   0.615 0.539620    
# `Perceptions of corruption`    -0.75876    0.29579  -2.565 0.011343 *  
# Residual standard error: 0.5592 on 143 degrees of freedom
# Multiple R-squared:  0.738,     Adjusted R-squared:  0.7289 

### square root of R-squared
sqrt(0.738)
# [1] 0.8590693

### Regressions of supplementary variables on PCs to get biplot coordinates
supp.reg <- matrix(0, 2, 2)
rownames(supp.reg) <- c("happy","logGDP")
supp.reg[1,] <- lm(scale(ladder) ~ happy.rpc[,1] + happy.rpc[,2])$coefficients[2:3]
supp.reg[2,] <- lm(scale(logGDP) ~ happy.rpc[,1] + happy.rpc[,2])$coefficients[2:3]
summary(lm(scale(ladder) ~ happy.rpc[,1] + happy.rpc[,2]))  # R2 = 0.728, sign. dim1 (p<0.0001)
summary(lm(scale(logGDP) ~ happy.rpc[,1] + happy.rpc[,2]))  # R2 = 0.756, sign, dim1 and dim2 (p<0.0001)

### Percentages of variance again, saved
happy.perc <- 100*happy.pca$sv^2/sum(happy.pca$sv^2)
# [1] 46.964540 24.512744 14.050280  9.566051  4.906384

happy.rpc <- -happy.pca$rowpcoord
happy.rcc <- -happy.pca$rowcoord * sqrt(happy.pca$rowmass)
happy.ccc  <- -happy.pca$colcoord
happy.ctr <- ( happy.rcc[,1]^2 > 1/nrow(happy) ) | ( happy.rcc[,2]^2 > 1/nrow(happy) )
sum(happy.ctr)
# [1] 82

### Average regional points
happy.region <- aggregate(-happy.pca$rowpcoord ~ factor(HAPPY[,2]), FUN=mean)
rownames(happy.region) <- happy.region[,1]
happy.region <- happy.region[,-1]

### Regional colours and 50%-reduced alpha colours
### colorspace or RColorBrewer colours
### (RColorBrewer used in Review)
require(RColorBrewer)
region.col <- c(brewer.pal(8, "Dark2"),"forestgreen","blue" )
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}
region.col.alpha <- add.alpha(region.col, 0.5)
region.num <- as.numeric(factor(HAPPY[,2]))

### Figure 1A
### Example of pdf() and dev.off() commands for saving figure
# pdf(file="Fig1A.pdf", width=8.5, height=6, useDingbats=FALSE, family="ArialMT")
par(mar=c(4.2,4,1,1), mgp=c(2.0,0.7,0), font.lab=2, cex.axis=0.8, las=1)
plot(happy.rpc, type="n", xlab="PCA dimension 1 (47.0%)", ylab="PCA dimension 22 (24.5%)", main="", asp=1)
abline(v=0, h=0, lty=2, col="gray")
text(happy.rpc[!happy.ctr,], labels=happy.pca$rownames[!happy.ctr], cex=0.6, 
     col=region.col.alpha[region.num][!happy.ctr], font=2)
text(happy.rpc[happy.ctr,], labels=happy.pca$rownames[happy.ctr], cex=0.6, 
     col=region.col[region.num][happy.ctr], font=2)
text(happy.region, labels=rownames(happy.region), cex=0.9, font=2, col=region.col)
# dev.off()

### Figure 1B
### Scree plot
par(mar=c(4.5,2,3,1), mgp=c(2,0.7,0), font.lab=2, mfrow=c(1,1))
barplot(happy.perc, ylim=c(0,50), width=0.5, space=0.5, xlim=c(0,5.5), names=c("dim1","dim2","dim3","dim4","dim5"), col="lightblue")

### Figure 1C
col <- c("blue","red")
region.pch <- 1:10
region.cex <- c(1,rep(0.8,8),1)
rescale <- 1
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8, las=1, mfrow=c(1,1))
plot(rbind(happy.rpc, rescale*happy.ccc), type="n", xlab="PC1 (47.0%)", ylab="PC2 (24.5%)", main="", asp=1)
#axis(1)
#axis(2)
#axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 
#                                        2), col = "black", col.ticks = col[2], col.axis = col[2])
#axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 
#                                        2), col = "black", col.ticks = col[2], col.axis = col[2])
abline(v=0, h=0, lty=2, col="gray")
arrows(0,0,0.95*rescale*happy.ccc[,1],0.95*rescale*happy.ccc[,2], col="pink", angle=10, length=0.15, lwd=2)
arrows(0,0,0.95*supp.reg[,1],0.95*supp.reg[,2], col="grey", angle=10, length=0.15, lwd=2)
points(happy.rpc, pch=region.pch[region.num], cex=region.cex[region.num], 
     col=region.col[region.num], font=2)
# text(happy.rpc[happy.ctr,], labels=happy.pca$rownames[happy.ctr], cex=0.6, 
#      col=region.col[region.num][happy.ctr], font=2)
text(happy.region, labels=rownames(happy.region), cex=0.9, font=2, col=region.col)
text(rescale*happy.ccc, labels=colnames(happy), col=col[2], font=4, cex=0.8)
text(1.03*supp.reg, labels=rownames(supp.reg), col="black", font=4, cex=0.8)
legend("topleft", legend=rownames(happy.region), col=region.col, pch=region.pch, bty="n", pt.cex=region.cex,
       text.col=region.col, text.font=2, cex=0.75)#, title="Regions", title.col="black")

### correlations and row and columns sums of squared correlations
happy.cor <- cor(cbind(happy.st, happy.rpc))[1:5,6:10]
### Table 1
round(happy.cor, 3)                                                              
# Social support                0.825 -0.295 0.303  0.183  0.328
# Healthy life expectancy       0.862 -0.269 0.002  0.252 -0.347
# Freedom to make life choices  0.764  0.285 0.178 -0.549 -0.050
# Generosity                   -0.007  0.884 0.380  0.268 -0.038
# Perceptions of corruption    -0.584 -0.451 0.659 -0.091 -0.114
rowSums(happy.cor^2,1)
# Social support  Healthy life expectancy  Freedom to make life choices  Generosity  Perceptions of corruption 
#              1                        1                             1           1                          1 
colSums(happy.cor^2,1)
# 2.3482270 1.2256372 0.7025140 0.4783026 0.2453192 


################################
### Khan non-sparse and sparse #
################################

### Package ISLR2 needs installation
require(ISLR2)
data(Khan)
khan.x <- Khan$xtrain
khan.y <- Khan$ytrain
table(khan.y)
 1  2  3  4 
 8 23 12 20 

### colours and labels definition
greys <- rep("",81)
for(i in 1:81) greys[i] <- paste("grey",i+15,sep="")
greys <- rev(greys)
col <- c("blue","red")
require(RColorBrewer)
khan.cols <- brewer.pal(4, "Dark2")
khan.labs <- c("BL","EW","NB","RM")

### Khan non-sparse
khan.svd <- svd(sqrt(1/nrow(khan.x))*scale(khan.x, scale=FALSE)*sqrt(1/ncol(khan.x)))
khan.rpc <- sqrt(nrow(khan.x))*khan.svd$u %*% diag(khan.svd$d)
khan.ccc <- khan.svd$v  
khan.ctr <- rowSums(khan.ccc[,1:2]^2) 
khan.ctr <- khan.ctr/max(khan.ctr)
khan.grey <- round(khan.ctr*80)+1

100*khan.svd$d^2/sum(khan.svd$d^2)
#  [1] 1.537229e+01 1.313125e+01 7.678043e+00 6.628233e+00 5.473905e+00 4.890792e+00 3.913989e+00 3.456066e+00 2.844549e+00 2.244159e+00 2.043301e+00
# [12] 1.898516e+00 1.818857e+00 1.815152e+00 1.550291e+00 1.446684e+00 1.430365e+00 1.163250e+00 1.110249e+00 1.079673e+00 9.813216e-01 9.216966e-01
# ...
# ...
# [56] 2.430630e-01 2.365833e-01 2.286081e-01 2.168095e-01 2.040763e-01 1.646078e-01 9.753056e-03 3.262413e-30

### FIG 4A 
rescale <- 10
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.7, las=1)
plot(rbind(khan.rpc, rescale*khan.ccc), type="n", xlab="PC1 (15.4%)", ylab="PC2 (13.1%)", main="", asp=1)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
abline(v=0, h=0, lty=2, col="gray")
points(rescale*khan.ccc, col=greys[khan.grey], pch=21, bg=greys[khan.grey], cex=0.5)
### convex hulls and group labels
for(g in 1:4) {
  foo <- khan.rpc[khan.y==g,]
  hpts <- chull(foo)
  hpts <- c(hpts, hpts[1])
  lines(foo[hpts, ], col=khan.cols[g], lty=2)
 }
text(khan.rpc, labels=khan.labs[khan.y], col=khan.cols[khan.y], font=2, cex=0.4)
khan.aggr <- aggregate(khan.rpc ~ factor(khan.y), FUN=mean)
rownames(khan.aggr) <- khan.labs
khan.aggr <- khan.aggr[,-1]
text(khan.aggr, labels=rownames(khan.aggr), col=khan.cols, font=2, cex=0.8)
### confidence ellipses
require(easyCODA)
set.seed(123)
CIplot_biv(khan.rpc[,1], khan.rpc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE, shade=TRUE)
set.seed(123)
CIplot_biv(khan.rpc[,1], khan.rpc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE)

### Khan centroids non-sparse
khan.means <-  aggregate(khan.x ~ factor(khan.y), FUN=mean)
khan.wts <- table(khan.y)/length(khan.y) 
rownames(khan.means) <- khan.labs
khan.means <- khan.means[,-1]
khan.centroids <- sweep(as.matrix(khan.means), 2, colMeans(khan.x))
# khan.centers.foo <- diag(sqrt(khan.wts)) %*% sweep(as.matrix(khan.centers), 2, khan.means)
khan.ctrs.svd <- svd(diag(sqrt(khan.wts)) %*% khan.centroids * sqrt(1/ncol(khan.centers)))
khan.ctrs.rpc <- khan.ctrs.svd$u %*% diag(khan.ctrs.svd$d)
khan.ctrs.ccc <- khan.ctrs.svd$v  
khan.ctrs.rpc <- -khan.ctrs.rpc
khan.ctrs.ccc <- -khan.ctrs.ccc 

khan.ctrs.ctr <- rowSums(khan.ctrs.ccc[,1:2]^2) 
khan.ctrs.ctr <- khan.ctrs.ctr/max(khan.ctrs.ctr)
khan.ctrs.grey <- round(khan.ctrs.ctr*80)+1

100*khan.ctrs.svd$d^2/sum(khan.ctrs.svd$d^2)
# [1] 4.530604e+01 3.030479e+01 2.438917e+01 3.685463e-29

khan.ctrs.spc <- ( scale(khan.x, scale=FALSE)/sqrt(nrow(khan.x)) ) %*% khan.ctrs.svd$v
khan.ctrs.spc <- -khan.ctrs.spc

### FIG 4B 
rescale <- 25
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.7, las=1)
plot(rbind(khan.ctrs.spc, rescale*khan.ctrs.ccc), type="n", xlab="PC1 (45.3%)", ylab="PC2 (30.3%)", main="", asp=1)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
abline(v=0, h=0, lty=2, col="gray")
points(rescale*khan.ctrs.ccc, col=greys[khan.ctrs.grey], pch=21, bg=greys[khan.ctrs.grey], cex=0.5)
### convex hulls
for(g in 1:4) {
  foo <- khan.ctrs.spc[khan.y==g,]
  hpts <- chull(foo)
  hpts <- c(hpts, hpts[1])
  lines(foo[hpts, ], col=khan.cols[g], lty=2)
 }
text(khan.ctrs.spc, labels=khan.labs[khan.y], col=khan.cols[khan.y], font=2, cex=0.4)
khan.ctrs.aggr <- aggregate(khan.ctrs.spc ~ factor(khan.y), FUN=mean)
rownames(khan.ctrs.aggr) <- khan.labs
khan.ctrs.aggr <- khan.ctrs.aggr[,-1]
text(khan.ctrs.aggr, labels=rownames(khan.ctrs.aggr), col=khan.cols, font=2, cex=0.8)
require(easyCODA)
set.seed(123)
CIplot_biv(khan.ctrs.spc[,1], khan.ctrs.spc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE, shade=TRUE)
set.seed(123)
CIplot_biv(khan.ctrs.spc[,1], khan.ctrs.spc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE)

##########################
### Khan sparse analyses #
##########################

### Khan sparse
library(PMA)
sparsity <- 7
khan.sparse <- SPC(sqrt(1/nrow(khan.x))*scale(khan.x, scale=FALSE) * sqrt(1/ncol(khan.x)), 
                   sumabsv = sparsity, niter = 100, K = 2, trace = FALSE, 
                   orth = TRUE, center = FALSE)
sparse.rpc <- -khan.sparse$u %*% diag(khan.sparse$d) * sqrt(nrow(khan.x))
sparse.ccc <- -khan.sparse$v

### percentages of variance using rda() in vegan
### for first dimension and for dimensions 1&2 (not nested)
rda(khan.x ~ sparse.rpc[,1])
#                Inertia Proportion Rank
# Total         997.7739     1.0000     
# Constrained   135.4496     0.1358    1
# Unconstrained 862.3243     0.8642   61
rda(khan.x ~ sparse.rpc[,1:2])
#                Inertia Proportion Rank
# Total         997.7739     1.0000     
# Constrained   259.0440     0.2596    2
# Unconstrained 738.7299     0.7404   60

### Difference is percentage on second axis (not nested)
0.2596-0.1358
[1] 0.1238

### Figure 4C
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.lab=0.9, cex.axis=0.7, las=1)
rescale <- 1
plot(rbind(sparse.rpc, rescale*sparse.ccc), type="n", xlab="Sparse PC1 (13.6%)", ylab="Sparse PC2 (12.4%)", main="", asp=1)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 
                                        2), col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 
                                        2), col = "black", col.ticks = col[2], col.axis = col[2])
abline(v=0, h=0, lty=2, col="gray")
### convex hulls
for(g in 1:4) {
  foo <- sparse.rpc[khan.y==g,]
  hpts <- chull(foo)
  hpts <- c(hpts, hpts[1])
  lines(foo[hpts, ], col=khan.cols[g], lty=2)
 }
text(sparse.rpc, labels=khan.labs[khan.y], col=khan.cols[khan.y], font=2, cex=0.4)
points(rescale*sparse.ccc, col="gray40", pch=21, bg="gray", cex=0.5)
### centroids and confidence ellipses
sparse.aggr <- aggregate(sparse.rpc ~ factor(khan.y), FUN=mean)
rownames(sparse.aggr) <- khan.labs
sparse.aggr <- sparse.aggr[,-1]
text(sparse.aggr, labels=rownames(sparse.aggr), col=khan.cols, font=2, cex=0.8)
set.seed(123)
CIplot_biv(sparse.rpc[,1], sparse.rpc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE, shade=TRUE)
set.seed(123)
CIplot_biv(sparse.rpc[,1], sparse.rpc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE)
### Number of genes shrunk to 0 on each dimension
sum(sparse.ccc[,1] != 0)  
# [1] 103
sum(sparse.ccc[,2] != 0)  
# [1] 84

### Khan centroids sparse
khan.means <-  aggregate(khan.x ~ factor(khan.y), FUN=mean)
rownames(khan.means) <- khan.labs
khan.means <- khan.means[,-1]
khan.centroids <- as.matrix(sweep(khan.means , 2, colMeans(khan.x)))
sparsity <- 7
khan.centroids.sparse <- SPC(diag(khan.wts)%*% scale(khan.centroids, scale=FALSE) * sqrt(1/ncol(khan.centroids)), 
                             sumabsv = sparsity, niter = 100, K = 2, trace = FALSE, 
                             orth = TRUE, center = FALSE)
sparse.ctrs.rpc <- diag(sqrt(1/khan.wts)) %*% khan.centroids.sparse$u %*% diag(khan.centroids.sparse$d)
sparse.ctrs.ccc <- khan.centroids.sparse$v
sparse.ctrs.rpc[,1] <- -sparse.ctrs.rpc[,1]
sparse.ctrs.ccc[,1] <- -sparse.ctrs.ccc[,1]

### percetnages of variance using rda in vegan
rda(khan.centroids ~ sparse.ctrs.rpc[,1])
               Inertia Proportion Rank
Total         350.0980     1.0000     
Constrained   113.3542     0.3238    1
Unconstrained 236.7438     0.6762    2
rda(khan.centroids ~ sparse.ctrs.rpc[,1:2])
              Inertia Proportion Rank
Total         350.098      1.000     
Constrained   225.124      0.643    2
Unconstrained 124.974      0.357    1

### Difference is percentage on second axis (not nested)
0.643-0.324
[1] 0.319

### get coordinates of individuals
sparse.ctrs.spc <- khan.x %*% sparse.ctrs.ccc

### FIG 4D
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.7, las=1)
rescale <- 20
plot(rbind(sparse.ctrs.spc, rescale*sparse.ctrs.ccc), type="n", xlab="Sparse PC1 (32.4%)", ylab="Sparse PC2 (31.9%)", main="", asp=1)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 2), 
     col = "black", col.ticks = col[2], col.axis = col[2])
abline(v=0, h=0, lty=2, col="gray")
### convex hulls
for(g in 1:4) {
  foo <- sparse.ctrs.spc[khan.y==g,]
  hpts <- chull(foo)
  hpts <- c(hpts, hpts[1])
  lines(foo[hpts, ], col=khan.cols[g], lty=2)
 }
points(rescale*sparse.ctrs.ccc, col="gray40", pch=21, bg="gray", cex=0.6)
text(sparse.ctrs.spc, labels=khan.labs[khan.y], col=khan.cols[khan.y], font=2, cex=0.5)

### centroids of displayed points and confidence ellipses
sparse.ctrs.aggr <- aggregate(sparse.ctrs.spc ~ factor(khan.y), FUN=mean)
rownames(sparse.ctrs.aggr) <- khan.labs
sparse.ctrs.aggr <- sparse.ctrs.aggr[,-1]
text(sparse.ctrs.aggr, labels=rownames(sparse.aggr), col=khan.cols, font=2, cex=0.8)
set.seed(123)
CIplot_biv(sparse.ctrs.spc[,1], sparse.ctrs.spc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE, shade=TRUE)
set.seed(123)
CIplot_biv(sparse.ctrs.spc[,1], sparse.ctrs.spc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE)
sum(sparse.ctrs.ccc[,1] != 0)  # 114
sum(sparse.ctrs.ccc[,2] != 0)  # 111


######################################################
### Correspondence analysis of Barents Sea fish data #
######################################################
FISH <- read.csv("BarentsFish.csv", row.names=1)
dim(FISH)
[1] 600  83
fish <- FISH[,2:83]
year <- FISH[,1]
[1] 600  82
table(year)
1999 2000 2001 2002 2003 2004 
  88  107   84   99  100  122 
fish <- as.matrix(fish)
rownames(fish) <- year
year.lab <- 1999:2004

### remove species with no data (do this only once!)
remove <- which(colSums(fish>0) == 0)
fish <- fish[,-remove]
dim(fish)
[1] 600  66

### Permutation test between years using package vegan
require(vegan)
fish.cca <- cca(fish ~ factor(year))
#               Inertia Proportion Rank
# Total         1.32599    1.00000     
# Constrained   0.02297    0.01732    5  <- 1.7% of the variance is between-year
# Unconstrained 1.30303    0.98268   65

set.seed(123)
anova(fish.cca)
# Permutation test for cca under reduced model
# Permutation: free
# Number of permutations: 999
# Model: cca(formula = fish ~ factor(year))
#           Df ChiSquare     F Pr(>F)    
# Model      5   0.02297 2.094  0.001 ***  <- but highly significant (p<0.001)
# Residual 594   1.30303                 

### compute year sums, which will be used as centroids in correspondence analysis
### since CA relativizes the data and uses row sums as weights
fish.sums <- aggregate(fish ~ factor(year), FUN=sum)
rownames(fish.sums) <- fish.sums[,1]
fish.sums <- fish.sums[,-1]

### use ca() function in ca package or CA() function in easyCODA
library(easyCODA)
fish.ca <- ca(rbind(fish.sums,fish), suprow=7:606)
100*fish.ca$sv^2 / sum(fish.ca$sv^2)ç
# [1] 54.437441 25.997621 11.092084  5.997792  2.475062
fish.ca$rowpcoord <- fish.ca$rowcoord %*% diag(fish.ca$sv)
fish.ca$colpcoord <- fish.ca$colcoord %*% diag(fish.ca$sv)
fish.rpc <- fish.ca$rowpcoord[1:6,]
fish.spc <- fish.ca$rowpcoord[7:606,]
fish.ccc <- fish.ca$colcoord * sqrt(fish.ca$colmass)
fish.ctr <- (fish.ccc[,1]^2 + fish.ccc[,2]^2)/2  > 1/ncol(fish)
rescale <- 1
col <- c("blue","red")
library(RColorBrewer)
fish.col <- brewer.pal(6, "Dark2")

### FIG 5
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.8, las=1)
plot(rbind(fish.rpc, rescale*fish.ccc), type="n", xlab="CA dimension 1 (54.4%)", 
     ylab="CA dimension 2 (26.0%)", main="", asp=1)
axis(1)
axis(2)
axis(3, at = axTicks(3), labels = round(axTicks(3)/rescale, 
                                        2), col = "black", col.ticks = col[2], col.axis = col[2])
axis(4, at = axTicks(4), labels = round(axTicks(4)/rescale, 
                                        2), col = "black", col.ticks = col[2], col.axis = col[2])
abline(v=0, h=0, lty=2, col="gray")
arrows(0,0,0.95*rescale*fish.ccc[fish.ctr,1],0.95*rescale*fish.ccc[fish.ctr,2],col="pink",
       angle=10, length=0.1)
# text(fish.rpc, labels=c("1999","2000","2001","2002","2003","2004"), col=fish.col, font=2, cex=0.9)
text(rescale*fish.ccc[fish.ctr,], labels=colnames(fish)[fish.ctr], col="red", font=4, cex=0.7)

# bootstrap to get confidence ellipses
nboot <- 1000
fish.spc2 <- matrix(0,6*nboot,2)
set.seed(123)
for(iboot in 1:nboot) {
  foo <- sample(1:600, replace=TRUE)
  fish.boot <- fish[foo,]
  year.lab.boot <- year[foo]
  fish.year.boot <- aggregate(fish.boot ~ factor(year.lab.boot-1998), FUN=sum)
  fish.year.boot <- fish.year.boot[,-1]
  fish.year.boot.pro <- fish.year.boot/rowSums(fish.year.boot)
  fish.spc2[(iboot*6-5):(iboot*6),] <- (as.matrix(fish.year.boot.pro) %*% fish.ca$colcoord)[,1:2]
}
require(ellipse)
set.seed(123)
CIplot_biv(fish.spc2[,1], fish.spc2[,2], ellipse=-1, group=rep(1:6,1000), groupcols=fish.col, 
           add=TRUE, shade=TRUE, shownames=FALSE)
set.seed(123)
CIplot_biv(fish.spc2[,1], fish.spc2[,2], ellipse=-1, group=rep(1:6,1000), groupcols=fish.col, 
           add=TRUE, groupnames=year.lab)

}
