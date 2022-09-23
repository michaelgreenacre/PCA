### R script #2 for Nature Primer Review on PCA (Greenacre et al., 2022) 
### Required packages ISLR2, PMA and easyCODA (need installing)
### Some additional analyses in this script are not reported in the Review

### functions to compute RSS and explained variances on first two PC
rss <- function(X, decomp, r){
  if(r == 0) X.hat <- matrix(0, nrow(X), ncol(X))
  else{
    pcomp <-  X %*% decomp$v[, 1:r, drop = F]
    X.hat <- pcomp %*% solve(t(pcomp) %*% pcomp) %*% t(pcomp) %*% X
  }
  sum((X - X.hat)^2)
}

var.exp <- function(X, decomp){
  rss0 <- rss(X, decomp, 0)
  rss1 <- rss(X, decomp, 1)
  rss2 <- rss(X, decomp, 2)
  c((rss0 - rss1)/rss0,
    (rss1 - rss2)/rss0)
}

#############################################
### Khan cancer data: non-sparse and sparse #
#############################################

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
khan.X   <- sqrt(1/nrow(khan.x))*scale(khan.x, scale=FALSE)*sqrt(1/ncol(khan.x))
khan.svd <- svd(khan.X)
khan.rpc <- sqrt(nrow(khan.x))*khan.svd$u %*% diag(khan.svd$d)
khan.ccc <- khan.svd$v  
khan.ctr <- rowSums(khan.ccc[,1:2]^2) 
khan.ctr <- khan.ctr/max(khan.ctr)
khan.grey <- round(khan.ctr*80)+1

### Two ways of computing the proportions of explained variance
100*khan.svd$d^2/sum(khan.svd$d^2)
#  [1] 1.537229e+01 1.313125e+01 7.678043e+00 6.628233e+00 5.473905e+00 4.890792e+00 3.913989e+00 3.456066e+00 2.844549e+00 2.244159e+00 2.043301e+00
# [12] 1.898516e+00 1.818857e+00 1.815152e+00 1.550291e+00 1.446684e+00 1.430365e+00 1.163250e+00 1.110249e+00 1.079673e+00 9.813216e-01 9.216966e-01
# ...
# ...
# [56] 2.430630e-01 2.365833e-01 2.286081e-01 2.168095e-01 2.040763e-01 1.646078e-01 9.753056e-03 3.262413e-30
var.exp(khan.X, khan.svd)
# [1] 0.1537229 0.1313125

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
### group means and confidence ellipses
require(easyCODA)
set.seed(123)
CIplot_biv(khan.rpc[,1], khan.rpc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE, shade=TRUE)
set.seed(123)
CIplot_biv(khan.rpc[,1], khan.rpc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, add=TRUE)

### Khan centroids non-sparse
khan.means <-  aggregate(khan.x ~ factor(khan.y), FUN=mean)
khan.wts <- table(khan.y)/length(khan.y) 
rownames(khan.means) <- khan.labs
khan.means <- khan.means[,-1]
khan.centroids <- sweep(as.matrix(khan.means), 2, colMeans(khan.x))
khan.ctrs.X   <- diag(sqrt(khan.wts)) %*% khan.centroids * sqrt(1/ncol(khan.centroids))
khan.ctrs.svd <- svd(khan.ctrs.X)
khan.ctrs.rpc <- khan.ctrs.svd$u %*% diag(khan.ctrs.svd$d)
khan.ctrs.ccc <- khan.ctrs.svd$v  
khan.ctrs.rpc <- -khan.ctrs.rpc
khan.ctrs.ccc <- -khan.ctrs.ccc 

khan.ctrs.ctr <- rowSums(khan.ctrs.ccc[,1:2]^2) 
khan.ctrs.ctr <- khan.ctrs.ctr/max(khan.ctrs.ctr)
khan.ctrs.grey <- round(khan.ctrs.ctr*80)+1

100*khan.ctrs.svd$d^2/sum(khan.ctrs.svd$d^2)
# [1] 4.530604e+01 3.030479e+01 2.438917e+01 3.685463e-29
var.exp(khan.ctrs.X, khan.ctrs.svd)
# [1] 0.4530604 0.3030479

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
### group means and confidence ellipses
require(easyCODA)
set.seed(123)
CIplot_biv(khan.ctrs.spc[,1], khan.ctrs.spc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE, shade=TRUE)
set.seed(123)
CIplot_biv(khan.ctrs.spc[,1], khan.ctrs.spc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, add=TRUE)

##########################
### Khan sparse analyses #
##########################

### Khan sparse
library(PMA)
sparsity <- 7
khan.X <- sqrt(1/nrow(khan.x))*scale(khan.x, scale=FALSE) * sqrt(1/ncol(khan.x))
khan.sparse <- SPC(khan.X, 
                   sumabsv = sparsity, niter = 100, K = 2, trace = FALSE, 
                   orth = TRUE, center = FALSE)
sparse.rpc <- -khan.sparse$u %*% diag(khan.sparse$d) * sqrt(nrow(khan.x))
sparse.ccc <- -khan.sparse$v

### Proportions of explained variance
var.exp(khan.X, khan.sparse)
[1] 0.1357526 0.1238702

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
### group centroids and confidence ellipses
require(easyCODA)
set.seed(123)
CIplot_biv(sparse.rpc[,1], sparse.rpc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE, shade=TRUE)
set.seed(123)
CIplot_biv(sparse.rpc[,1], sparse.rpc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, add=TRUE)

### Number of genes not shrunk to 0 on each dimension
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
khan.ctrs.X  <- diag(sqrt(khan.wts))%*% khan.centroids * sqrt(1/ncol(khan.centroids))
khan.ctrs.sparse <- SPC(khan.ctrs.X, 
                        sumabsv = sparsity, niter = 100, K = 2, trace = FALSE, 
                        orth = TRUE, center = FALSE)
sparse.ctrs.rpc <- diag(sqrt(1/khan.wts)) %*% khan.ctrs.sparse$u %*% diag(khan.ctrs.sparse$d)
sparse.ctrs.ccc <- khan.ctrs.sparse$v
sparse.ctrs.rpc <- -sparse.ctrs.rpc
sparse.ctrs.ccc <- -sparse.ctrs.ccc
### get coordinates of individuals
sparse.ctrs.spc <- khan.x %*% sparse.ctrs.ccc
### variances explained
var.exp(khan.ctrs.X, khan.ctrs.sparse)
# [1] 0.3972605 0.3282668

### FIG 4D
par(mar=c(4.2,4,2,2.5), mgp=c(2,0.7,0), font.lab=2, cex.axis=0.7, las=1)
rescale <- 40
plot(rbind(sparse.ctrs.spc, rescale*sparse.ctrs.ccc), type="n", 
     xlab="Sparse PC1 (39.7%)", ylab="Sparse PC2 (32.8%)", main="", asp=1)
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

### group centroids and confidence ellipses
require(easyCODA)
set.seed(123)
CIplot_biv(sparse.ctrs.spc[,1], sparse.ctrs.spc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, shownames=FALSE, add=TRUE, shade=TRUE)
set.seed(123)
CIplot_biv(sparse.ctrs.spc[,1], sparse.ctrs.spc[,2], group=factor(khan.y), groupnames=khan.labs, 
           groupcols=khan.cols, add=TRUE)
sum(sparse.ctrs.ccc[,1] != 0)  
# 101
sum(sparse.ctrs.ccc[,2] != 0)  
# 115

}
