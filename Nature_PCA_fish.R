### R script #3 for Nature Primer Review on PCA (Greenacre et al., 2022) 
### Package easyCODA has all the functions needed, also loading 
### vegan, ca and ellipse

require(easyCODA)

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
fish.ca <- ca(rbind(fish.sums,fish), suprow=7:606)
100*fish.ca$sv^2 / sum(fish.ca$sv^2)
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
