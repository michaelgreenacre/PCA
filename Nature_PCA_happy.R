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
table(HAPPY[,2])
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


### Average regional points
happy.region <- aggregate(happy.rpc ~ factor(HAPPY[,2]), FUN=mean)
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
col <- c("blue","red")
region.pch <- 1:10
region.cex <- c(1,rep(0.8,8),1)
shortnames <- c("Social", "Life", "Choices", "Generosity", "Corruption")
rescale <- 1
# pdf(file="Fig1Bnew.pdf", width=9.5, height=7.5, useDingbats=FALSE, family="ArialMT")
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
text(happy.region, labels=rownames(happy.region), cex=0.9, font=2, col=region.col)
text(rescale*happy.ccc, labels=shortnames, col=col[2], font=4, cex=0.8)
text(1.03*supp.reg, labels=rownames(supp.reg), col="black", font=4, cex=0.8)
legend("topleft", legend=rownames(happy.region), col=region.col, pch=region.pch, bty="n", pt.cex=region.cex,
       text.col=region.col, text.font=2, cex=0.75)#, title="Regions", title.col="black")
# dev.off()

### Figure 1C
### Scree plot
par(mar=c(4.5,2,3,1), mgp=c(2,0.7,0), font.lab=2, mfrow=c(1,1))
barplot(happy.perc, ylim=c(0,50), width=0.5, space=0.5, xlim=c(0,5.5), names=c("dim1","dim2","dim3","dim4","dim5"), col="lightblue")

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

