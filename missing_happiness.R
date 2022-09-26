library(dplyr)
library(tidyr)
library(zoo)
library(ggplot2)

#####################
prop <- 0.1
nrep <- 100
#####################

df <- read.csv("world-happiness-report-2021.csv")
info <- df[,c(1:2)] %>%
  rename(country = 1, region = 2)
df <- df[,c(8:12)] %>% 
  rename(social = 1, life = 2, choice = 3, generosity = 4, corruption = 5) %>%
  mutate(life = as.numeric(gsub(",","",life)))
n <- nrow(df) * ncol(df)

rss <- function(X, A){
  sum((X - A)^2, na.rm = TRUE)
}

mPCA <- function(X, rank, maxit = 100, thresh = 1e-5){
  A <- outer(rep(1, nrow(X)), colMeans(X, na.rm = T))
  score0 <- rss(X, A)
  delta <- Inf
  it <- 0
  while(it < maxit & delta > thresh){
    #impute
    Xtilde <- X
    Xtilde[is.na(X)] <- A[is.na(X)]
    #center
    m <- colMeans(Xtilde)
    Xtilde <- scale(Xtilde, center = m, scale = F)
    #lrma
    SVD <- svd(Xtilde)
    A <- SVD$u[, 1:rank, drop = F] %*% diag(SVD$d[1:rank], rank, rank) %*% t(SVD$v[, 1:rank, drop = F])
    #unscale
    A <- scale(A, center = -m, scale = F)
    #loss
    score <- rss(X, A)
    delta <- abs((score - score0)/score0)
    #cat("iter", it , "rss", score, "delta", delta, "\n")
    score0 <- score
    it <- it + 1
  }
  return(score)
}

mPCAnull <- function(X){
  A <- outer(rep(1, nrow(X)), colMeans(X, na.rm = T))
  score0 <- rss(X, A)
  return(score0)
}

result <- c()
for(rep in 1:nrep){
  missing <- sample(1:n, n*prop)
  df_upd <- as.matrix(df)
  df_upd[missing] <- NA
  s <- apply(df_upd, 2, sd, na.rm = T)
  df_upd <- scale(df_upd, scale = s, center = F) 
  
  score0 <- mPCAnull(df_upd) 
  result <- rbind(result, data.frame(rss = score0, r2 = 0, rank = 0, rep = rep))
  for(rank in 1:4){
    score <- mPCA(df_upd, rank)
    result <- rbind(result, data.frame(rss = score, r2 = 1 - score/score0, rank = rank, rep = rep))
  }
}


resultsum <- result %>% group_by(rep) %>% 
  mutate(vardiff = r2 - lag(r2)) %>% 
  filter(rank > 0) %>% 
  group_by(rank) %>% 
  summarise(mean = mean(vardiff), upper = quantile(vardiff, 0.975), lower = quantile(vardiff, 0.025))
  
ggplot(resultsum, aes(x = rank, y = mean)) + 
  geom_bar(stat="identity", fill = "lightskyblue2", color = "black",
           position=position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = upper, ymax= lower), width=.2,
                position=position_dodge(.9)) +
  ylab("variance explained")+
  xlab("PC dimension")+
  theme_classic()
ggsave("pca.png", width = 3, height = 4)
