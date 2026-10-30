parameter<-readRDS("~/weihaowang/code/整理/数据生成/lognormal.urt2.rds") #方差用coat算的
corr_matrix = readRDS("~/weihaowang/podasuppsim/simulation/correlation_non_modify.rds")
# parameter<-readRDS("~/weihaowang/code/整理/数据生成/lognormal.urt20.rds") #完全linda的估计方法算urt
# parameter <- readRDS("~/weihaowang/code/整理/数据生成/log.normal.para.rds")
beta0<-parameter$beta0 
sigma2<-(parameter$sigma2)
# sigma2<-1
# sigma2<-readRDS("~/weihaowang/code/整理/数据生成/sigma2_URT.rds")
# sigma2<-runif(m,0.01,5)
# sigma2=rep(1,m)
sig.density.vec <- c(0.05, 0.2)
sample.size.vec <- c(50, 200)
sig.strength.vec <- seq(1.5, 3, length.out = 6)
# sig.strength.vec <- seq(1.05, 2, length.out = 6)
p_value<-rep(0,500)
s1 <- 2; s2 <- 2; s3 <- 6
s <- s1 * s2 * s3 
sample.size <- rep(sample.size.vec, each = s2 * s3)
sig.density <- rep(rep(sig.density.vec, each = s3), s1)
sig.strength <- rep(sig.strength.vec, s1 * s2)
setting <- cbind(sample.size, sig.density, sig.strength)
s <- s1 * s2 * s3 
S1.zero.prop <- 0.3 #zero-inflated的时候有用 现在暂时没用
para <- setting[2, ] #2or6
# mumu=2
m<-500 #taxon数
n<-nnum
gamma<-gammain
set.seed(q) ##############################################
taxon_frac = runif(m, min = 0.1, max = 1)
# beta1 = sample(beta0,m)
# order=match(beta1,beta0)
# sigma2 = sigma2[order]
set.seed(q) ##############################################
X0 <- matrix(exp(rnorm(m * n, beta0, sqrt((sigma2)))), nrow = m)

pi0 <- t(t(X0) / colSums(X0))
pi0.ave <- rowMeans(pi0)
if (any(pi0.ave == 0)) {
  ind <- which(pi0.ave == 0)
  pi0.ave[ind] <- min(pi0.ave[-ind]) / 10
  pi0.ave <- pi0.ave / sum(pi0.ave)
}
mu <- para[3]
ratio=0.01
tmp <- (pi0.ave > ratio)
mu=mumu
mu <- 1* mu * (n <= 50) + 1*mu * (n > 50)
mu.1 <- log(mu * tmp + mu * (ratio / pi0.ave) ^ (1 / 4) * (1 - tmp))
# mu.1<-sample(c(-1,-2,1,2), m,replace = TRUE)
# mu.1<-sample(c(1,3,2),m,replace=TRUE)
# mu.1<-sample(seq(1,2,by=0.5),m,replace=TRUE)
# mu.1<-rep(2,m)
set.seed(q) ##############################################
H <- rbinom(m, 1, gamma)
alpha <- mu.1 * H

# H=rep(0,m)
# abun <- order(pi0.ave, decreasing = TRUE)[1 : floor(m * 0.5)]
# H[sample(abun, m * gamma)] <- 1
# alpha <- mu.1 * H

H=rep(0,m)
abun <- order(pi0.ave, decreasing = TRUE)[1 : floor(m * 0.25)]
H[sample(abun, m * gamma)] <- 1
alpha <- mu.1 * H
# ######## SETTING 1 #########
# u <- rbinom(n, 1, 0.5)
# # while(length(unique(u)) == 1) {
# #   u <- rbinom(n, 1, 0.5)
# # }
# # u[which(u== 0)] <- -1
# # u <- rnorm(n, 0, sqrt(1))
# # u <- runif(n, -1, 1)
# Z <- as.matrix(u)
# beta <- alpha
# colnames(Z)<-'u'

##SETTING WITH confounders 1### I like that
# set.seed(q*1) ##############################################
# # z1 <- rbinom(n, 1, 0.5)
# # z1[which(z1 == 0)] <- -1
# z1 <- rnorm(n, 0, 1)
# set.seed(q*2) ##############################################
# # z2 <- rnorm(n, 0, 1)
# z2 <- runif(n, -1, 1)
# set.seed(q*3) ##############################################
# beta1 <- rnorm(m, 1, 1)
# set.seed(q*4) ##############################################
# beta2 <- rnorm(m, 1, 1)
# set.seed(q*5) ##############################################
# u <- rbinom(n, 1, 1 / ( 1+exp(-0.5 * z1 - 0.5 * z2)))
# confun<-cbind(z1,z2)
# Z <- cbind(u, confun)
# beta <- cbind(alpha, beta1, beta2)

# # #SETTING WITH confounders 2###  ############flag
# set.seed(q) ##############################################
# # z1 <- rbinom(n, 1, 0.5)
# # z1[which(z1 == 0)] <- -1
# z1 <- rnorm(n, 0, 1)
# set.seed(q*2) ##############################################
# # z2 <- rnorm(n, 0, 1)
# z2 <- runif(n, -1, 1)
# set.seed(q*3) ##############################################
# # beta1 <- rnorm(m, 2, 1)
# beta1 <- rnorm(m, 0, 1)
# set.seed(q*4) ##############################################
# # beta2 <- rnorm(m, 2, 1)
# beta2 <- rnorm(m, 0, 1)
# set.seed(q*5) ##############################################
# # u <- rbinom(n, 1, 1 / ( 1+exp(-0.5 * z1 - 0.5 * z2)))
# u = rbinom(n, 1, 0.5)
# confun<-cbind(z1,z2)
# Z <- cbind(u, confun)
# beta <- cbind(alpha, beta1, beta2)

# ##SETTING WITH confounders 3###
# set.seed(q*1) ##############################################
# z1 <- rbinom(n, 1, 0.5)
# z1[which(z1 == 0)] <- -1
# # z1 <- rnorm(n, 0, 1)
# set.seed(q*2) ##############################################
# # z2 <- rnorm(n, 0, 1)
# z2 <- runif(n, -1, 1)
# set.seed(q*3) ##############################################
# beta1 <- rnorm(m, 2, 1)
# set.seed(q*4) ##############################################
# beta2 <- rnorm(m, 2, 1)
# set.seed(q*5) ##############################################
# coef1=runif(m,-1,1)
# set.seed(q*6) ##############################################
# coef2=runif(m,-1,1)
# set.seed(q*7) ##############################################
# u=coef1*z1+coef2*z2+rnorm(n,1,1)
# confun<-cbind(z1,z2)
# Z <- cbind(u, confun)
# beta <- cbind(alpha, beta1, beta2)



#SETTING WITH confounders 4######################################### #test
set.seed(q) ##############################################
u = rbinom(n,1,0.5)
r <- 2
coef <-as.matrix( rep(sqrt(2)/sqrt(r),r))
set.seed(q) ##############################################
zz= as.matrix(u) %*% t(coef) + matrix(rnorm(n * r), n, r)
z1=zz[,1]
z2=zz[,2]
set.seed(q) ##############################################
beta1 <- rnorm(m, 0, 1)
set.seed(q*1234) ##############################################
beta2 <- rnorm(m, 0, 1)
confun<-cbind(z1,z2)
Z <- cbind(u, confun)
beta <- cbind(alpha, beta1,beta2)
# #

set.seed(q) ##############################################
# seq_eff = log(runif(m, min = 0.1, max = 1))
tmp <- beta %*% t(Z) + beta0
X <- matrix(exp(rnorm(m * n, tmp, rep(sqrt((sigma2)), n))), nrow = m)
X_true=X

# # ########### setting A5 for correlated taxa ################
# tmp <- beta %*% t(Z) + beta0
# eigen_corr=eigen(corr_matrix)
# vl = eigen_corr$values
# Sig <- eigen_corr$vectors %*% diag(sqrt(vl)) %*% t(eigen_corr$vectors)
# set.seed(q)
# X <- exp(Sig %*% matrix(rnorm(m * n), nrow = m) + tmp)
# X_true=X
# ##########################################################


X=X*taxon_frac
pi <- t(t(X) / colSums(X))

set.seed(q) ##############################################
N <- 1*rnbinom(n, size = 6.5, mu =8000)
# N <- 1*rnbinom(n, size = 6.5, mu =48000)
# N <- 1*rnbinom(n, size = 5.3, mu =7645)
# N <- 1*rnbinom(n, size = 6.5, mu =5000)
set.seed(q) ##############################################
# N[which(u==1)]=rnbinom(length(which(u==1)), size = 6.5, mu = 48000)
# N[which(u==1)]=rnbinom(length(which(u==1)), size = 6.5, mu = 25000)
set.seed(q) ##############################################
Y <- sapply(1 : n, function(s)rmultinom(1, N[s], pi[, s]))
Y.ld=Y
Z <- as.data.frame(Z)
Z.ld=Z
#Z$u <- factor(Z$u)


Y <- as.data.frame(Y)
colnames(Y) <- rownames(Z) <- paste0('sample', 1 : n)
rownames(Y) <- paste0('taxon', 1 : m)