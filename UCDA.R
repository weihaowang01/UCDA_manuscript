

UCDA = function(Z,Y,type="hard-modi",formulaa,seed="1234"){
  ###################################################################################
  ############################# zero-imputation #####################################
  ###################################################################################
  
  O=Y
  # Y=as.matrix(log(O+0.5))
  Y =log(.imputation(O,Z,formulaa))
  # X_true=as.matrix(Z)
  # X=cbind(as.matrix(rep(1,n)),X_true[,1])
  Z=as.matrix(Z)
  m = nrow(Y) # row: taxa
  n = ncol(Y) # col: sample
  X=cbind(as.matrix(rep(1,n)),Z)
  aa=ginv(t(X)%*%X)[2,2]
  # I have count matrix O, and group variables X (contain intercept)
  # take log Y=log(O+0.5) (or imputation)
  # Y: m*n, X: n*p
  
  
  
  ###################################################################################
  ################################## rotation #######################################
  ###################################################################################
  # Step 1: rotation
  # m = nrow(Y) # row: taxa
  # n = ncol(Y) # col: sample
  H1 = diag(rep(1, n))-X%*%solve(t(X)%*%X)%*%t(X) #确认并标记维数
  Y1 = Y%*%X%*%solve(t(X)%*%X) # m*p
  Y2 = Y%*%H1 # m*n
  
  ###################################################################################
  ##################################### QMLE ########################################
  ###################################################################################
  # Step 2.1: 使用Bai et al. 提出的 QMLE 得到Phi和confounders C的估计。先直接使用CATE里的fa.em看看效果
  ## 输入：
  ## Y2
  ## r
  

  
  ## 然后将Y2输入到fa.em中
  # r=4
  # XX=data.frame(X1=X[,2],X2=X[,3]) 
  XX <- data.frame(X[, 2:ncol(X), drop = FALSE])
  colnames(XX) <- paste0("X", 1:(ncol(X)-1))
  formula.r <- reformulate(colnames(XX), response = NULL)
  r=est.confounder.num(formula.r , XX, t(Y), method = "ed",bcv.plot = FALSE) 
  print(c("r=",r))
 
  
  if (r == 0){
    Gamma=NULL
    Factor_coef=NULL
    SSE = (Y2-rowMeans(Y2))%*%H1
    # SSE = Y2- t(Factor_coef) %*% t(Factor)
    SSigma = diag( SSE %*% t(SSE) )/(n-ncol(X)-r)
  }else{
    QMLE.result = fa.em(Y2,r) # r是因子数，要看怎么估计
    
    Gamma = QMLE.result$Gamma # n*r
    Sigma = QMLE.result$Sigma # vector: n*1
    
    # Step 2.2: 使用GLS给出factor的估计，请注意这里估计出来的东西是中心化以后的
    H2 = t(Gamma)%*%solve(diag(Sigma))%*%Gamma
    Factor_coef = solve(H2)%*%t(Gamma)%*%solve(diag(Sigma))%*%t(Y2-rowMeans(Y2)) # r*m
    Factor = Gamma
    # Step 2.3: 方差估计
    # 我们先试一下plug-in的点估计
    # SSE = Y2%*%H1- t(Factor_coef) %*% t(Factor)%*%H1
    SSE = (Y2-rowMeans(Y2))%*%H1- t(Factor_coef) %*% t(Factor)%*%H1
    # SSE = Y2- t(Factor_coef) %*% t(Factor)
    # SSigma = diag( SSE %*% t(SSE) )/(n-ncol(X)-r) %%%%%%%%% 原flag
    SSigma = diag( SSE %*% t(SSE) )/n
    # print(length(SSigma))
    # SSigma = sigma2[keepOTUs]
    
  }
  
  ###################################################################################
  ############################ outlier detection ####################################
  ###################################################################################
  
  # Step 3 代入估计出来的L，然后进行outlier detection 
  # 在这一步的基础上 进行代码补充 进行inference
  # 现在第二列是response，估计出来的L(也就是Factor_coef)是predictor，权重是SSigma
  Y1_2 = as.matrix(Y1[,2]) # m*1 #ncol(X)
  if(is.null(Factor_coef)){
    X_2=NULL
  }else{
  X_2 = t(Factor_coef) # m*r
  }
  # X_2 = beta[keepOTUs,-1]
  # outlier detction model: Y1_2 = B + (1,X_2) * beta + E, B是shift vector, E异方差，predictor中包含一个截距项
  
  # X_IPOD = cbind(matrix(1,nrow(X_2)),X_2)
  X_IPOD = cbind(matrix(1,nrow(Y1_2)),X_2)
  Y_IPOD = Y1_2
  
  if(type=="hard-modi"){
    # H = X_IPOD %*% solve(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    H = X_IPOD %*% ginv(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    weight = SSigma
    res = IPOD_modi(X_IPOD, Y_IPOD, H, weight, method = "hard", TOL = 1e-04, length.out = 50)
  }else if(type=="soft-modi"){
    H = X_IPOD %*% solve(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    weight = SSigma
    res = IPOD_modi(X_IPOD, Y_IPOD, H, weight, method = "soft", TOL = 1e-04, length.out = 50)
  }else if(type=="SCAD"){
    H = X_IPOD %*% solve(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    weight = SSigma
    res = IPOD_modi(X_IPOD, Y_IPOD, H, weight, method = "SCAD", TOL = 1e-04, length.out = 50)
  }else if(type=="MCP"){
    H = X_IPOD %*% solve(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    weight = SSigma
    res = IPOD_modi(X_IPOD, Y_IPOD, H, weight, method = "MCP", TOL = 1e-04, length.out = 50)
  }else if(type=="hard"){
    H = X_IPOD %*% solve(t(X_IPOD)  %*% (X_IPOD)) %*% t(X_IPOD) 
    res = IPODNEW(X_IPOD, Y_IPOD, H, method = "hard", TOL = 1e-04, length.out = 50)
  }else{
    H = X_IPOD %*% solve(t(X_IPOD)  %*% (X_IPOD)) %*% t(X_IPOD) 
    res = IPODNEW(X_IPOD, Y_IPOD, H, method = "soft", TOL = 1e-04, length.out = 50)
  }
  gamma = res$gamma # original data中fit出来的参数
  best_lambda = res$lmdOpt # 得到这个estimate用的lmd
  # res中输出什么？res$gamma是 待估B的值 res$beta0是 待估的tildeC的值
  
  
  ###################################################################################
  ############################ null data 生成数据 ###################################
  ###################################################################################
  # Step 3.1 从Y1_2 = (1,X_2) * beta + E 中生成null data
  residuals = Y_IPOD - X_IPOD%*%as.matrix(res$beta0) - res$gamma
  n_res = length(Y_IPOD)
  # residuals = min(n_res^(1/4), sqrt(n_res/(n_res - sum(res$gamma!= 0) - nrow(as.matrix(res$beta0)))))*residuals ###################### scale
  # print(min(n_res^(1/4), sqrt(n_res/(n_res - sum(res$gamma!= 0) - nrow(as.matrix(res$beta0))))))
  # aa = ginv(t(X_IPOD)%*%X_IPOD)[2,2]
  print(ncol(X))
  set.seed(seed)
  # Y_knockoff = X_IPOD%*%res$beta0 + sample(residuals/sqrt(SSigma), replace = TRUE)*sqrt(SSigma) ###################### flag2
  Y_knockoff = X_IPOD%*%res$beta0 + sample(residuals, replace = TRUE) ###################### flag
  # Y_knockoff = X_IPOD%*%res$beta0 + MASS::mvrnorm(n = 1, mu = rep(0, length(SSigma)), Sigma = aa*diag(SSigma))
  # Y_knockoff = X_IPOD%*%res$beta0 + residuals*sample(c(-1, 1), length(residuals), replace = TRUE)
  # Y_knockoff = X_IPOD%*%res$beta0 + sqrt(aa*SSigma)*rnorm(length(residuals))
  
  
  
  
  # Step 3.2 用Y_knockoff X_IPOD 重新估计beta0和gamma
  if(type=="hard-modi"){
    weight = SSigma
    betaInit = ginv(t(X_IPOD) %*% solve(diag(weight)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(weight)) %*% Y_knockoff
    # betaInit = bic_select_lambda(Y_knockoff,X_IPOD,weight)$best_fit$C_hat
    betaInit = estimate_C(Y_knockoff,X_IPOD,weight,lambda = NULL)$C_hat
    # H = X_IPOD %*% solve(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    H = X_IPOD %*% ginv(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    res_null = IPODFUN_modi(X_IPOD, Y_knockoff, H, best_lambda, weight, betaInit, method = "hard", TOL = 1e-04) #估计值hat{beta}^b
  }else if(type=="soft-modi"){
    weight = SSigma
    betaInit = solve(t(X_IPOD) %*% solve(diag(weight)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(weight)) %*% Y_knockoff
    betaInit = bic_select_lambda(Y_knockoff,X_IPOD,weight)$best_fit$C_hat
    H = X_IPOD %*% solve(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    res_null = IPODFUN_modi(X_IPOD, Y_knockoff, H, best_lambda, weight, betaInit, method = "soft", TOL = 1e-04)
  }else if(type=="SCAD"){
    betaInit = solve(t(X_IPOD) %*% solve(diag(weight)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(weight)) %*% Y_knockoff
    H = X_IPOD %*% solve(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    weight = SSigma
    res_null = IPODFUN_modi(X_IPOD, Y_knockoff, H, best_lambda, weight, betaInit, method = "SCAD", TOL = 1e-04)
  }else if(type=="MCP"){
    betaInit = solve(t(X_IPOD) %*% solve(diag(weight)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(weight)) %*% Y_knockoff
    H = X_IPOD %*% solve(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
    weight = SSigma
    res_null = IPODFUN_modi(X_IPOD, Y_knockoff, H, best_lambda, weight, betaInit, method = "MCP", TOL = 1e-04)
  }else if(type=="hard"){
    betaInit = rlm(Y_knockoff ~ X_IPOD - 1)$coefficients
    H = X_IPOD %*% solve(t(X_IPOD)  %*% (X_IPOD)) %*% t(X_IPOD) 
    res_null = IPODFUNNEW(X_IPOD, Y_knockoff, H, best_lambda,betaInit, method = "hard", TOL = 1e-04 )
  }else{
    betaInit = rlm(Y_knockoff ~ X_IPOD - 1)$coefficients
    H = X_IPOD %*% solve(t(X_IPOD)  %*% (X_IPOD)) %*% t(X_IPOD) 
    res_null = IPODNEW(X_IPOD, Y_IPOD, H, method = "soft", TOL = 1e-04, length.out = 50)
  }
  
  
  
  gamma_knockoff = res_null$gamma
  # H = X_IPOD %*% solve(t(X_IPOD) %*% solve(diag(SSigma)) %*% (X_IPOD)) %*% t(X_IPOD)  %*% solve(diag(SSigma))
  # weight = SSigma
  # res_null = IPOD_modi(X_IPOD, Y_knockoff, H, weight, method = "soft", TOL = 1e-04, length.out = 50) 
  #这里也应该是相同的lmd
  # gamma_knockoff = res_null$gamma #null data中fit出来的参数
  
  
  
  ###################################################################################
  ########################### data-driven bias term #################################
  ###################################################################################
  # Step 3.3 Find the data-driven correction factor gamma
  nu_vec = c()

  B <- 5
  if (is.null(B)) {
    B <- 5
  }
  for (b in 1:B) {
    #print(b)
    set.seed(b)
    # nu = correction_factor_estimation(X, y, coef_corr, hat_sigma, residuals, best_lmd, fdr_value) 
    nu = correction_factor_estimation(X_IPOD, Y_IPOD, res$beta0, res$gamma, SSigma, residuals, best_lambda, 0.05, type) 
    nu_vec = c(nu_vec, nu)
  }
  
  print(nu_vec)
  # fdr_control = max(nu_vec) #计算矫正参数
  fdr_control = median(nu_vec)
  
  # fdr_control = quantile(nu_vec, probs = 0.75)
  
  ###################################################################################
  ######################## 校正null估计量 并inference ###############################
  ###################################################################################
  # Step 3.4 校正null估计量
  # gamma_knockoff = abs(gamma_knockoff)+ fdr_control * hat_sigma * (best_lambda + sqrt(log(m+r)) / sqrt(m)) 
  # gamma_knockoff = abs(gamma_knockoff)+ fdr_control * (best_lambda + sqrt(log(m+r)) / sqrt(m)) ########
  gamma_knockoff = abs(gamma_knockoff)+ fdr_control * (sqrt(log(m)) / sqrt(n))
  # gamma_knockoff = abs(gamma_knockoff)+ fdr_control * (best_lambda + sqrt(log(m+r)) ) 
  threshold = binary_search(gamma, gamma_knockoff, 0.05) #找你的threshold
  reject = which(gamma > threshold)
  return(reject)
}
