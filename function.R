# r-score
r.score=function(Xt,X0,lambda=NULL){
  if(is.null(lambda)){lambda = 1}
  n0 = nrow(X0); nt = nrow(Xt)
  Xt = cbind(rep(1,nt),Xt) #加上平均
  X0 = cbind(rep(1,n0),X0) #加上平均
  A=eigenMapMatMult(t(Xt),solve(eigenMapMatMult(Xt,t(Xt))+lambda*diag(nt)))
  B=eigenMapMatMult(eigenMapMatMult(t(X0),(diag(n0)-matrix(1/n0,n0,n0))),X0)
  C=eigenMapMatMult(eigenMapMatMult(t(A),B),A)
  q12=tr(eigenMapMatMult(eigenMapMatMult(B,A),Xt))
  q1=(n0-1)+tr(B)
  q2=tr(C)+tr(eigenMapMatMult(eigenMapMatMult(t(Xt),C),Xt))
  return(q12/sqrt(q1*q2))
}
# MSPE
MSPE.score = function(Xt,X0,lambda=NULL){
  if(is.null(lambda)){lambda = 1}
  n0 = nrow(X0); nt = nrow(Xt);  p= ncol(X0)
  Xt = cbind(rep(1,nt),Xt) #加上平均
  X0 = cbind(rep(1,n0),X0) #加上平均
  A=eigenMapMatMult(t(Xt),solve(eigenMapMatMult(Xt,t(Xt))+lambda*diag(nt)))
  B=diag(p+1)-eigenMapMatMult(A,Xt)
  C=eigenMapMatMult(X0,A)
  D=eigenMapMatMult(X0,B)
  MSPE_score = 1+1/n0*(tr(eigenMapMatMult(C,t(C)))+tr(eigenMapMatMult(D,t(D))))
  return(MSPE_score)
}
# G-Score version1 with SNP matrix
G.score = function(Xt,X0,lambda=NULL){
  if(is.null(lambda)){lambda = 1}
  n0 = nrow(X0); nt = nrow(Xt);  p= ncol(X0)
  Xt = Xt #加上平均 ,Xt Training
  X0 = X0 #加上平均, X0 Testing
  K1 = eigenMapMatMult(Xt,t(Xt))/p
  K2 = eigenMapMatMult(X0,t(X0))/p
  K12 = eigenMapMatMult(Xt,t(X0))/p
  K21 = eigenMapMatMult(X0,t(Xt))/p
  M = diag(1,nt,nt)-matrix(1/nt,nt,nt)
  K1_inv = solve(K1)
  H = eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(K21,K1_inv),solve(M+lambda*K1_inv)),M)
  G_score = tr(eigenMapMatMult(H,K12))/(tr(K2)*tr(eigenMapMatMult(eigenMapMatMult(H,(K1+lambda*diag(1,nt,nt))),t(H))))^0.5
  return(G_score)
}
# G-Score new version2 with SNP matrix
G.score.new = function(Xt,X0,lambda=NULL){
  if(is.null(lambda)){lambda = 1}
  n0 = nrow(X0); nt = nrow(Xt);  p= ncol(X0)
  Xt = Xt # Xt Training
  X0 = X0 # X0 Testing
  K1 = eigenMapMatMult(Xt,t(Xt))/p
  K2 = eigenMapMatMult(X0,t(X0))/p
  K12 = eigenMapMatMult(Xt,t(X0))/p
  K21 = eigenMapMatMult(X0,t(Xt))/p
  M = diag(1,nt,nt)-matrix(1/nt,nt,nt)
  K1_inv = solve(K1)
  H = eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(K21,K1_inv),solve(M+lambda*K1_inv)),M)
  A = eigenMapMatMult(eigenMapMatMult(H,(K1+lambda*diag(1,nt,nt))),t(H))
  G_score =sum(diag(eigenMapMatMult(H,K12))^2/(diag(K2)*diag(A)))*(1/n0)
  return(G_score)
}

# G-Score  version3 with SNP matrix
G.score.v3 = function(Xt,X0,lambda=NULL){
  if(is.null(lambda)){lambda = 1}
  n0 = nrow(X0); nt = nrow(Xt);  p= ncol(X0)
  Xt = Xt 
  X0 = X0 
  K1 = eigenMapMatMult(Xt,t(Xt))/p
  K2 = eigenMapMatMult(X0,t(X0))/p
  K12 = eigenMapMatMult(Xt,t(X0))/p
  K21 = eigenMapMatMult(X0,t(Xt))/p
  M = diag(1,nt,nt)-matrix(1/nt,nt,nt)
  K1_inv = solve(K1)
  H = eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(K21,K1_inv),solve(M+lambda*K1_inv)),M)
  A = eigenMapMatMult(eigenMapMatMult(H,(K1+lambda*diag(1,nt,nt))),t(H))
  G_score = sum((eigenMapMatMult(H,K12)))^2/(sum(K2)*sum(A))
  return(G_score)
}
# G-Score_all with Kernel
G.score.all = function(Kinship_matrix,testing,training,lambda=NULL){
  if(is.null(lambda)){lambda = 1}
  nt = length(training)
  n0 = length(testing)
  K1 = Kinship_matrix[training,training]
  K2 = Kinship_matrix[testing,testing]
  K12 = Kinship_matrix[training,testing]
  K21 = Kinship_matrix[testing,training]
  M = diag(1,nt,nt)-matrix(1/nt,nt,nt)
  K1_inv = solve(K1)
  H = eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(K21,K1_inv),solve(M+lambda*K1_inv)),M)
  A = eigenMapMatMult(eigenMapMatMult(H,(K1+lambda*diag(1,nt,nt))),t(H))
  G_score_1 = tr(eigenMapMatMult(H,K12))/(tr(K2)*tr(eigenMapMatMult(eigenMapMatMult(H,(K1+lambda*diag(1,nt,nt))),t(H))))^0.5
  G_score_2 = sum(diag(eigenMapMatMult(H,K12))^2/(diag(K2)*diag(A)))*(1/n0)
  G_score_3 = sum((eigenMapMatMult(H,K12)))^2/(sum(K2)*sum(A))
  G_score_0 = sum(diag(K1-lambda*solve(M+lambda*K1_inv))/diag(K1))*(1/n0)
  G_score_all = c( G_score_1, G_score_2, G_score_3,G_score_0)
  return(G_score_all)
}

# 四捨五入
Rounding = function(dat,sum_total){
  a=round(dat)
 if(sum(a) > sum_total){
   if(sum(a)-sum_total >1 ){
     d = sum(a)-sum_total
     max_value <- max(a) # 找到向量 a 中的最大值
     max_value_minus_one <- max_value - d # 將最大值减一
     max_index <- which(a == max_value)[1]  # 找到最大值的索引
     a[max_index] <- max_value_minus_one # 替换原向量中的最大值元素為新的减一後的值
   }else{
   max_value <- max(a) # 找到向量 a 中的最大值
   max_value_minus_one <- max_value - 1 # 將最大值减一
   max_index <- which(a == max_value)[1]  # 找到最大值的索引
   a[max_index] <- max_value_minus_one # 替换原向量中的最大值元素為新的减一後的值
   }
 }else if(sum(a)==sum_total){
   a=round(dat)
 }else{
   if(sum_total-sum(a) > 1  ){
     d = sum_total-sum(a)
     min_value <- min(a) # 找到向量 a 中的最小值
     min_value_plus_one <- min_value + d # 將最小值加一
     min_index <- which(a == min_value)[1] # 找到最小值的索引
     a[min_index] <- min_value_plus_one # 替换原向量中的最大值元素為新的减一後的值
   }else{
   min_value <- min(a) # 找到向量 a 中的最小值
   min_value_plus_one <- min_value + 1 # 將最小值加一
   min_index <- which(a == min_value)[1] # 找到最小值的索引
   a[min_index] <- min_value_plus_one # 替换原向量中的最小小值元素為新的加一後的值
   }
 }
  return(a)
}

# 分層抽樣 改 targeted
STRATE<-   function(Total, test, nt) {
  a <- c()
  cand <- setdiff(Total, test)
  number_total <- classified_testgroup(c(1:413)) #1:413是44k 用其他資料要改
  number_training <- Rounding(classified_testgroup(test) / length(test) * nt, nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  
  return(a)
}
# 分層抽樣 改 targeted Proportion
STRATE_367<- function(Total, test, nt) {
  a <- c()
  cand <- setdiff(Total, test)
  number_total <- classified_testgroup_44k_367(c(1:367)) #1:413是44k 用其他資料要改
  number_training <- Rounding(classified_testgroup_44k_367(test) / length(test) * nt, nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}
# 分層抽樣 改 targeted Deffiency
STRATE_367_Deffiency<- function(Total, test, nt) {
  a <- c()
  cand <- setdiff(Total, test)
  number_total <- classified_testgroup_44k_367(c(1:367)) #1:413是44k 用其他資料要改
  d_eff = c() # 算不同次族群的 deffiency
  d_n = classified_testgroup(Testing[[i]]) # 將Total內的不同次族群分類挑選出來
  for (k in 1:6) {
    if(k==1){
      KA = geno_std[Total[1:d_n[k]],]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
        d_eff[k] = det(KA%*%t(KA))^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }else{
      
      b = sum(d_n[1:(k-1)])+1
      c = sum(d_n[1:k])
      KA = geno_std[Total[b:c],]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
      d_eff[k] = det(KA%*%t(KA))^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }
  }
  number_training <-  Rounding(d_eff/sum(d_eff)*nt,nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}
# 分層抽樣 改 targeted Deffiency
STRATE_367_Deffiency<- function(Total, test, nt) {
  a <- c()
  cand <- setdiff(Total, test)
  number_total <- classified_testgroup_44k_367(c(1:367)) #1:413是44k 用其他資料要改
  d_eff = c() # 算不同次族群的 deffiency
  d_n = classified_testgroup(test) # 將Total內的不同次族群分類挑選出來
  for (k in 1:6) {
    if(k==1){
      KA = geno_std[Total[1:d_n[k]],]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
        d_eff[k] = det(KA%*%t(KA))^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }else{
      b = sum(d_n[1:(k-1)])+1
      c = sum(d_n[1:k])
      KA = geno_std[Total[b:c],]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
      d_eff[k] = det(KA%*%t(KA))^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }
  }
  number_training <-  Rounding(d_eff/sum(d_eff)*nt,nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}
# 分層抽樣  untargeted
STRATE.untarget = function(Total,cand,nt){
  a <- c()
  number_total <- classified_testgroup(Total)
  number_training <- Rounding(classified_testgroup(cand) / length(cand) * nt, nt)# 分層抽樣 改 targeted
  n_subgroups <- length(number_training)
  for (i in 1:n_subgroups) {
    if (i == 1) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
    }else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      a <- c(a, sample(subgroup_samples, number_training[i]))
    }
  }
  return(a)
}

# Targeted Deffiency 44k
  Deffiency_44K_targeted<- function(K,Total, test, nt) {
  a <- c()
  cand <- setdiff(Total, test)
  number_total <- classified_testgroup_44k_367(c(1:367)) #1:413是44k 用其他資料要改
  d_eff = c() # 算不同次族群的 deffiency
  d_n = classified_testgroup(test) # 將Total內的不同次族群分類挑選出來
  for (k in 1:6) {
    if(k==1){
      KA = K[cand[1:d_n[k]],cand[1:d_n[k]]]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
        d_eff[k] = det(KA)^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }else{
      b = sum(d_n[1:(k-1)])+1
      c = sum(d_n[1:k])
      KA = K[cand[1:d_n[k]],cand[1:d_n[k]]]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
        d_eff[k] = det(KA)^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }
  }
  number_training <-  Rounding(d_eff/sum(d_eff)*nt,nt)
  return(number_training)
  }
  # Targeted Deffiency Sorghum
  Deffiency_Sorghum_targeted<- function(K,Total, test, nt) {
    a <- c()
    cand <- setdiff(Total, test)
    number_total <- classified_testgroup_Sorghum(c(1:451)) #1:413是44k 用其他資料要改
    d_eff = c() # 算不同次族群的 deffiency
    d_n = classified_testgroup(test) # 將Total內的不同次族群分類挑選出來
    for (k in 1:4) {
      if(k==1){
        KA = K[cand[1:d_n[k]],cand[1:d_n[k]]]
        if(is.null(dim(KA))==T){
          d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
        }else{
          d_eff[k] = det(KA)^(1/d_n[k]) # 算不同次族群的 deffiency
        }
      }else{
        b = sum(d_n[1:(k-1)])+1
        c = sum(d_n[1:k])
        KA = K[cand[1:d_n[k]],cand[1:d_n[k]]]
        if(is.null(dim(KA))==T){
          d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
        }else{
          d_eff[k] = det(KA)^(1/d_n[k]) # 算不同次族群的 deffiency
        }
      }
    }
    number_training <-  Rounding(d_eff/sum(d_eff)*nt,nt)
    return(number_training)
  }
  # Targeted Deffiency Soybean
  Deffiency_Soybean_targeted<- function(K,Total, test, nt) {
    a <- c()
    cand <- setdiff(Total, test)
    number_total <- classified_testgroup_Soybean(c(1:401)) #1:413是44k 用其他資料要改
    d_eff = c() # 算不同次族群的 deffiency
    d_n = classified_testgroup(test) # 將Total內的不同次族群分類挑選出來
    for (k in 1:5) {
      if(k==1){
        KA = K[cand[1:d_n[k]],cand[1:d_n[k]]]
        if(is.null(dim(KA))==T){
          d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
        }else{
          d_eff[k] = det(KA)^(1/d_n[k]) # 算不同次族群的 deffiency
        }
      }else{
        b = sum(d_n[1:(k-1)])+1
        c = sum(d_n[1:k])
        KA = K[cand[1:d_n[k]],cand[1:d_n[k]]]
        if(is.null(dim(KA))==T){
          d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
        }else{
          d_eff[k] = det(KA)^(1/d_n[k]) # 算不同次族群的 deffiency
        }
      }
    }
    number_training <-  Rounding(d_eff/sum(d_eff)*nt,nt)
    return(number_training)
  }
  # UnTargeted Deffiency 44k
  Deffiency_44K_Untargeted<- function(K,Total,test, nt) {
    a <- c()
    cand <- setdiff(Total, test)
    number_total <- classified_testgroup_44k_367(c(1:367)) #1:413是44k 用其他資料要改
    d_eff = c() # 算不同次族群的 deffiency
    d_n = classified_testgroup(cand) # 將Total內的不同次族群分類挑選出來
    for (k in 1:6) {
      if(k==1){
        KA = K[cand[1:d_n[k]],cand[1:d_n[k]]]
        if(is.null(dim(KA))==T){
          d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
        }else{
          d_eff[k] = det(KA)^(1/d_n[k]) # 算不同次族群的 deffiency
        }
      }else{
        b = sum(d_n[1:(k-1)])+1
        c = sum(d_n[1:k])
        KA = K[cand[1:d_n[k]],cand[1:d_n[k]]]
        if(is.null(dim(KA))==T){
          d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
        }else{
          d_eff[k] = det(KA)^(1/d_n[k]) # 算不同次族群的 deffiency
        }
      }
    }
    number_training <-  Rounding(d_eff/sum(d_eff)*nt,nt)
    return(number_training)
  }
# 分層抽樣  untargeted
STRATE.untarget = function(Total,cand,nt){
  a <- c()
  number_total <- classified_testgroup(Total)
  number_training <- Rounding(classified_testgroup(cand) / length(cand) * nt, nt)# 分層抽樣 改 targeted
  n_subgroups <- length(number_training)
  for (i in 1:n_subgroups) {
    if (i == 1) {
      a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
    }else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      a <- c(a, sample(subgroup_samples, number_training[i]))
    }
  }
  return(a)
}


# 分層抽樣 改soybean targeted Deffiency
STRATE_367_soybean_Deffiency<- function(Total, test, nt) {
  a <- c()
  cand <- setdiff(Total, test)
  number_total <- classified_testgroup_Soybean(c(1:401)) #1:413是44k 用其他資料要改
  d_eff = c() # 算不同次族群的 deffiency
  d_n = classified_testgroup(test) # 將Total內的不同次族群分類挑選出來
  for (k in 1:5) {
    if(k==1){
      KA = geno_std[Total[1:d_n[k]],]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
        d_eff[k] = det(KA%*%t(KA))^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }else{
      b = sum(d_n[1:(k-1)])+1
      c = sum(d_n[1:k])
      KA = geno_std[Total[b:c],]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
        d_eff[k] = det(KA%*%t(KA))^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }
  }
  number_training <-  Rounding(d_eff/sum(d_eff)*nt,nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}

# 分層抽樣 改 soybean targeted Proportion
STRATE_367_soybean<- function(Total, test, nt) {
  a <- c()
  cand <- setdiff(Total, test)
  number_total <- classified_testgroup_Soybean(c(1:401)) #1:401是soybean 用其他資料要改
  number_training <- Rounding(classified_testgroup_Soybean(test) / length(test) * nt, nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}

# 分層抽樣 改 sorghum targeted Proportion
STRATE_367_sorghum<- function(Total, test, nt) {
  a <- c()
  cand <- setdiff(Total, test)
  number_total <- classified_testgroup_Sorghum(c(1:451)) #1:401是soybean 用其他資料要改
  number_training <- Rounding(classified_testgroup_Sorghum(test) / length(test) * nt, nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}

# 分層抽樣 改 sorghum targeted Deffiency
STRATE_367_sorghum_Deffiency<- function(Total, test, nt) {
  a <- c()
  cand <- setdiff(Total, test)
  number_total <- classified_testgroup_Sorghum(c(1:451)) #1:413是44k 用其他資料要改
  d_eff = c() # 算不同次族群的 deffiency
  d_n = classified_testgroup(test) # 將Total內的不同次族群分類挑選出來
  for (k in 1:4) {
    if(k==1){
      KA = geno_std[Total[1:d_n[k]],]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
        d_eff[k] = det(KA%*%t(KA))^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }else{
      b = sum(d_n[1:(k-1)])+1
      c = sum(d_n[1:k])
      KA = geno_std[Total[b:c],]
      if(is.null(dim(KA))==T){
        d_eff[k] = sum(KA^2) # 算不同次族群的 deffiency
      }else{
        d_eff[k] = det(KA%*%t(KA))^(1/d_n[k]) # 算不同次族群的 deffiency
      }
    }
  }
  number_training <-  Rounding(d_eff/sum(d_eff)*nt,nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}


# 分層抽樣 改 44k untargeted Proportion
STRATE_367_44k_untargeted<- function(Total, nt) {
  a <- c()
  cand <- Total
  number_total <- classified_testgroup_44k_367(c(1:367)) #1:367是44K用其他資料要改
  number_training <- Rounding(classified_testgroup_44k_367(cand) / length(cand) * nt, nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}

# 分層抽樣 改 sorghum untargeted Proportion
STRATE_367_sorghum_untargeted<- function(Total, nt) {
  a <- c()
  cand <- Total
  number_total <- classified_testgroup_Sorghum(c(1:451)) #1:451是sorghum 用其他資料要改
  number_training <- Rounding(classified_testgroup_Sorghum(cand) / length(cand) * nt, nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}

# 分層抽樣 改 soybean untargeted(Total-Testing = Candidate) Proportion
STRATE_367_soybean_untargeted<- function(Total, nt) {
  a <- c()
  cand <-  Total
  number_total <- classified_testgroup_Soybean(c(1:401)) #1:401是soybean 用其他資料要改
  number_training <- Rounding(classified_testgroup_Soybean(cand) / length(cand) * nt, nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}

# 分層抽樣 改 soybean untargeted(All Candidate Set) Proportion 
STRATE_367_soybean_untargeted_Total<- function(Total, nt) {
  a <- c()
  cand <- Total
  number_total <- classified_testgroup_Soybean(c(1:401)) #1:401是soybean 用其他資料要改
  number_training <- Rounding(classified_testgroup_Soybean(Total) / length(Total) * nt, nt)
  n_subgroups <- length(number_training)
  
  for (i in 1:n_subgroups) {
    if (i == 1) {
      if (length(cand[which(cand <= number_total[i])]) >= number_training[i]) {
        a <- c(a, sample(cand[which(cand <= number_total[i])], number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, cand[which(cand <= number_total[i])])
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 從個數最多的子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    } else {
      b <- sum(number_total[1:(i - 1)]) + 1
      c <- sum(number_total[1:i])
      subgroup_samples <- cand[which(cand >= b & cand <= c )]
      
      if (length(subgroup_samples) >= number_training[i]) {
        a <- c(a, sample(subgroup_samples, number_training[i]))
      } else {
        # 将剩下的全部样本作为训练集
        a <- c(a, subgroup_samples)
        # 计算需要从比例最高的子群中补充的样本数量
        remaining <- number_training[i] - length(a)
        if (remaining > 0) {
          # 隨機從其他子群中抽取足够数量的样本来补充
          max_group <- which.max(number_total)
          ml <- sum(number_total[1:(max_group - 1)]) + 1
          mh <- sum(number_total[1:max_group])
          a <- c(a, sample(cand[which(cand >= ml & cand <= mh )], remaining, replace = F))
        }
      }
    }
  }
  a <- sort(a)
  a <- a[duplicated(a)==FALSE]
  b = nt-length(a)
  c = sample(setdiff(cand,a),b)
  a = sort(c(a,c))
  return(a)
}

# target method training set proportion 跟 testing set 有關係
classified_testgroup = function(n0){
  n_Aromatic=0;n_Aus=0;n_Indica=0;n_Temj=0;n_Troj=0;n_Admixed=0
  for (i in 1:length(n0)) {
    if(n0[i]>=1 & n0[i]<=14){
      n_Aromatic= n_Aromatic+1
    }else if(n0[i]>=15 & n0[i]<=71){
      n_Aus= n_Aus+1
    }else if(n0[i]>=72 & n0[i]<=158){
      n_Indica= n_Indica+1
    }else if(n0[i]>=159 & n0[i]<=254){
      n_Temj= n_Temj+1
    }else if(n0[i]>=255 & n0[i]<=351){
      n_Troj= n_Troj+1
    }else{
      n_Admixed=n_Admixed+1
    }
  }
  number=c(n_Aromatic,n_Aus,n_Indica,n_Temj,n_Troj,n_Admixed)
  return(number)
}
# 用來分n=367的44k rice
classified_testgroup_44k_367 = function(n0){
  n_Aromatic=0;n_Aus=0;n_Indica=0;n_Temj=0;n_Troj=0;n_Admixed=0
  for (i in 1:length(n0)) {
    if(n0[i]>=1 & n0[i]<=11){
      n_Aromatic= n_Aromatic+1
    }else if(n0[i]>=12 & n0[i]<=64){
      n_Aus= n_Aus+1
    }else if(n0[i]>=65 & n0[i]<=134){
      n_Indica= n_Indica+1
    }else if(n0[i]>=135 & n0[i]<=221){
      n_Temj= n_Temj+1
    }else if(n0[i]>=222 & n0[i]<=310){
      n_Troj= n_Troj+1
    }else{
      n_Admixed=n_Admixed+1
    }
  }
  number=c(n_Aromatic,n_Aus,n_Indica,n_Temj,n_Troj,n_Admixed)
  return(number)
}
# 用來Soybean
classified_testgroup_Soybean = function(n0){
  Pop1=0;Pop2=0;Pop3=0;Pop4=0;Pedigree=0
  for (i in 1:length(n0)) {
    if(n0[i]>=1 & n0[i]<=84){
      Pop1= Pop1+1
    }else if(n0[i]>=85 & n0[i]<=168){
      Pop2= Pop2+1
    }else if(n0[i]>=169 & n0[i]<=250){
      Pop3= Pop3+1
    }else if(n0[i]>=251 & n0[i]<=334){
      Pop4= Pop4+1
    }else{
      Pedigree=Pedigree+1
    }
  }
  number=c(Pop1,Pop2,Pop3,Pop4,Pedigree)
  return(number)
}
# 用來Sorghum
classified_testgroup_Sorghum = function(n0){
  Pop1=0;Pop2=0;Pop3=0;Pop4=0
  for (i in 1:length(n0)) {
    if(n0[i]>=1 & n0[i]<=193){
      Pop1= Pop1+1
    }else if(n0[i]>=194 & n0[i]<=290){
      Pop2= Pop2+1
    }else if(n0[i]>=291 & n0[i]<=411){
      Pop3= Pop3+1
    }else{
      Pop4= Pop4+1
    }
  }
  number=c(Pop1,Pop2,Pop3,Pop4)
  return(number)
}
# 當Candidate用CD比完順序，再分到Subpopulation
CD_rank = function(n0){
  Total = list()
  n_Aromatic=c();n_Aus=c();n_Indica=c();n_Temj=c();n_Troj=c();n_Admixed=c()
  for (i in 1:length(n0)) {
    if(n0[i]>=1 & n0[i]<=14){
      n_Aromatic = c(n_Aromatic,n0[i])
    }else if(n0[i]>=15 & n0[i]<=71){
      n_Aus = c(n_Aus,n0[i])
    }else if(n0[i]>=72 & n0[i]<=158){
      n_Indica = c(n_Indica,n0[i])
    }else if(n0[i]>=159 & n0[i]<=254){
      n_Temj = c(n_Temj,n0[i])
    }else if(n0[i]>=255 & n0[i]<=351){
      n_Troj = c(n_Troj,n0[i])
    }else{
      n_Admixed = c(n_Admixed,n0[i])
    }
  }
  Total[[1]] = n_Aromatic; Total[[2]] = n_Aus; Total[[3]] = n_Indica; Total[[4]] = n_Temj; 
  Total[[5]] = n_Troj; Total[[6]] = n_Admixed; 
  return(Total)
}
# 當Candidate用CD比完順序，再分到Subpopulation
CD_rank_44k_367 = function(n0){
  Total = list()
  n_Aromatic=c();n_Aus=c();n_Indica=c();n_Temj=c();n_Troj=c();n_Admixed=c()
  for (i in 1:length(n0)) {
    if(n0[i]>=1 & n0[i]<=11){
      n_Aromatic = c(n_Aromatic,n0[i])
    }else if(n0[i]>=12 & n0[i]<=64){
      n_Aus = c(n_Aus,n0[i])
    }else if(n0[i]>=65 & n0[i]<=134){
      n_Indica = c(n_Indica,n0[i])
    }else if(n0[i]>=135 & n0[i]<=221){
      n_Temj = c(n_Temj,n0[i])
    }else if(n0[i]>=222 & n0[i]<=310){
      n_Troj = c(n_Troj,n0[i])
    }else{
      n_Admixed = c(n_Admixed,n0[i])
    }
  }
  Total[[1]] = n_Aromatic; Total[[2]] = n_Aus; Total[[3]] = n_Indica; Total[[4]] = n_Temj; 
  Total[[5]] = n_Troj; Total[[6]] = n_Admixed; 
  return(Total)
}
# 當Candidate用CD比完順序，再分到Subpopulation Soybean
CD_rank_Soybean = function(n0){
  Total = list()
  Pop1=c();Pop2=c();Pop3=c();Pop4=c();Pedigree=c()
  for (i in 1:length(n0)) {
    if(n0[i]>=1 & n0[i]<=84){
      Pop1 = c(Pop1,n0[i])
    }else if(n0[i]>=85 & n0[i]<=168){
      Pop2 = c(Pop2,n0[i])
    }else if(n0[i]>=169 & n0[i]<=250){
      Pop3 = c(Pop3,n0[i])
    }else if(n0[i]>=251 & n0[i]<=334){
      Pop4 = c(Pop4,n0[i])
    }else{
      Pedigree = c(Pedigree,n0[i])
    }
  }
  Total[[1]] = Pop1; Total[[2]] = Pop2; Total[[3]] = Pop3; Total[[4]] = Pop4; 
  Total[[5]] = Pedigree
  return(Total)
}
# 當Candidate用CD比完順序，再分到Subpopulation Sorghum
CD_rank_Sorghum = function(n0){
  Total = list()
  Pop1=c();Pop2=c();Pop3=c();Pop4=c()
  for (i in 1:length(n0)) {
    if(n0[i]>=1 & n0[i]<=193){
      Pop1 = c(Pop1,n0[i])
    }else if(n0[i]>=194 & n0[i]<=290){
      Pop2 = c(Pop2,n0[i])
    }else if(n0[i]>=291 & n0[i]<=411){
      Pop3 = c(Pop3,n0[i])
    }else{
      Pop4 = c(Pop4,n0[i])
    }
  }
  Total[[1]] = Pop1; Total[[2]] = Pop2; Total[[3]] = Pop3; Total[[4]] = Pop4; 
  return(Total)
}
