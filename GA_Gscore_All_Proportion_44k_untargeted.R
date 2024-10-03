#### OPT-train ####
opt.train.gscore.all.proportion.untargeted  = function(X, n , test, th = 10^-4, lambda = NULL,CD_Criteria){
  ## Initial
  N = nrow(X); n0 = length(test); Nc=n0
  cand = test
  if(is.null(lambda)){lambda = 1}
  K = armaMatMult(X,t(X))/ncol(X)
  
  ## Set solution size
  #if(n >= 30){
  #  ss = n
  #}else{
  #  ss = 30
  #}
  ss=10
  
  
  ##solution matrix
  sol = c()
  for (i in 1:ss) {
    a = sort(STRATE_367_44k_untargeted(test,n))  #分層抽樣
    b = matrix(0,N,1)
    b[a,1] = 1
    test_real = setdiff(1:367,test)
    b=b[-test_real,]
    sol = cbind(sol,b)
  }
  random.set = cand[which(sol[,1]==1)]
  number_train = classified_testgroup_44k_367(random.set)
  ## Set th
  
  
  ## Start Genetic Algorithm
  iter = 0; stop = 0
  while (stop == 0) {
    iter = iter+1
    cat(iter, "...", sep="")
    ## G-score
    score = c()
    if(iter==1){max.score = c()}
    for (i in 1:ss) {
      score = c(score, G.score.all(K,test,cand[which(sol[,i]==1)], lambda)[CD_Criteria])
    }
    max.score = c(max.score, max(score)[1])
    
    
    ## Stop criteria
    if(iter > 1000){
      if(max(score)-max.score[iter-100] < th)
      {stop = stop+1}
    }
    
    ## Elite solution
    elite = which(rank(-score)<=5)
    #elite = which(rank(score)<=floor(n/10))
    
    ## Delete solution
    del =sample(seq(ss)[-elite], 2, prob = (1/score[-elite])/sum(1/score[-elite]))
    ## Crossover
    for (i in del) {
      ct = 1
      while (ct==1) {
        chr = sample(seq(ss)[-del], 2)
        pos = sample(seq(Nc-1), 1)
        sol[,i] = c(sol[1:pos, chr[1]], sol[(pos+1):Nc, chr[2]])
        for (j in 1:length(number_train)) {
          if(classified_testgroup_44k_367(cand[which(sol[,i]==1)])[j] !=number_train[j]){
            ct = 0
            a = sort(STRATE_367_44k_untargeted(test,n))  #分層抽樣
            b = matrix(0,N,1)
            b[a,1] = 1
            test_real = setdiff(1:367,test)
            b=b[-test_real,]
            sol[,i] = b
            break;
          }else{ct = 0}
        }  
      }
    }
    ## Mutation
    #max_attempts <- 100  # 設置最大嘗試次数
    
    for (i in seq(ss)) {
      ct <- 1
      attempts <- 0  # 記錄嘗試次数
      max_attempts = 500
      while (ct == 1 && attempts < max_attempts) {
        sol.new <- sol[, i]
        n.sol <- length(which(sol[, i] == 1))
        if(number_train[1]==0){ 
          pos <- c(sample(which(sol[, i] == 0)[-(1:11)], 2))
        }else{pos <- c(sample(which(sol[, i] == 0), 2))}
        number_total <- classified_testgroup_44k_367(c(1:N))
        for (k in 1:length(pos)) {
          a = which(classified_testgroup_44k_367(pos[k])==1)
          if(a==1){
            b = 1
            c = 11
          }else{
            b <- sum(number_total[1:(a - 1)]) + 1
            c <- sum(number_total[1:a])
          }
          if(length(which(cand >= b & cand <= c )[which((sol[which(cand >= b & cand <= c ), i] == 1))])==0){
            break
          }
          d = sample(which(cand >= b & cand <= c )[which((sol[which(cand >= b & cand <= c ), i] == 1))],1)
          pos[k+2] = d
        } 
        sol.new[pos] <- abs(sol.new[pos] - 1)
        for (j in 1:length(number_train)) {
          if (classified_testgroup_44k_367(cand[which(sol.new == 1)])[j] != number_train[j]) {
            ct = 0
            a = sort(STRATE_367_44k_untargeted(test,n))  #分層抽樣
            b = matrix(0,N,1)
            b[a,1] = 1
            test_real = setdiff(1:367,test)
            b=b[-test_real,]
            sol.new = b
            break
          } else {
            ct <- 0
          }
        }
        
        attempts <- attempts + 1  # 更新嘗試次数
        
        
        if (attempts >= max_attempts) {
          cat("Reached maximum attempts for solution", i, ". Trying next solution.\n")
        }
      } 
      # 其他部分保持不變
      if (!(i %in% elite)) {
        sol[, i] <- sol.new
      } else {
        old <- G.score.all(K,test,cand[which(sol[,i]==1)], lambda)[CD_Criteria]
        new <-  G.score.all(K,test,cand[which(sol.new == 1)], lambda)[CD_Criteria]
        if (new > old) {
          sol[, i] <- sol.new
        }
      }
      
    }
    if(stop != 0){cat("\nGenetic Algorithm ended!\n")}
  } ## GA End
  
  
  ## OPT Training Set
  sol = sol[,which(score==max(score))]
  opt.set = cand[which(sol==1)]
  
  fin=list(
    opt = opt.set,
    opt.G.score = max.score,
    rdn = random.set)
  
  return(
    fin
  )
  
}
opt.train.gscore.all.proportion.untargeted  = cmpfun(opt.train.gscore.all.proportion.untargeted) 

