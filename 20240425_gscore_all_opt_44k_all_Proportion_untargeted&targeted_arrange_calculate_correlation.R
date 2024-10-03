###########################################################################################################
##########This program is for automatically anylasis optimal training set by Genomic Selection#############
###########################################################################################################
setwd("E:/論文程式")
rm(list=ls())
setwd("C:/Users/User/OneDrive - 國立台灣大學/桌面/論文程式")
### Source functions & installes packaged
library(BGLR)
library(MASS) # g-inverse
library(psych)  #trace
library(Rcpp);library(RcppArmadillo);library(RcppEigen) 
library(MASS)
require(compiler)

# source function
source("function.R")
sourceCpp("matrix_mutiplication.cpp")
### Import Data
# na是和最完整的413比所少的部分
#na = c(8,12,14 , 64 , 69 , 70 , 71 , 73 , 79 , 83,  85 , 93 ,103 ,109, 110 ,112, 119, 129, 144, 154, 155 ,156 ,157 ,158, 162, 167 ,169, 170 ,171 ,175, 176,
#       253, 254 ,281, 295, 327 ,329, 330, 349 ,350, 351, 355 ,369 ,407 ,412 ,413)
# Import data
geno = read.csv("geno_44k_367_11047.csv")
geno = geno[,-1] # marker score matrix
geno_std = scale(geno) # standardlize marker score matrix
pheno = read.csv("pheno_44k_7traits_367.csv")
K = armaMatMult(geno_std,t(geno_std))/ncol(geno_std)
### Parameter Settin
lambda = 1 # lambda = sigma_e^2/sigma_g^2
#
load("44K_testing.rda")

# CD_opt_all_sort 中的每個元素都是一個列表，包含原始列表中餘數相同的子列表

# Create train & ETA & calculate correlation
setwd("E:/論文程式/20240502_44k_correlation_4CD/BRSA")
CD_list = c(1,2,3,4)
nt_list = c(50,75,100,150,200)
cr_list = c("opt")
# tg_list = c("targeted","untargeted")
tg_list = c("untargeted")
K = armaMatMult(geno_std,t(geno_std))/ncol(geno_std)
load("20240507_gscore_all_opt_44k_Proportion_untargeted.rda")

p_true = pheno$BRSA
X_train_all = list()
Y_train_all = list()
K_train_all = list()
ETA_train_all = list()
fm_all = list(list())
yhat_all = list()
correla_all = list()
correla = c()
for (i in 1:30) {
  for (nt in nt_list) {
    for (CD in CD_list) {
      for (tg in tg_list) {
        for (cr in cr_list) {
          # Choose data
          CD_name <- paste("CD_opt_",CD,"_44k_",nt,"_",i,"_Proportion_",tg, sep = "")
          correla =c()
          # Create X_train_all
          x_train_name <- paste("x_train_", CD_name, sep = "")
          xtrain <- eval(parse(text = paste("sort(unlist(CD_opt_all$",CD_name,"$",cr,"))", sep = "")))
          X_train_all[[x_train_name]] <- xtrain
          # Create Y_train_all
          y_train_name <- paste("y_train_",  CD_name, sep = "")
          ytrain <- p_true[eval(parse(text =  paste("sort(unlist(CD_opt_all$",CD_name,"$",cr,"))", sep = "")))]
          Y_train_all[[y_train_name]] <- ytrain
          # K_train_all
          k_train_name <- paste("K_train_", CD_name, sep = "")
          ktrain <- eval(parse(text = paste("K", sep = "")))[xtrain,xtrain]
          K_train_all[[k_train_name]] <- ktrain
          # Create ETA_all
          ETA_train_name <- paste("ETA_train_", CD_name, sep = "")
          ETAtrain <- list(list(K = ktrain ,model="RKHS"))
          ETA_train_all[[ETA_train_name]] <- ETAtrain
          # fm
          fm_name <- paste("fm.",CD_name, sep = "") 
          fm <- BGLR(ytrain,ETA=ETAtrain,verbose = F)
          fm_all[[fm_name]] <- fm
          # y.hat
          yhat_name <- paste("y.",CD_name, sep = "") 
          yhat <- fm$mu+eval(parse(text = paste("K",sep = "")))[Testing[[i]],xtrain]%*%solve(ktrain)%*%fm[["ETA"]][[1]][["u"]]
          yhat_all[[yhat_name]] <- yhat
          # correlation 
          correla_name <- paste("corla_",CD_name,"_",cr , sep = "")
          correla <-cor(yhat,p_true[Testing[[i]]])
          correla_all[[correla_name]] <- correla
        }
      }
    }
  }
}
# save("correla_all",file = "20240425_correlation_BRSA_CD123_ut&t.rda")
save("correla_all",file = "20240502_correlation_BRSA_CD1234_ut.rda")

# new correct untargeted 
load("20240502_correlation_BRSA_CD1234_ut.rda")
correla_all_untargetd <- correla_all
correla_all_sort_untargeted <- split(correla_all_untargetd, rep(1:20))
###### untargeted random repeat 50 times
set.seed(123)
Candidate <- list()
for (i in 1:30) {
  Candidate[[i]] <-setdiff(1:367,Testing[[i]])
}
random_t =list()
random_ut =list()

# # 初始化 random_t 和 random_ut 列表
# random_t_50 <- list()
# random_ut_50 <- list()
# random_t_75 <- list()
# random_ut_75 <- list()
# random_t_100 <- list()
# random_ut_100 <- list()
# random_t_150 <- list()
# random_ut_150 <- list()
# random_t_200 <- list()
# random_ut_200 <- list()
# 
# # 創建 random_t 和 random_ut 列表中的元素
# for (i in 1:30) {
#   random_t_50[[i]] <- vector("list", length = 50)
#   random_ut_50[[i]] <- vector("list", length = 50)
#   random_t_75[[i]] <- vector("list", length = 50)
#   random_ut_75[[i]] <- vector("list", length = 50)
#   random_t_100[[i]] <- vector("list", length = 50)
#   random_ut_100[[i]] <- vector("list", length = 50)
#   random_t_150[[i]] <- vector("list", length = 50)
#   random_ut_150[[i]] <- vector("list", length = 50)
#   random_t_200[[i]] <- vector("list", length = 50)
#   random_ut_200[[i]] <- vector("list", length = 50)
# }
# 
# for (i in 1:50) {
#   for (j in 1:30) {
#     random_t_50[[j]][[i]] = STRATE_367(1:367,Testing[[j]],50)
#     random_ut_50[[j]][[i]] =  STRATE_367_44k_untargeted(Candidate[[j]],50)
#     random_t_75[[j]][[i]] = STRATE_367(1:367,Testing[[j]],75)
#     random_ut_75[[j]][[i]] =  STRATE_367_44k_untargeted(Candidate[[j]],75)
#     random_t_100[[j]][[i]] = STRATE_367(1:367,Testing[[j]],100)
#     random_ut_100[[j]][[i]] =  STRATE_367_44k_untargeted(Candidate[[j]],100)
#     random_t_150[[j]][[i]] = STRATE_367(1:367,Testing[[j]],150)
#     random_ut_150[[j]][[i]] =  STRATE_367_44k_untargeted(Candidate[[j]],150)
#     random_t_200[[j]][[i]] = STRATE_367(1:367,Testing[[j]],200)
#     random_ut_200[[j]][[i]] = STRATE_367_44k_untargeted(Candidate[[j]],200)
#   }
# }
# random_t <- c(random_t_50,random_t_75,random_t_100,random_t_150,random_t_200)
# random_ut <- c(random_ut_50,random_ut_75,random_ut_100,random_ut_150,random_ut_200)
# 
# save(random_t,file = "random_targeted_all.rda")
# save(random_ut,file = "random_untargeted_all.rda")
load("random_targeted_all.rda")
load("random_untargeted_all.rda")

# p_true = pheno$BRSA
# X_train_all = list()
# Y_train_all = list()
# K_train_all = list()
# ETA_train_all = list()
# fm_all = list(list())
# yhat_all = list()
# correla_all = list()
# correla = c()
# for (i in 1:150) {
#   for (j in 1:50) {
#     # Create X_train_all
#     x_train_name <- paste("x_train_",i,"_",j, sep = "")
#     xtrain <- eval(parse(text = paste("sort(unlist(random_ut[[",i,"]][[",j,"]]))", sep = "")))
#     # Create Y_train_all
#     y_train_name <- paste("y_train_", i,"_",j, sep = "")
#     ytrain <- p_true[eval(parse(text = paste("sort(unlist(random_ut[[",i,"]][[",j,"]]))", sep = "")))]
#     # K_train_all
#     k_train_name <- paste("K_train_", i,"_",j, sep = "")
#     ktrain <- eval(parse(text = paste("K", sep = "")))[xtrain,xtrain]
#     # Create ETA_all
#     ETA_train_name <- paste("ETA_train_",i,"_",j, sep = "")
#     ETAtrain <- list(list(K = ktrain ,model="RKHS"))
#     # fm
#     fm_name <- paste("fm.",i,"_",j, sep = "") 
#     fm <- BGLR(ytrain,ETA=ETAtrain,verbose = F)
#     if(i<31){
#       # y.hat
#       yhat_name <- paste("y.",i,"_",j, sep = "") 
#       yhat <- fm$mu+eigenMapMatMult(eigenMapMatMult(eval(parse(text = paste("K",sep = "")))[Testing[[i]],xtrain],solve(ktrain)),fm[["ETA"]][[1]][["u"]])
#       # correlation 
#       correla_name <- paste("corla_",i,"_",j, sep = "")
#       correla <-cor(yhat,p_true[Testing[[i]]])
#       correla_all[[correla_name]] <- correla
#     }else if(i>=31&i<61){
#       # y.hat
#       yhat_name <- paste("y.",i,"_",j, sep = "") 
#       yhat <- fm$mu+eigenMapMatMult(eigenMapMatMult(eval(parse(text = paste("K",sep = "")))[Testing[[i-1*30]],xtrain],solve(ktrain)),fm[["ETA"]][[1]][["u"]])
#       # correlation 
#       correla_name <- paste("corla_",i,"_",j, sep = "")
#       correla <-cor(yhat,p_true[Testing[[i-1*30]]])
#       correla_all[[correla_name]] <- correla
#     }else if(i>=61&i<91){
#       # y.hat
#       yhat_name <- paste("y.",i,"_",j, sep = "") 
#       yhat <- fm$mu+eigenMapMatMult(eigenMapMatMult(eval(parse(text = paste("K",sep = "")))[Testing[[i-2*30]],xtrain],solve(ktrain)),fm[["ETA"]][[1]][["u"]])
#       # correlation 
#       correla_name <- paste("corla_",i,"_",j, sep = "")
#       correla <-cor(yhat,p_true[Testing[[i-2*30]]])
#       correla_all[[correla_name]] <- correla
#     }else if(i>=91&i<121){
#       # y.hat
#       yhat_name <- paste("y.",i,"_",j, sep = "") 
#       yhat <- fm$mu+eigenMapMatMult(eigenMapMatMult(eval(parse(text = paste("K",sep = "")))[Testing[[i-3*30]],xtrain],solve(ktrain)),fm[["ETA"]][[1]][["u"]])
#       # correlation 
#       correla_name <- paste("corla_",i,"_",j, sep = "")
#       correla <-cor(yhat,p_true[Testing[[i-3*30]]])
#       correla_all[[correla_name]] <- correla
#     }else if(i>=121&i<151){
#       # y.hat
#       yhat_name <- paste("y.",i,"_",j, sep = "") 
#       yhat <- fm$mu+eigenMapMatMult(eigenMapMatMult(eval(parse(text = paste("K",sep = "")))[Testing[[i-4*30]],xtrain],solve(ktrain)),fm[["ETA"]][[1]][["u"]])
#       # correlation 
#       correla_name <- paste("corla_",i,"_",j, sep = "")
#       correla <-cor(yhat,p_true[Testing[[i-4*30]]])
#       correla_all[[correla_name]] <- correla
#     }
#   }
# }
# 
# 
# save("correla_all",file = "20240509_random_targeted_all.rda")
# save("correla_all",file = "20240509_random_untargeted_all.rda")

load("20240509_random_targeted_all.rda")
correla_all_t_all = correla_all
load("20240509_random_untargeted_all.rda")
correla_all_ut_all = correla_all

correla_all_t_all_mean <- split(correla_all_t_all, rep(1:ceiling(length(correla_all_t_all)/50), each=50, length.out=length(correla_all_t_all)))
correla_random_t <- c()
for (i in 1:150) {
  correla_random_t[i] <-mean(unlist(correla_all_t_all_mean[[i]]))
}

correla_all_ut_all_mean <- split(correla_all_ut_all, rep(1:ceiling(length(correla_all_ut_all)/50), each=50, length.out=length(correla_all_ut_all)))
correla_random_ut <- c()
for (i in 1:150) {
  correla_random_ut[i] <-mean(unlist(correla_all_ut_all_mean[[i]]))
}

# load("44k_CD0_correlation_BRSA.rda")
# load("44k_CD0_correlation_BRSA_untargeted.rda")
# 將correla_all中的子列表按照 30 的餘數重新整理
load("20240425_correlation_BRSA_CD123_ut&t.rda")
correla_all_sort <- split(correla_all, rep(1:60, length.out = length(correla_all)))

correla_CD = c()
for (i in 1:30) {
  correla_CD <- c(correla_CD,correla_all_sort[[2*i-1]])
}
new_names <- gsub("_opt", "", names(correla_CD))
# 將新的名稱賦值給列表
names(correla_CD) <- new_names

# CD0 
load("44k_CD0_correlation_BRSA.rda")
load("44k_CD0_correlation_BRSA_untargeted.rda")

correla_final =  c( #random n=50 t & ut
                    correla_random_t[1:30],correla_random_ut[1:30],
                    #CD0
                    correla_50_CDv0, corla_Design_all[1:30],
                    #CD1
                    correla_all_sort[[1]], correla_all_sort_untargeted[[1]],
                    #CD2
                    correla_all_sort[[5]], correla_all_sort_untargeted[[2]],
                    #CD3
                    correla_all_sort[[9]], correla_all_sort_untargeted[[3]],
                    #CD4 only untargeted
                                           #correla_all_sort_untargeted[[4]],
                   #random  n=75 t & ut
                   correla_random_t[61:90],correla_random_ut[31:60], 
                   #CD0
                   correla_75_CDv0, corla_Design_all[31:60],
                   #CD1
                   correla_all_sort[[13]], correla_all_sort_untargeted[[5]],
                   #CD2
                   correla_all_sort[[17]], correla_all_sort_untargeted[[6]],
                   #CD3
                   correla_all_sort[[21]], correla_all_sort_untargeted[[7]],
                   #CD4 only untargeted
                                           #correla_all_sort_untargeted[[8]],
                   #random  n=100 t & ut
                   correla_random_t[61:90],correla_random_ut[61:90],
                   #CD0
                   correla_100_CDv0, corla_Design_all[61:90],
                   #CD1
                   correla_all_sort[[25]], correla_all_sort_untargeted[[9]],
                   #CD2
                   correla_all_sort[[29]], correla_all_sort_untargeted[[10]],
                   #CD3
                   correla_all_sort[[33]], correla_all_sort_untargeted[[11]],
                   #CD4 only untargeted
                                           #correla_all_sort_untargeted[[12]],
                   #random n=150 t & ut
                   correla_random_t[91:120],correla_random_ut[91:120],
                   #CD0
                   correla_150_CDv0, corla_Design_all[91:120],
                   #CD1
                   correla_all_sort[[37]], correla_all_sort_untargeted[[13]],
                   #CD2
                   correla_all_sort[[41]], correla_all_sort_untargeted[[14]],
                   #CD3
                   correla_all_sort[[45]], correla_all_sort_untargeted[[15]],
                   #CD4 only untargeted
                                          #correla_all_sort_untargeted[[16]],
                   #random n=200 t & ut
                   correla_random_t[121:150],correla_random_ut[121:150],
                   #CD0
                   correla_200_CDv0, corla_Design_all[121:150],
                   #CD1
                   correla_all_sort[[49]], correla_all_sort_untargeted[[17]],
                   #CD2
                   correla_all_sort[[53]], correla_all_sort_untargeted[[18]],
                   #CD3
                   correla_all_sort[[57]], correla_all_sort_untargeted[[19]]
                   #CD4 only untargeted
                                           #correla_all_sort_untargeted[[20]]
                   )
correla_final <- split(correla_final, rep(1:50, each = 30))
library(ggplot2)
library(tidyr)

# 創建包含所有組合的資料框
combinations <- data.frame(
  CD = as.factor(rep(c("random","random",0,0,1,1,2,2,3,3),5)),
  ntrain = as.factor(rep(c(50, 75, 100, 150, 200),each=10)),
  tg = as.factor(rep(c("targeted", "untargeted","targeted", "untargeted","targeted", "untargeted","targeted", "untargeted","targeted", "untargeted"),5))
)

# 創建空的資料框用於存儲相關係數
correlation_df <- data.frame()

# 循環遍歷 correlation_list 中的每個列表
for (i in 1:length(correla_final)) {
  # 將相關係數提取出來，並轉換為資料框
  corr_data <- data.frame(mean = mean(unlist(correla_final[[i]])),
                          sd = sd(unlist(correla_final[[i]])))
  
  # 添加 CD、ntrain、cr 和 tg 列
  corr_data$CD <- combinations$CD[i]
  corr_data$ntrain <- combinations$ntrain[i]
  corr_data$tg <- combinations$tg[i]
  
  # 將每個列表的相關係數資料合併到 correlation_df 中
  correlation_df <- rbind(correlation_df, corr_data)
}


# # 使用 ggplot2 繪製長條圖
# ggplot(correlation_df, aes(x = factor(ntrain), y = mean, fill = CD)) +
#   geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +  # 畫出長條圖
#   geom_errorbar(aes(ymin = mean, ymax = mean + sd), 
#                 position = position_dodge(width = 0.9), width = 0.25) +  # 添加誤差棒，將底部設置為平均值，頂部設置為平均值加上標準差
#   facet_wrap(~ tg) +  # 按照 'tg' 列分面
#   labs(title = "BRSA")  # 添加圖的標題le = "BRSA")  # 添加圖的標題

# 绘制长条图
ggplot(correlation_df, aes(x = factor(ntrain), y = mean, fill = CD)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +  # 画出长条图，并调整间距
  geom_errorbar(aes(ymin = mean, ymax = mean + sd), 
                position = position_dodge(width = 0.7), width = 0.25) +  # 添加误差棒
  facet_wrap(~ tg) +  # 按照 'tg' 列分面
  labs(title = "BRSA", x = "Training set size", y = "Mean") +  # 添加标题和轴标签
  scale_fill_brewer(palette = "Set3") +  # 更改填充颜色
  theme_minimal()  # 使用更简洁的主题

