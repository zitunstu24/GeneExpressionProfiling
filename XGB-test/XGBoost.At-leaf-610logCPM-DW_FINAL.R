# Goal: Use Arabidopsis 610 N-DEGs to predict NUE
# Input: gene expression (logCPM) and trait (biomass)
# Output: model performance (correlation between actural biomass and predicted biomass) and
#         feature importance ranking
# Conclusion
# Author: Chia-Yi Cheng 
# Last updated: 2020-12-02

# 0. environment ----------------------------------------------------------------------------------------

library(xgboost)
library(data.table)
library(mlr)
library(ggpubr)
library(ggplot2)

rm(list = setdiff(ls(), lsf.str()))

# 1. Load data -----------------------------------------------------------------------------------

data<-read.table(file = "./ML-input/NutriNet-At-veg.610logCPM-dryweight.input-for-machine-learning.tsv",header = T) 

## Remove unwanted features
data=data[,-c(2:4)] 
data$trait = 1000*(data$trait)

# 2. XGBoost -------------------------------------------------------------------------------------

## Setting parameters

### About the data set
n=18 # number of genotypes
c=6  # number of samples per genotype

### About the run
k=100 # number of iteration
jj=k-1

### Hyperparameters
r=50 # number of rounds
colsample=0.33
eta=0.075 #0.075 #0.1
num_parallel_tree=1
subsample=0.25

### About the output structure
## Need to reset the following otherwise you'll get 'out of subscript' error
y=0
obs = matrix(nrow = k*n, ncol=c)
pred = matrix(nrow = k*n, ncol=c)
rmse = vector()
train.rmse = vector()
impt.out=vector('list',k*n)

for (i in 1:n){
  for (j in 123:(123+jj)) {
    
    test.index <-c(i,i+18,i+36,i+54,i+72,i+90)
    testing<-data[test.index,]
    training<-data[-test.index,]
    
    #convert data frame to data table
    setDT(training) 
    setDT(testing)
    
    #using one hard encoding 
    train.trait <- training$trait 
    test.trait <- testing$trait
    new_training <- model.matrix(~.+0,data = training[,-c("trait"),with=F]) 
    new_testing <- model.matrix(~.+0,data = testing[,-c("trait"),with=F])
    
    #preparing matrix 
    dtrain <- xgb.DMatrix(data = new_training,label = train.trait) 
    dtest <- xgb.DMatrix(data = new_testing,label=test.trait)
    watchlist <- list(train=dtrain, test=dtest)
    
    #user defined evaluation metric
    cor_metric <- function(preds, dtrain) {
      labels <- getinfo(dtrain, "label")
      cor <- cor(preds,labels)
      list(metric = "cor", value = cor)}
    
    params <- list(booster = "gbtree", 
                   objective = "reg:linear", 
                   eta= eta,
                   gamma= 50,  
                   max_depth=6,
                   min_child_weight=1, 
                   eval_metric=cor_metric,
                   subsample=subsample, 
                   colsample_bytree=colsample,
                   num_parallel_tree=num_parallel_tree) 
    
    set.seed(j)
    bst.val<-xgb.train( params = params, 
                        data = dtrain, 
                        nrounds = r,
                        nfold = 5, 
                        showsd = T, 
                        stratified = T, 
                        print_every_n = 20, 
                        early_stop_round = 5, 
                        watchlist = watchlist,
                        maximize = F,
                        verbose = 0)
    
    y=y+1
    
    pred[y,1:c]<- predict(bst.val, dtest)
    obs[y,1:c]<-test.trait
    
    rmse[y]<- as.numeric(bst.val$evaluation_log[r,3])
    train.rmse[y]<-as.numeric(bst.val$evaluation_log[r,2])
    
    # extract important features
    importance_matrix <- xgb.importance(model = bst.val)
    impt.out[[y]]<-paste(importance_matrix$Feature,importance_matrix$Gain, sep = ",")
  }}

save(data,obs, pred, impt.out, rmse, train.rmse,k,n,jj,r,c,colsample,eta,num_parallel_tree,subsample,
     file="./Arabidopsis/XGBoost.At-leaf-610logCPM-DW-r50.RData")

## Organize cor

load(file = "XGBoost.At-leaf-610logCPM-DW-r50.RData")

n=18

pred.mat = matrix(pred, nrow=k)
obs.mat = matrix(obs, nrow=k)

# For each accession, calculate COR for each iteration 

cor.mat=matrix(nrow = jj+1,ncol = n)

for (i in 1:n){
  for (j in 1:(jj+1)){
    O=c(obs.mat[j,c(i,i+n,i+2*n,i+3*n,i+4*n,i+5*n)])
    P=c(pred.mat[j,c(i,i+n,i+2*n,i+3*n,i+4*n,i+5*n)])
    cor.mat[j,i]=cor(P,O)
  }}

COR=vector()
COR=colMeans(cor.mat)

mean(COR)   # 0.6811985

## Plot the output
apply(cor.mat,2,sd)
sd=apply(cor.mat,2,sd)
upper<-COR+sd
lower<-COR-sd

genotype<-factor(unlist(strsplit(rownames(data)[1:18],"_",fixed=T
))[seq(1,by=3,54)])

df <- data.frame(genotype=genotype,
                 Correlation=COR)

par(mar=c(5,6,4,2)+0.1)
par(mgp=c(5,6,4))

ggplot(data=df, aes(x=genotype, y=Correlation)) +
  geom_bar(stat="identity",fill="#af8dc3",
           width =0.4,position=position_dodge())+
  geom_errorbar(aes(ymin=Correlation, ymax=upper), width=.2,colour="#af8dc3",
                position=position_dodge(.9)) +
  theme(axis.text.x = element_text(angle = 90,vjust=1))+
  theme_minimal()

# 3. Extract importantgenes  ----------------------------------------------------------------------

library(tidyverse)
library(data.table)

load(file="XGBoost.At-leaf-610logCPM-DW-r50.RData")

annotation<-read.csv("C:/Users/chiayi/Desktop/genome.fasta/Arabidopsis/Arabidopsis-gene.symbol-description-202005.tsv",
                     header = F,sep = "\t")
names(annotation)=c("Gene","Symbol","Description")

s=18 # number of genotypes, n
t=100 # number of iteration, jj+1
impt=vector("list",s)

# Convert list to data frame

for (i in 1:s){
  m=(i-1)*t+1
  n=t*i
  
  d2 <- as.data.frame(do.call(rbind,flatten(impt.out[m:n]))) %>% 
    separate(.,V1, into =c("Gene","Importance"), sep=",")
  
  ## convert importance from character to numeric
  d2$Importance=as.numeric(d2$Importance)
  
  ## convert data frame to data table for easy calculation
  setDT(d2)
  d2[,sum(Importance),by=Gene]
  output=d2[,.(SUM=sum(Importance)),by=Gene][order(-SUM)]
  ## calculate impt
  impt[[i]]=paste(output$Gene,output$SUM, sep = ",")
  }

# Combine the impt
d3 <- as.data.frame(do.call(rbind,flatten(impt))) %>% 
  separate(.,V1, into =c("Gene","Importance"), sep=",")

d3$Importance=as.numeric(d3$Importance)

# Calculate composite score of importance
setDT(d3)
d3[,sum(Importance),by=Gene]
output=d3[,.(SUM=sum(Importance)),by=Gene][order(-SUM)]
nrow(output) #610

output2=merge(output,annotation,by="Gene",all.x = T)[order(-SUM)]

# Count the frequency of each feature
output.freq=d3[,.N,by=Gene][order(-N)]

table(output.freq$N==18) #609/610
table(output.freq$N>10)  #610/610

output3=merge(output2,output.freq,by="Gene",all.x = T)[order(-SUM)]

write.table(output3,quote = F,sep = "\t",row.names = F, col.names = T,
            file="./NutriNet-At610logCPM-predict-DW.XGBoost-importantgene-frequency-r50.tsv")

# 4. Session information --------------------------------------------------------------------------
sessionInfo()

#R version 4.0.3 (2020-10-10)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)
#
#Matrix products: default
#
#locale:
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
#[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.1252    
#
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     
#
#other attached packages:
# [1] forcats_0.5.0     stringr_1.4.0     dplyr_1.0.2       purrr_0.3.4       readr_1.4.0      
# [6] tidyr_1.1.2       tibble_3.0.4      tidyverse_1.3.0   ggpubr_0.4.0      ggplot2_3.3.2    
#[11] mlr_2.18.0        ParamHelpers_1.14 data.table_1.13.2 xgboost_1.2.0.1  
#
#loaded via a namespace (and not attached):
# [1] Rcpp_1.0.5        lubridate_1.7.9   lattice_0.20-41   assertthat_0.2.1  digest_0.6.27    
# [6] R6_2.5.0          cellranger_1.1.0  backports_1.1.10  reprex_0.3.0      httr_1.4.2       
#[11] pillar_1.4.6      rlang_0.4.8       curl_4.3          readxl_1.3.1      rstudioapi_0.11  
#[16] car_3.0-10        Matrix_1.2-18     checkmate_2.0.0   labeling_0.4.2    splines_4.0.3    
#[21] foreign_0.8-80    munsell_0.5.0     broom_0.7.2       compiler_4.0.3    modelr_0.1.8     
#[26] pkgconfig_2.0.3   BBmisc_1.11       tidyselect_1.1.0  rio_0.5.16        fansi_0.4.1      
#[31] crayon_1.3.4      dbplyr_2.0.0      withr_2.3.0       grid_4.0.3        jsonlite_1.7.1   
#[36] gtable_0.3.0      lifecycle_0.2.0   DBI_1.1.0         magrittr_1.5      scales_1.1.1     
#[41] zip_2.1.1         cli_2.1.0         stringi_1.5.3     carData_3.0-4     farver_2.0.3     
#[46] ggsignif_0.6.0    fs_1.5.0          parallelMap_1.5.0 xml2_1.3.2        ellipsis_0.3.1   
#[51] generics_0.1.0    vctrs_0.3.4       openxlsx_4.2.3    fastmatch_1.1-0   tools_4.0.3      
#[56] glue_1.4.2        hms_0.5.3         abind_1.4-5       parallel_4.0.3    survival_3.2-7   
#[61] colorspace_1.4-1  rstatix_0.6.0     rvest_0.3.6       haven_2.3.1     

