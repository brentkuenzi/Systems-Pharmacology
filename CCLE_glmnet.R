# Load libraries
library(mlbench)
library(caret)
library(caretEnsemble)
library(glmnet)
library(glmnetcr)
library(plyr)
library(heatmap3)
library(rgl)
library(ggrepel)
setwd("~/Documents/Drug SVM/R_analysis")
# read in CCLE
data <- read.table("CCLE.txt",sep="\t",header=TRUE); data$Cell.Line <- as.character(data$Cell.Line)
data <- unique(data)
rownames(data) <- data$Cell.Line
# define functions
S2N <- function(CCLE,sensitive, resistant){
  data_s <- subset(CCLE, Cell.Line %in% sensitive)
  data_s <- unique(data_s)
  data_r <- subset(CCLE, Cell.Line %in% resistant)
  data_r <- unique(data_r)
  sens_genes_mean <- c(); sens_genes_sd <- c()
  resi_genes_mean <- c(); resi_genes_sd <- c()
  genes_p <- c()
  for(i in 2:length(colnames(data_s))){
    sens_genes_mean <- c(sens_genes_mean, mean(data_s[,i]))
    sens_genes_sd <- c(sens_genes_sd, sd(data_s[,i]))
    resi_genes_mean <- c(resi_genes_mean, mean(data_r[,i]))
    resi_genes_sd <- c(resi_genes_sd, sd(data_r[,i]))
    genes_p <- c(genes_p, t.test(data_s[,i], data_r[,i])$p.value)
  }
  S2N <- (sens_genes_mean - resi_genes_mean) / (sens_genes_sd + resi_genes_sd)
  S2N_df <- data.frame(S2N = S2N, pvalue = genes_p)
  S2N_df$genes <- colnames(data)[2:length(colnames(data))]
  return(S2N_df)
}
volcano.plot <- function(x,y,xcutoffs = c(-2,2),ycutoff,cols = c("black","red"),xlab = "",ylab = ""){
  data <- data.frame(x=x,y=y)
  x1 <- subset(data,x >= xcutoffs[2]) ; x1 <- subset(x1, y <= ycutoff)
  x2 <- subset(data, x <= xcutoffs[1]); x2 <- subset(x2, y <= ycutoff)
  plot(data$x,-log10(data$y),xlab=xlab,ylab=ylab,col=cols[1])
  points(x1$x,-log10(x1$y),col=cols[2])
  points(x2$x,-log10(x2$y),col=cols[2])
}
# feature selection
S2N_test <- S2N(CCLE = data, sensitive = c("NCIH1155_LUNG", "NCIH460_LUNG", "NCIH1299_LUNG", "RS411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "KOPN8_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "NCIH2172_LUNG"),
                resistant=c("NCIH650_LUNG", "NCIH23_LUNG", "NCIH661_LUNG", "HCC2935_LUNG", "NCIH1395_LUNG","NCIH82_LUNG","NCIH526_LUNG","NCIH524_LUNG","NCIH146_LUNG","NCIH226_LUNG"))
#####
"You need to find what top 1% or top 0.1% corresponds to in ABS(S2N) for plot"
#####
volcano.plot(x=S2N_test$S2N,y=S2N_test$pvalue,xcutoffs = c(-1.49,1.49),ycutoff = 0.1,xlab="S2N",ylab="-log10(p-value)") # volcano plot

S2N_test$ABS.S2N <- abs(S2N_test$S2N)
S2N_df <- arrange(S2N_test,desc(ABS.S2N))
S2N_df <- S2N_df[1:(0.01*length(S2N_df$genes)),]
S2N_df <- subset(S2N_df, pvalue < 0.1)

# Define all your cell lines
cells <- c("NCIH1155_LUNG", "NCIH460_LUNG", "NCIH1299_LUNG", "RS411_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "KOPN8_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", "NCIH2172_LUNG", "NCIH650_LUNG", "NCIH23_LUNG", "NCIH661_LUNG", "HCC2935_LUNG", "NCIH1395_LUNG","NCIH82_LUNG","NCIH526_LUNG","NCIH524_LUNG","NCIH146_LUNG","NCIH226_LUNG")
data <- unique(data)
row.names(data) <- data$Cell.Line
data1 <- subset(data, Cell.Line %in% cells)
data1 <- unique(data1)
data1 <- data1[,colnames(data1) %in% S2N_df$genes] # may need to pick better cutoffs for feature selection

heatmap3(as.matrix(data1),scale="column")#,labCol = FALSE) # check if selected genes can cluster sensitive vs. resistant

###
"I do Leave-group out cross validation (Monte Carlo cross validation) but LOOCV or repeatedcv should also work.
You just need to look at your lambda / alpha plots to see which CV strategy results in the best model."
###
control <- trainControl(method="adaptive_LGOCV", number=8, repeats=3, savePredictions=TRUE, classProbs=TRUE,p=0.7)
set.seed(0)
alpha.grid <- seq(0,1, length=10)
lambda.grid <- 10^seq(2,-2,length=100)
srchgrd <- expand.grid(.alpha=alpha.grid, .lambda = lambda.grid)
# if your model keeps overfitting you can assign your own alpha
#srchgrd2 <- expand.grid(.alpha=0.33333333333333, .lambda = lambda.grid)
myElastic <- train(x=data1,y=c("s","r","s","r","r","r","s","r","r","s","r","r","r","r","s","s"), method="glmnet",trControl = control,tuneGrid = srchgrd, standardize=FALSE,maxit=1000000)

###
"Evaluate your model with plots below and by looking at the accuracy / alphas / lambdas"

"Once you get your gene predictors, you should look at the Empirical Cumulative Distribution Functions of them
across your two cell populations. You can get a p-value for how different they are by either a kolmogorov smirnov test
or a ranked kolmogorov smirnov test if they are of vastly different n's.

A ranked KS test is distinct from a standard KS test in several ways: First, when
comparing a large distribution to a small distribution in a standard KS test, the NULL
hypothesis is biased towards being rejected, and thus introduce false positives. Second,
a ranked KS test allows for the preferential ranking of sets that are separated from the
background at the tails of the distribution. I'm working on python code to perform a ranked
KS test (not available in any package).
"
###
lbs_fun <- function(fit, ...) { # labels lines in lambda plot
  L <- length(fit$lambda)
  x <- log(fit$lambda[L])
  y <- fit$beta[, L]
  labs <- names(y)
  text(x, y, labels=labs, ...)
}
plot(myElastic,scales = list(x = list(log = 2)),xlim=c(2^-2,2^5))
plot(myElastic$finalModel,xvar="lambda")
lbs_fun(myElastic$finalModel)
coefs <- coef(myElastic$finalModel,myElastic$bestTune$.lambda)
imps <- arrange(subset(varImp(myElastic$finalModel),Overall !=0),Overall)

###
"CCLE Predictions"
###
test <- data[!(data$Cell.Line %in% cells),] # not working
rownames(test) <- test$Cell.Line
#test <- unique(test)
test <- test[,colnames(test) %in% S2N_df$genes]
test$Cell.Line <- NULL
predictions <- predict(myElastic, newdata = test)
predictions_probs <- predict(myElastic, newdata = test,type="prob")
output <- data.frame(Cell.line = rownames(test),Prediction = predictions,Resistant.Prob = predictions_probs$resistant,Sensitive.Prob = predictions_probs$sensitive)
write.table(output,"glmnet_predictions.txt",sep="\t",quote=FALSE,row.names=FALSE)
