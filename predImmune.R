#' @name predImmune
#' @description A R script to determine tumor samples as one out of three immune subtypes using gene expression profile
#' @param expr A properly normalized expression data (e.g., FPKM or TPM)
#' @param prob.cutoff A numeric cutoff for classification; 0.5 by default
#' @return A data frame indicates immune subtype for each sample
TGCT-Immune <- function(expr        = NULL,
                       prob.cutoff = 0.5) {
  
  # load R packages
  library(varSelRF)
  library(dplyr)
  
  # load internal data
  load("subgroup.RData")
  
  # initial check
  if(max(expr) >= 25){
    message("--please make sure a properly normalized expression data has been provided (e.g., FPKM or TPM); count data is not suitable because it does not consider gene length.")
  }
  
  if(max(expr) < 25 | (max(expr) >= 25 & min(expr) < 0)) {
    message("--expression profile seems to have been standardised (z-score or log transformation), no more action will be performed.")
    expr <- expr
  }
  if(max(expr) >= 25 & min(expr) >= 0){
    message("--log2 transformation done for expression data.")
    expr <- log2(expr + 1)
  }
  # data normalization
  data1 <- scale(t(expr[IM_NonIM.rf$selected.vars,colnames(expr)]))
  
  # prediction
  pred1 <- predict(IM_NonIM.rf$rf.model,
                   newdata = subset(data1, select = IM_NonIM.rf$selected.vars),
                   type = "prob")
  
  train.pred1 <- data.frame(Immune1 = pred1[,1])
  train.pred1$Immune1 <- ifelse(train.pred1$Immune1 > prob.cutoff, "Non-Immune","Immune")
  train.pred1$ID <- rownames(train.pred1)
  Immunesub <- rownames(train.pred1[which(train.pred1$Immune == "Immune"),,drop = F])
  
  data2 <- scale(t(expr[Act_Exh.rf$selected.vars,Immunesub]))
  pred2 <- predict(Act_Exh.rf$rf.model,
                   newdata = subset(data2, select = Act_Exh.rf$selected.vars),
                   type = "prob")
  train.pred2 <- data.frame(Immune2 = pred2[,1])
  train.pred2$Immune2 <- ifelse(train.pred2$Immune2 > prob.cutoff, "Exhausted","Activated")
  train.pred2$ID <- rownames(train.pred2)
  
  Group <- left_join(train.pred1,train.pred2,by="ID")
  rownames(Group) <- Group$ID
  Group$Group <- ifelse(Group$Immune1 == "Non-Immune","Non-Immune",
                        ifelse(Group$Immune2 == "Exhausted","Exhausted","Activated"))
  Group <- Group[,c(2,4)]
  
  # return predicted immune group
  return(Group)
}
