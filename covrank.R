# A composite zero-dose vulnerability index for equity assessment
# and spatial prioritization in low and middle-income countries


# Description:
# This script contains the function covrank() used to calculate the 
# regression-based weighting method vulnerability index. 


covrank <- function(Data) {

  # Remove clusters where TotChild is zero

  zero.clust <- which(Data$TotChild == 0)
  if (length(zero.clust)>0) Data <- Data[-zero.clust,]
  
  # Remove TotChild for model
  
  Data.1 <- Data[,-1]
  covlist <- names(Data.1)
  covlist <- covlist[-grep("TotVax", covlist)]

  # Single covariate models
  
  AICs = dev = nulldev = n.par = r2 = pr2 = iter = numeric()
  model = response=character()
  n.iter = 5
  propsub = 0.8
  resp = "cbind(TotVax, TotChild-TotVax)"
  resp1 = "TotVax"
  weights = "TotChild"

  for(i in 1:n.iter) {
    
    subind = sample(1:nrow(Data), propsub*nrow(Data), replace = F)
    subdat = Data[subind,]
    preddat = Data[-subind,]
    
    # Null model for comparison
    
    nullmod = glm(formula(paste(resp, "~1")), data = Data, family = binomial(logit))
  
    for(j in covlist){
    
      form = formula(paste(resp, "~", j))
      pmod = glm(form, data = subdat, family = binomial(logit))
      fullmod = glm(form, data = Data, family = binomial(logit))
      
      # pr2 checks the predictive power of the model 
      # against a 'new' subset of the data
      
      pr2 = c(pr2,cor(pmod$fitted, subdat[[resp1]]/subdat[[weights]])^2) 
      
      # AIC for the model on the full data
      
      AICs = c(AICs, AIC(fullmod))
      dev = c(dev, deviance(fullmod))
      n.par = c(n.par,1)
      iter = c(iter,i)
      r2 = c(r2,cor(fullmod$fitted, Data[[resp1]]/Data[[weights]])^2) 
      model = c(model,j)
      response = c(response, resp1)
      
    }
  }
  
  op = data.frame(model, response, dev, AICs, r2, pr2, n.par)
  singles = plyr::ddply(op, .(model), summarise, model = model[1],
                        response = response[1], dev = dev[1], 
                        pr2 = mean(pr2), AICs = AICs[1], r2 = r2[1])

  # Add deviance reduction
  
  singles$devred = singles$dev/deviance(nullmod)

  # Order by pr2
  
  singles = singles[order(singles$pr2, decreasing = T),]
  singles$rank <- nrow(singles):1
  return(singles)

}
