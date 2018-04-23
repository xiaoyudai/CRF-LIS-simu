source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/rdata1.hmm.R")


simu_Mis <- function(alter_comp=1, pc=1, f0, f1, trainsition_mat, n_test, n_node, size_train){
  # alter_comp = 2, pc <- (p_1,...,p_k), f1 <- matrix(c(mu_1,sd_1,mu_2,sd_2,...,mu_k,sd_k), col=2)
  
  # f_0: N(0.5, 1)
  
  # f_1: 
  
  # node feature: log(dnorm(x,0,1))
  
  n_state = 2
  
  nf_f <- function(x){
    # node feature function: f_0/f_1
    return(log(dnorm(x,0,1)))
  }
  
  
  # generate training data
  
  s_train = as.vector(NULL, mode = "numeric") # hidden status
  o_train = as.vector(NULL, mode = "numeric") # observations
  if(alter_comp == 1){
    for(i in 1:size_train){
      dt = rdata.hmm(n_test, pii = c(1, 0), trainsition_mat, f0, 1, f1)  
      s_train = c(s_train, dt$s)
      o_train = c(o_train, dt$o)
    }
  } else {
    for(i in 1:size_train){
      dt = rdata.hmm(n_test, pii = c(1, 0), trainsition_mat, f0, pc, f1)  
      s_train = c(s_train, dt$s)
      o_train = c(o_train, dt$o)
    }
  }
  s_train = matrix(s_train, ncol = n_node, byrow = T) + 1 # status: 1 or 2
  o_train = matrix(o_train, ncol = n_node, byrow = T) 
  
  
  # Training CRF
  adj_train <- matrix(0, n_node, n_node)
  for (i in 1:(n_node-1))
  {
    adj_train[i, i+1] <- 1
  }
  
  crf_train <- make.crf(adj_train, n_state)
  crf_train <- make.features(crf_train, 2, 1)
  crf_train <- make.par(crf_train, 5)
  crf_train$node.par[1,1,1] <- 1
  for (i in 1:crf_train$n.edges)
  {
    crf_train$edge.par[[i]][1,1,] <- 2
    crf_train$edge.par[[i]][1,2,] <- 3
    crf_train$edge.par[[i]][2,1,] <- 4
  }
  crf_train$node.par[,1,2] <- 5
  
  
  train_nf <- lapply(1:dim(o_train)[1], function(i) matrix(1, crf_train$n.nf, crf_train$n.nodes))
  for (i in 1:dim(o_train)[1])
  {
    train_nf[[i]][2, ] <- nf_f(o_train[i,]) # log(f0/f1)
  }
  train_ef <- lapply(1:dim(o_train)[1], function(i) matrix(1, crf_train$n.ef, crf_train$n.edges))
  crf_train <- train.crf(crf_train, s_train, train_nf, train_ef)
  
  
  # convert crf_train to an CRF with the same coefficients and n_node = n_tests
  adj_new <- matrix(0, n_test, n_test)
  for (i in 1:(n_test-1))
  {
    adj_new[i, i+1] <- 1
  }
  crf_new <- make.crf(adj_new, n_state)
  crf_new <- make.features(crf_new, 2, 1)
  crf_new <- make.par(crf_new, 5)
  crf_new$node.par[1,1,1] <- 1
  for (i in 1:crf_new$n.edges)
  {
    crf_new$edge.par[[i]][1,1,] <- 2
    crf_new$edge.par[[i]][1,2,] <- 3
    crf_new$edge.par[[i]][2,1,] <- 4
  }
  crf_new$node.par[,1,2] <- 5
  
  crf_new$par = crf_train$par # set the new CRF with trained coefficients
  
  
  # Testing
  if(alter_comp == 1){
    dt = rdata.hmm(n_test, pii = c(1, 0), trainsition_mat, f0, 1, f1)  
    s_test = dt$s
    o_test = dt$o
  } else {
    dt = rdata.hmm(n_test, pii = c(1, 0), trainsition_mat, f0, pc, f1)  
    s_test = dt$s
    o_test = dt$o
  }
  em.res<-em1.hmm(o_test, maxiter=500) # always use one-comp hmm model to fit
  lsi.hmm<-em.res$lf # HMM procedure
  
  
  o_test = matrix(o_test, nrow = 1, ncol = n_test, byrow = T)
  test.nf <- lapply(1:dim(o_test)[1], function(i) matrix(1, crf_new$n.nf, crf_new$n.nodes))
  for (i in 1:dim(o_test)[1])
  {
    test.nf[[i]][2, ] <- nf_f(o_test[i,]) # log(f0/f1)
  }
  test.ef <- lapply(1:dim(o_test)[1], function(i) matrix(1, crf_new$n.ef, crf_new$n.edges))
  hmm.infer <- matrix(0, nrow=dim(o_test)[1], ncol=dim(o_test)[2])
  for (i in 1:dim(o_test)[1])
  {
    crf_new <- crf.update(crf_new, test.nf[[i]], test.ef[[i]])
    hmm.infer[i,] <- infer.chain(crf_new)$node.bel[,1]
  }
  lsi.crf <- as.vector(t(hmm.infer))
  
  
  # performance comparison
  id1 = which(s_test==0)
  level = seq(0,0.2,by=0.01)
  RES <- cbind(level, fdr_hmm=0, fdr_crf=0, fnr_hmm=0, fnr_crf=0)
  for(i in 1:nrow(RES)){
    q = RES[i,'level']
    # the decision rule
    crf.res<-mt.hmm(lsi.crf, q)
    crf.de<-crf.res$de 
    hmm.res<-mt.hmm(lsi.hmm, q)
    hmm.de<-hmm.res$de
    
    RES[i,'fdr_hmm'] = 1 - sum(hmm.de*s_test)/sum(hmm.de)
    RES[i,'fdr_crf'] = 1 - sum(crf.de*s_test)/sum(crf.de)
    RES[i,'fnr_hmm'] = sum(hmm.de==0 & s_test==1)/(n_test-sum(hmm.de))
    RES[i,'fnr_crf'] = sum(crf.de==0 & s_test==1)/(n_test-sum(crf.de))
  }  
  
  return(RES)
}