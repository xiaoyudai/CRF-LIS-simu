# HMM vs CRF, two sample Poisson test


prior.prob <- c(0.8, 0.2)
trainsition_mat = matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow=T)
adj <- matrix(0, n_test, n_test)
for (i in 1:(n_test-1))
{
  adj[i, i+1] <- 1
}
mc <- make.crf(adj, n.states=2)
mc$node.pot[1,] <- prior.prob
for (i in 1:mc$n.edges)
{
  mc$edge.pot[[i]] <- trainsition_mat
}

genData <- function(n_test, mc, lambda1, lambda2){
  s <- as.vector(sample.chain(mc, 1))-1 # hidden status
  o1 <- rep(0,n_test)  
  o2 <- rep(0,n_test)
  id <- which(s==0)
  o1 <- rpois(n_test, lambda1)
  o2[id] <- rpois(length(id), lambda1)
  o2[-id] <- rpois(n_test-length(id), lambda2)
  p <- rep(0,n_test)
  for(i in 1:n_test){
    p[i] <- poisson.test(c(o1[i],o2[i]), c(1,1), alternative = "two.sided")$p.value
  }
  return(list(s=s,o1=o1,o2=o2,p=p))
}




simu_pois2 <- function(mc, lambda1, lambda2, n_test, n_node, size_train, j){
  # alter_comp = 2, pc <- (p_1,...,p_k), f1 <- matrix(c(mu_1,sd_1,mu_2,sd_2,...,mu_k,sd_k), col=2)
  
  # f_0: N(0.5, 1)
  
  # f_1: 
  
  # node feature: log(dnorm(x,0,1))
  
  n_state = 2
  
  # generate training data
  
  s_train = as.vector(NULL, mode = "numeric") # hidden status
  o1_train = as.vector(NULL, mode = "numeric") # observations
  o2_train = as.vector(NULL, mode = "numeric")
  
  for(i in 1:size_train){
    dt = genData(n_test, mc, lambda1, lambda2)  
    s_train = c(s_train, dt$s)
    o1_train = c(o1_train, dt$o1)
    o2_train = c(o2_train, dt$o2)
  }
  
  s_train = matrix(s_train, ncol = n_node, byrow = T) + 1 # status: 1 or 2
  o1_train = matrix(o1_train, ncol = n_node, byrow = T) 
  o2_train = matrix(o2_train, ncol = n_node, byrow = T)
  
  # Training CRF
  adj_train <- matrix(0, n_node, n_node)
  for (i in 1:(n_node-1))
  {
    adj_train[i, i+1] <- 1
  }
  
  crf_train <- make.crf(adj_train, n_state)
  crf_train <- make.features(crf_train, 3, 1)
  crf_train <- make.par(crf_train, 6)
  crf_train$node.par[1,1,1] <- 1
  for (i in 1:crf_train$n.edges)
  {
    crf_train$edge.par[[i]][1,1,] <- 2
    crf_train$edge.par[[i]][1,2,] <- 3
    crf_train$edge.par[[i]][2,1,] <- 4
  }
  crf_train$node.par[,1,2] <- 5
  crf_train$node.par[,1,3] <- 6
  
  
  train_nf <- lapply(1:dim(o1_train)[1], function(i) matrix(1, crf_train$n.nf, crf_train$n.nodes))
  for (i in 1:dim(o1_train)[1])
  {
    train_nf[[i]][2, ] <- o1_train[i,]
    train_nf[[i]][3, ] <- o2_train[i,]
  }
  train_nf[[i]][2, which(train_nf[[i]][2, ]<(-7))] = -7 # avoid over-flow (too large number )
  
  train_ef <- lapply(1:dim(o1_train)[1], function(i) matrix(1, crf_train$n.ef, crf_train$n.edges))
  crf_train <- train.crf(crf_train, s_train, train_nf, train_ef)
  
  
  # convert crf_train to an CRF with the same coefficients and n_node = n_tests
  adj_new <- matrix(0, n_test, n_test)
  for (i in 1:(n_test-1))
  {
    adj_new[i, i+1] <- 1
  }
  crf_new <- make.crf(adj_new, n_state)
  crf_new <- make.features(crf_new, 3, 1)
  crf_new <- make.par(crf_new, 6)
  crf_new$node.par[1,1,1] <- 1
  for (i in 1:crf_new$n.edges)
  {
    crf_new$edge.par[[i]][1,1,] <- 2
    crf_new$edge.par[[i]][1,2,] <- 3
    crf_new$edge.par[[i]][2,1,] <- 4
  }
  crf_new$node.par[,1,2] <- 5
  crf_new$node.par[,1,3] <- 6
  
  crf_new$par = crf_train$par # set the new CRF with trained coefficients
  
  
  # Testing
  
  dt = genData(n_test, mc, lambda1, lambda2)  
  s_test = dt$s
  o1_test = dt$o1
  o2_test = dt$o2
  p_test = dt$p
  p_test[which(p_test==1)] = 0.99999 # avoid z_test = Inf
  z_test = qnorm(p_test)

  em.res<-em1.hmm(z_test, maxiter=500)
  
  lsi.hmm<-em.res$lf # HMM procedure
  
  
  o1_test = matrix(o1_test, nrow = 1, ncol = n_test, byrow = T)
  o2_test = matrix(o2_test, nrow = 1, ncol = n_test, byrow = T)
  
  test.nf <- lapply(1:dim(o1_test)[1], function(i) matrix(1, crf_new$n.nf, crf_new$n.nodes))
  for (i in 1:dim(o1_test)[1])
  {
    test.nf[[i]][2, ] <- o1_test[i,]
    test.nf[[i]][3, ] <- o2_test[i,]
  }
  test.ef <- lapply(1:dim(o1_test)[1], function(i) matrix(1, crf_new$n.ef, crf_new$n.edges))
  hmm.infer <- matrix(0, nrow=dim(o1_test)[1], ncol=dim(o1_test)[2])
  for (i in 1:dim(o1_test)[1])
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








