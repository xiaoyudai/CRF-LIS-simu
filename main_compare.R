library(ggplot2)
library(gridExtra)
library(grid)
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


library(scales) # https://stackoverflow.com/questions/35511951/r-ggplot2-collapse-or-remove-segment-of-y-axis-from-scatter-plot
squish_trans <- function(from, to, factor) {
  trans <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)
  }
  inv <- function(x) {
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    return(x)
  }
  # return the transformation
  return(trans_new("squished", trans, inv))
}


library(CRF)

source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/rdata1.hmm.R")
source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/rdata.hmm.R")
source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/bwfw1.hmm.R")
source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/bwfw.hmm.R")
source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/em1.hmm.R")
source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/em.hmm.R")
source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/mt.hmm.R")
source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/mt.gw.R")
source("~/Documents/Paper3/Rcode/FDR-HMM/HMM-FDR/mt.bh.R")


source('~/Documents/Paper3/Rcode/CRF/simu_OR.R')
source('~/Documents/Paper3/Rcode/CRF/simu_Mis.R')
source('~/Documents/Paper3/Rcode/CRF/simu_pMis.R')
source('~/Documents/Paper3/Rcode/CRF/simu_mrf2.R')
source('~/Documents/Paper3/Rcode/CRF/simu_pois2.R')


##################### Oracle one-component HMM model    ###############################


alter_comp=1
pc = 1
f0 = c(0,1)
f1 = c(1,1)
# pc=c(0.8,0.2)
# f1 = matrix(c(1,1,2,1), 2, 2, byrow = T) 
trainsition_mat = matrix(c(0.98, 0.02, 0.26, 0.74), 2, 2, byrow=T)
n_test = 10000
n_node = 10000
size_train = 1


hmm.fdr005 <- as.vector(NULL,mode='numeric')
hmm.fnr005 <- as.vector(NULL,mode='numeric')
hmm.fdr01 <- as.vector(NULL,mode='numeric')
hmm.fnr01 <- as.vector(NULL,mode='numeric')

crf.fdr005 <- as.vector(NULL,mode='numeric')
crf.fnr005 <- as.vector(NULL,mode='numeric')
crf.fdr01 <- as.vector(NULL,mode='numeric')
crf.fnr01 <- as.vector(NULL,mode='numeric')

n_rep = 16

for(j in 1:n_rep){
  RES = simu_OR(alter_comp, pc, f0, f1, trainsition_mat, n_test, n_node, size_train)
  hmm.fdr005 = c(hmm.fdr005, as.numeric(RES[6,'fdr_hmm']))
  crf.fdr005 = c(crf.fdr005, as.numeric(RES[6,'fdr_crf']))
  hmm.fdr01 = c(hmm.fdr01, as.numeric(RES[11,'fdr_hmm']))
  crf.fdr01 = c(crf.fdr01, as.numeric(RES[11,'fdr_crf']))
  hmm.fnr005 = c(hmm.fnr005, as.numeric(RES[6,'fnr_hmm']))
  crf.fnr005 = c(crf.fnr005, as.numeric(RES[6,'fnr_crf']))
  hmm.fnr01 = c(hmm.fnr01, as.numeric(RES[11,'fnr_hmm']))
  crf.fnr01 = c(crf.fnr01, as.numeric(RES[11,'fnr_crf']))
  print(j)
  # print(c(hmm.fdr005[j],crf.fdr005[j],hmm.fdr01[j],crf.fdr01[j],hmm.fnr005[j],crf.fnr005[j],hmm.fnr01[j],crf.fnr01[j]))
}

# names = c('hmm.fdr005','crf.fdr005','hmm.fdr01','crf.fdr01','hmm.fnr005','crf.fnr005','hmm.fnr01','crf.fnr01')
# values = c(hmm.fdr005,crf.fdr005,hmm.fdr01,crf.fdr01,hmm.fnr005,crf.fnr005,hmm.fnr01,crf.fnr01)
# dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
# dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,5,2,6,3,7,4,8)])
# 
# pdf('/Users/xiaoyudai/Documents/Paper3/Simu/HMM_CRF/OR/p2.pdf',width=7,height=5.4)
# dt.plt <- ggplot(dat1, aes(x=name, y=value)) + 
#   geom_boxplot() +
#   coord_cartesian(ylim=c(0,0.45)) + 
#   xlab(paste(c('n_test=', n_test, ', n_node=', n_node, ', size_train=', size_train, ', n_rep=',n_rep), collapse = "")) + 
#   ylab('')
#   # xlab(expression(paste(gamma,' = 0.05, b = 50'))) 
#   # geom_hline(yintercept=0.1)
# dt.plt
# dev.off()




names = c('HMM','CRF')
values = c(hmm.fdr005,crf.fdr005)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p1 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot() + # geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,0.15)) + 
  geom_hline(yintercept = 0.05,linetype="dotted") +
  xlab('FDR under level 0.05') +
  ylab('') #+
  #theme(axis.text.x=element_blank(),
        #axis.ticks.x=element_blank())
p1

names = c('HMM','CRF')
values = c(hmm.fnr005,crf.fnr005)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p2 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot() + # geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0.05,0.09)) +
  xlab('FNR under level 0.05') +
  ylab('') 
p2

names = c('HMM','CRF')
values = c(hmm.fdr01,crf.fdr01)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p3 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot() + # geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,0.2)) +
  geom_hline(yintercept = 0.1,linetype="dotted") +
  xlab('FDR under level 0.1') +
  ylab('') 
p3

names = c('HMM','CRF')
values = c(hmm.fnr01,crf.fnr01)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p4 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot() + # geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0.05,0.08)) +
  xlab('FNR under level 0.1') +
  ylab('') 
p4

# pdf('/Users/xiaoyudai/Documents/Paper3/Simu/HMM_CRF/OR/p2.pdf',width=9,height=5.4)
pdf('/Users/xiaoyudai/Documents/Paper3/Latex/Bioinformatics/CRFLIS/pic/or.pdf',width=9,height=2.4)
multiplot(p1, p2, p3, p4, cols=4)
dev.off()




##################### mis-specified    ###############################
alter_comp=3
f0 = c(0.5,1)
# pc = 1
# f1 = c(1,1)
pc=c(0.5,0.2,0.3)
f1 = matrix(c(1,1,2,1,3,1), 3, 2, byrow = T) 
trainsition_mat = matrix(c(0.98, 0.02, 0.26, 0.74), 2, 2, byrow=T)
n_test = 10000
n_node = 10000
size_train = 1


hmm.fdr005 <- as.vector(NULL,mode='numeric')
hmm.fnr005 <- as.vector(NULL,mode='numeric')
hmm.fdr01 <- as.vector(NULL,mode='numeric')
hmm.fnr01 <- as.vector(NULL,mode='numeric')

crf.fdr005 <- as.vector(NULL,mode='numeric')
crf.fnr005 <- as.vector(NULL,mode='numeric')
crf.fdr01 <- as.vector(NULL,mode='numeric')
crf.fnr01 <- as.vector(NULL,mode='numeric')

n_rep = 50

for(j in 1:n_rep){
  RES = simu_Mis(alter_comp, pc, f0, f1, trainsition_mat, n_test, n_node, size_train)
  hmm.fdr005 = c(hmm.fdr005, as.numeric(RES[6,'fdr_hmm']))
  crf.fdr005 = c(crf.fdr005, as.numeric(RES[6,'fdr_crf']))
  hmm.fdr01 = c(hmm.fdr01, as.numeric(RES[11,'fdr_hmm']))
  crf.fdr01 = c(crf.fdr01, as.numeric(RES[11,'fdr_crf']))
  hmm.fnr005 = c(hmm.fnr005, as.numeric(RES[6,'fnr_hmm']))
  crf.fnr005 = c(crf.fnr005, as.numeric(RES[6,'fnr_crf']))
  hmm.fnr01 = c(hmm.fnr01, as.numeric(RES[11,'fnr_hmm']))
  crf.fnr01 = c(crf.fnr01, as.numeric(RES[11,'fnr_crf']))
  print(j)
  # print(c(hmm.fdr005[j],crf.fdr005[j],hmm.fdr01[j],crf.fdr01[j],hmm.fnr005[j],crf.fnr005[j],hmm.fnr01[j],crf.fnr01[j]))
}


names = c('HMM','CRF')
values = c(hmm.fdr005,crf.fdr005)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p1 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot() + 
  scale_y_continuous(trans = squish_trans(0.1, 0.7, 100),
                     breaks = c(0.01, 0.05, 0.1, 0.7, 0.8, 0.9)) + 
  geom_hline(yintercept = 0.05,linetype="dotted") +
  xlab('FDR under level 0.05') +
  ylab('') 
p1

names = c('HMM','CRF')
values = c(hmm.fnr005,crf.fnr005)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p2 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot() + # geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,0.07)) +
  xlab('FNR under level 0.05') +
  ylab('') 
p2

names = c('HMM','CRF')
values = c(hmm.fdr01,crf.fdr01)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p3 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot() +
  # coord_cartesian(ylim=c(0,0.9)) + 
  scale_y_continuous(trans = squish_trans(0.2, 0.8, 100),
                     breaks = c(0.05, 0.1, 0.2, 0.8, 0.9)) +
  geom_hline(yintercept = 0.1,linetype="dotted") +
  xlab('FDR under level 0.1') +
  ylab('') 
p3

names = c('HMM','CRF')
values = c(hmm.fnr01,crf.fnr01)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p4 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot() + # geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,0.07)) +
  xlab('FNR under level 0.1') +
  ylab('') 
p4

# pdf('/Users/xiaoyudai/Documents/Paper3/Simu/HMM_CRF/Mis/p4.pdf',width=9,height=5.4)
pdf('/Users/xiaoyudai/Documents/Paper3/Latex/Bioinformatics/CRFLIS/pic/mis.pdf',width=9,height=2.4)
multiplot(p1, p2, p3, p4, cols=4)
dev.off()



##################### partially mis-specified    ###############################
alter_comp=3
f0 = c(0,1)
# pc = 1
# f1 = c(1,1)
pc=c(0.5,0.2,0.3)
f1 = matrix(c(1,1,2,1,3,1), 3, 2, byrow = T) 
trainsition_mat = matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow=T)
n_test = 1000
n_node = 100
size_train = 1


hmm.fdr005 <- as.vector(NULL,mode='numeric')
hmm.fnr005 <- as.vector(NULL,mode='numeric')
hmm.fdr01 <- as.vector(NULL,mode='numeric')
hmm.fnr01 <- as.vector(NULL,mode='numeric')

crf.fdr005 <- as.vector(NULL,mode='numeric')
crf.fnr005 <- as.vector(NULL,mode='numeric')
crf.fdr01 <- as.vector(NULL,mode='numeric')
crf.fnr01 <- as.vector(NULL,mode='numeric')

n_rep = 30

for(j in 1:n_rep){
  RES = simu_Mis(alter_comp, pc, f0, f1, trainsition_mat, n_test, n_node, size_train)
  hmm.fdr005 = c(hmm.fdr005, as.numeric(RES[6,'fdr_hmm']))
  crf.fdr005 = c(crf.fdr005, as.numeric(RES[6,'fdr_crf']))
  hmm.fdr01 = c(hmm.fdr01, as.numeric(RES[11,'fdr_hmm']))
  crf.fdr01 = c(crf.fdr01, as.numeric(RES[11,'fdr_crf']))
  hmm.fnr005 = c(hmm.fnr005, as.numeric(RES[6,'fnr_hmm']))
  crf.fnr005 = c(crf.fnr005, as.numeric(RES[6,'fnr_crf']))
  hmm.fnr01 = c(hmm.fnr01, as.numeric(RES[11,'fnr_hmm']))
  crf.fnr01 = c(crf.fnr01, as.numeric(RES[11,'fnr_crf']))
  print(j)
  print(c(hmm.fdr005[j],crf.fdr005[j],hmm.fdr01[j],crf.fdr01[j],hmm.fnr005[j],crf.fnr005[j],hmm.fnr01[j],crf.fnr01[j]))
}

names = c('hmm.fdr005','crf.fdr005','hmm.fdr01','crf.fdr01','hmm.fnr005','crf.fnr005','hmm.fnr01','crf.fnr01')
values = c(hmm.fdr005,crf.fdr005,hmm.fdr01,crf.fdr01,hmm.fnr005,crf.fnr005,hmm.fnr01,crf.fnr01)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,5,2,6,3,7,4,8)])

pdf('/Users/xiaoyudai/Documents/Paper3/Simu/HMM_CRF/pMis/p6.pdf',width=7,height=5.4)
dt.plt <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot() +
  coord_cartesian(ylim=c(0,0.3)) + 
  xlab(paste(c('n_test=', n_test, ', n_node=', n_node, ', size_train=', size_train, ', n_rep=',n_rep), collapse = "")) + 
  ylab('')
# xlab(expression(paste(gamma,' = 0.05, b = 50'))) 
# geom_hline(yintercept=0.1)
dt.plt
dev.off()





##################### (2-order MC) mis-specified: theta ~ CRF   ###############################
f0 = c(0.5,1)
# pc = 1
# f1 = c(1,1)
pc=c(0.5,0.2,0.3)
f1 = matrix(c(1,1,2,1,3,1), 3, 2, byrow = T) 
trainsition_mat = matrix(c(0.98, 0.02, 0.26, 0.74), 2, 2, byrow=T)
n_test = 10000
n_node = 10000 # n_test/n_node has to be an integer
size_train = 1

prior.prob <- c(0.8, 0.2)
adj <- matrix(0, n_test, n_test)
for (i in 1:(n_node-1))
{
  adj[i, i+1] <- 1
}
for (i in 1:(n_node-2))
{
  adj[i, i+2] <- 1
}

mc <- make.crf(adj, n.states=2)
mc$node.pot[1,] <- prior.prob
for (i in 1:mc$n.edges)
{
  mc$edge.pot[[i]] <- trainsition_mat
}



hmm.fdr005 <- as.vector(NULL,mode='numeric')
hmm.fnr005 <- as.vector(NULL,mode='numeric')
hmm.fdr01 <- as.vector(NULL,mode='numeric')
hmm.fnr01 <- as.vector(NULL,mode='numeric')

crf.fdr005 <- as.vector(NULL,mode='numeric')
crf.fnr005 <- as.vector(NULL,mode='numeric')
crf.fdr01 <- as.vector(NULL,mode='numeric')
crf.fnr01 <- as.vector(NULL,mode='numeric')

n_rep = 50

for(j in 1:n_rep){
  
  RES = simu_mrf2(mc, pc, f0, f1, trainsition_mat, n_test, n_node, size_train, j)
  hmm.fdr005 = c(hmm.fdr005, as.numeric(RES[6,'fdr_hmm']))
  crf.fdr005 = c(crf.fdr005, as.numeric(RES[6,'fdr_crf']))
  hmm.fdr01 = c(hmm.fdr01, as.numeric(RES[11,'fdr_hmm']))
  crf.fdr01 = c(crf.fdr01, as.numeric(RES[11,'fdr_crf']))
  hmm.fnr005 = c(hmm.fnr005, as.numeric(RES[6,'fnr_hmm']))
  crf.fnr005 = c(crf.fnr005, as.numeric(RES[6,'fnr_crf']))
  hmm.fnr01 = c(hmm.fnr01, as.numeric(RES[11,'fnr_hmm']))
  crf.fnr01 = c(crf.fnr01, as.numeric(RES[11,'fnr_crf']))
  print(j)
  print(c(hmm.fdr005[j],crf.fdr005[j],hmm.fdr01[j],crf.fdr01[j],hmm.fnr005[j],crf.fnr005[j],hmm.fnr01[j],crf.fnr01[j]))
}


for(ii in 1:20){
  n_rep = 10
  for(j in 1:n_rep){
    RES = simu_mrf2(mc, pc, f0, f1, trainsition_mat, n_test, n_node, size_train, j)
    hmm.fdr005 = c(hmm.fdr005, as.numeric(RES[6,'fdr_hmm']))
    crf.fdr005 = c(crf.fdr005, as.numeric(RES[6,'fdr_crf']))
    hmm.fdr01 = c(hmm.fdr01, as.numeric(RES[11,'fdr_hmm']))
    crf.fdr01 = c(crf.fdr01, as.numeric(RES[11,'fdr_crf']))
    hmm.fnr005 = c(hmm.fnr005, as.numeric(RES[6,'fnr_hmm']))
    crf.fnr005 = c(crf.fnr005, as.numeric(RES[6,'fnr_crf']))
    hmm.fnr01 = c(hmm.fnr01, as.numeric(RES[11,'fnr_hmm']))
    crf.fnr01 = c(crf.fnr01, as.numeric(RES[11,'fnr_crf']))
    print(j)
    print(c(hmm.fdr005[j],crf.fdr005[j],hmm.fdr01[j],crf.fdr01[j],hmm.fnr005[j],crf.fnr005[j],hmm.fnr01[j],crf.fnr01[j]))
  }
}







names = c('HMM','CRF')
values = c(hmm.fdr005,crf.fdr005)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p1 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,0.9)) + 
  geom_hline(yintercept = 0.05,linetype="dotted") +
  xlab('FDR under level 0.05') +
  ylab('') #+
#theme(axis.text.x=element_blank(),
#axis.ticks.x=element_blank())
p1

names = c('HMM','CRF')
values = c(hmm.fnr005,crf.fnr005)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p2 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,0.9)) +
  xlab('FNR under level 0.05') +
  ylab('') 
p2

names = c('HMM','CRF')
values = c(hmm.fdr01,crf.fdr01)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p3 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,0.9)) +
  geom_hline(yintercept = 0.1,linetype="dotted") +
  xlab('FDR under level 0.1') +
  ylab('') 
p3

names = c('HMM','CRF')
values = c(hmm.fnr01,crf.fnr01)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p4 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,0.9)) +
  xlab('FNR under level 0.1') +
  ylab('') 
p4

pdf('/Users/xiaoyudai/Documents/Paper3/Simu/HMM_CRF/mc2/p8.pdf',width=9,height=5.4)
multiplot(p1, p2, p3, p4, cols=4)
dev.off()





#################### Pois2: theta ~ CRF   ###############################
n_test = 10000
n_node = 10000 # n_test/n_node has to be an integer
size_train = 1



hmm.fdr005 <- as.vector(NULL,mode='numeric')
hmm.fnr005 <- as.vector(NULL,mode='numeric')
hmm.fdr01 <- as.vector(NULL,mode='numeric')
hmm.fnr01 <- as.vector(NULL,mode='numeric')

crf.fdr005 <- as.vector(NULL,mode='numeric')
crf.fnr005 <- as.vector(NULL,mode='numeric')
crf.fdr01 <- as.vector(NULL,mode='numeric')
crf.fnr01 <- as.vector(NULL,mode='numeric')

n_rep = 30

for(j in 1:n_rep){
  
  RES = simu_pois2(mc, lambda1=10, lambda2=15, n_test, n_node, size_train, j)
  hmm.fdr005 = c(hmm.fdr005, as.numeric(RES[6,'fdr_hmm']))
  crf.fdr005 = c(crf.fdr005, as.numeric(RES[6,'fdr_crf']))
  hmm.fdr01 = c(hmm.fdr01, as.numeric(RES[11,'fdr_hmm']))
  crf.fdr01 = c(crf.fdr01, as.numeric(RES[11,'fdr_crf']))
  hmm.fnr005 = c(hmm.fnr005, as.numeric(RES[6,'fnr_hmm']))
  crf.fnr005 = c(crf.fnr005, as.numeric(RES[6,'fnr_crf']))
  hmm.fnr01 = c(hmm.fnr01, as.numeric(RES[11,'fnr_hmm']))
  crf.fnr01 = c(crf.fnr01, as.numeric(RES[11,'fnr_crf']))
  print(j)
  print(c(hmm.fdr005[j],crf.fdr005[j],hmm.fdr01[j],crf.fdr01[j],hmm.fnr005[j],crf.fnr005[j],hmm.fnr01[j],crf.fnr01[j]))
}

names = c('HMM','CRF')
values = c(hmm.fdr005,crf.fdr005)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p1 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,1)) + 
  geom_hline(yintercept = 0.05,linetype="dotted") +
  xlab('FDR under level 0.05') +
  ylab('') #+
#theme(axis.text.x=element_blank(),
#axis.ticks.x=element_blank())
p1

names = c('HMM','CRF')
values = c(rep(1,length(hmm.fnr005)),crf.fnr005)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p2 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,1)) +
  xlab('FNR under level 0.05') +
  ylab('') 
p2

names = c('HMM','CRF')
values = c(hmm.fdr01,crf.fdr01)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p3 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,1)) +
  geom_hline(yintercept = 0.1,linetype="dotted") +
  xlab('FDR under level 0.1') +
  ylab('') 
p3

names = c('HMM','CRF')
values = c(rep(1,length(hmm.fnr01)),crf.fnr01)
dat1 <- data.frame(name = rep(names, each = n_rep), value = values)
dat1$name=factor(dat1$name , levels=levels(dat1$name)[c(1,2)])
p4 <- ggplot(dat1, aes(x=name, y=value)) + 
  geom_boxplot(colour=c("red", "black")) +
  coord_cartesian(ylim=c(0,1)) +
  xlab('FNR under level 0.1') +
  ylab('') 
p4

pdf('/Users/xiaoyudai/Documents/Paper3/Simu/HMM_CRF/pois2/p1.pdf',width=9,height=5.4)
multiplot(p1, p2, p3, p4, cols=4)
dev.off()






