# load packages
#==========================#
require(fmriutils)
require(VCA)
require(stringr)
require(lme4)
require(RColorBrewer)
require(oro.nifti)
require(reshape2)
require(latex2exp)
require(ggplot2)
require(Rmisc)
require(subgraphing)
require(abind)

# functions
#==========================#

xscale.log2 <- function(a=3, b=6, n=10) {
  return(round(2^(seq(a, b, length.out=n))))
}

xscale.log10 <- function(a=1, b=3, n=10) {
  return(round(10^(seq(a, b, length.out=n))))
}

# accepts a matrix and thresholds/binarizes it
thresh_matrix = function(matrix, thresh=0.5) {
  thr = quantile(matrix, thresh)
  return(ifelse(matrix > thr, 1, 0))
}

# Sum of Squared Residual Test Statistic
#
# local [n, n, d]: is the sum of all dataset ssgs estimated
# global [n, n]: is the globally estimated ssg using all of the data points
tstat.ssr <- function(local, global) {
  local <- apply(local, c(1,2), sum)
  return(sum((local/max(local) - global)^2))
}

# Jaccard Index Test Statistic
# returns mean jaccard index btwn the ssgs within dataset and global
tstat.jaccard <- function(local, global) {
  d <- dim(local)[3]
  idx.glob <- which(global > 0)
  jaccard <- array(0, dim=c(d))
  for (i in 1:dim(local)[3]) {
    idx.loc <- which(local[,,i] > 0)
    jaccard[i] <- length(intersect(idx.loc, idx.glob))/length(union(idx.loc, idx.glob))
  }
  return(mean(jaccard))
}

parse_class <- function(basepath, dsets, subjects) {
  sex = list()
  # disease = list()
  age = list()
  include = c()
  for (dset in dsets) {
    path_to_file = file.path(basepath, paste(dset, "_phenotypic_data.csv", sep=""))
    tab = read.csv(path_to_file)
    tab$SEX[tab$SEX == '#' | is.na(tab$SEX) | is.nan(tab$SEX)] = NaN
    if (dset == 'KKI2009') {
      sexm <- 'M'
    } else if (dset == 'Templeton114') {
      sexm <- 1
    } else {
      sexm <- 2
    }
    tab$AGE_AT_SCAN_1 <- factor(tab$AGE_AT_SCAN_1)
    tab$SEX = tab$SEX == sexm
    tab = tab[complete.cases(tab$SEX),]
    tab$AGE_AT_SCAN_1 = as.numeric(levels(tab$AGE_AT_SCAN_1))[tab$AGE_AT_SCAN_1]
    for (idx in 1:dim(tab)[1]) {
      subid = toString(tab[idx,]$SUBID)
      sex[[subid]] = tab[idx,]$SEX
      age[[subid]] = tab[idx,]$AGE_AT_SCAN_1
      # disease[[subid]] = tab[idx,]$DSM_IV_TR
    }
  }

  sclass = array(NaN, dim=c(length(subjects)))
  ageclass = sclass
  # diseaseclass = sclass
  for (i in 1:length(subjects)) {
    subject = subjects[i]
    subid = sub('^0+(?=[1-9])', '', str_extract(subject, '(?<=sub-).*'), perl=TRUE)
    idx = which(names(sex) == subid)
    if (length(idx) >= 1) {
      sclass[i] <- sex[[subid]]
      ageclass[i] <- age[[subid]]
      # diseaseclass[i] <- disease[[subid]]
    }
  }
  return(list(sex=sclass, age=ageclass))#, disease=diseaseclass))
}

model.alternative <- function(samp, Y, Z, global.ssg, tstat=tstat, sig=.01, verb=FALSE, classify=FALSE) {
  zset = unique(Z)
  global.sz <- sum(global.ssg > 0)
  # alternate consists of the tstat computed between the subgraphs estimated per-z
  # with the global sub-graph
  zset.sz <- array(NaN, dim=c(length(zset)))
  if (classify) {
    lhat <- array(NaN, dim=c(length(zset)))
  }
  alt.ssgs <- array(0, dim=c(nroi, nroi, length(zset)))  # save the ssgs per dataset
  for (i in 1:length(zset)) {
    ss = Z == zset[i]
    subset.gr = samp[,,ss]
    subset.Y = Y[ss]
    sg_ests <- sg.bern.subgraph_estimator(subset.gr, subset.Y)
    # compute test statistics per edge from the contingency matrix
    tstats <- sg.bern.edge_test(sg_ests$cont_matrix)
    sg_edge <- sg.bern.subgraph_edge_estimation(tstats, global.sz)
    sg <- array(0, dim=c(nroi, nroi))
    sg[sg_edge] <- 1
    alt.ssgs[,,i] <- sg
    if (classify) {
      Yhat <- sg.bern.subgraph_classifier(subset.gr, sg_edge, sg_ests$p, sg_ests$pi, sg_ests$classes)
      lhat[i] <- (sum(ss) - sum(Yhat$Yhat == Y[ss]))/sum(ss)
    }
  }

  tstat.alt <- do.call(tstat, list(alt.ssgs, global.ssg))
  ret <- list(tstat.alt=tstat.alt, alt.ssgs=alt.ssgs)
  if (classify) {
    ret$lhat <- lhat
  }
  return(ret)
}

model.null.permutation <- function(samp, Y, Z, global.ssg, global.est, nrep=100, tstat=tstat, sig=.01, verb=FALSE) {
  nroi = dim(samp)[2]
  zset = unique(Z)
  tstat.null <- array(NaN, dim=c(nrep))
  null.counts <- array(0, dim=c(length(zset) + 1))

  for (i in 1:nrep) {
    null.ssgs <- array(0, dim=c(nroi, nroi, length(zset)))
    yclass <- list()
    for (k in 1:length(global.est$classes)) {
      yclass[[k]] <- which(Y == global.est$classes[k])  # indices of where each class is
    }
    for (j in 1:length(zset)) {
      ss = Z == zset[j]  # compute which indices correspond to our current study
      subset.Y = Y[ss]  # partition the class abels for the study we want to construct synthetic data under the null for
      ns <- array(0, dim=c(length(global.est$classes)))
      dats <- list()  # synthetic null data
      labs <- list()  # labels for synthetic null data
      for (k in 1:length(global.est$classes)) {
        ns[k] <- sum(subset.Y == global.est$classes[k])  # compute the number of each class we will sample for this synthetic study
        idx.samp <- sample(yclass[[k]], ns[k], replace=FALSE)  # sample ns[k] graphs from our graphs from all studies
        yclass[[k]] <- setdiff(yclass[[k]], idx.samp)  # remove already-chosen indices from our particular sample
        dats[[k]] <- samp[,,idx.samp]  # and sample from the original collection of graphs
        labs[[k]] <- Y[idx.samp]  # along w the appropriate labels
      }
      dset.null.dat <- abind(dats, along=3)
      dset.null.labs <- abind(labs, along=1)
      sg_ests <- sg.bern.subgraph_estimator(dset.null.dat, dset.null.labs)
      # compute test statistics per edge from the contingency matrix
      tstats <- sg.bern.edge_test(sg_ests$cont_matrix)
      sg_edge <- sg.bern.subgraph_edge_estimation(tstats, global.sz)
      sg <- array(0, dim=c(nroi, nroi))
      sg[sg_edge] <- 1
      null.ssgs[,,j] = sg
    }
  }
}

subgraph.dset_lhat <- function(samp, Y, Z, global_sg, alt.counts, sg_ests, tstat=tstat) {
  zset = unique(Z)
  un_cts <- sort(unique(c(unique(as.numeric(alt.counts)), 0)))
  niter <- (length(un_cts))*length(zset)
  count = array(0, dim=c(niter))
  dataset = array(NaN, dim=c(niter))
  lhat = array(0, dim=c(niter))
  for (i in un_cts) {
    if (i > 0) {
      sg <- which(alt.counts == i)  # get the appropriate edges
    } else if (i == 0) {
      sg <- which(global_sg > 0)
    }

    for (j in 1:length(zset)) {
      ss = Z == zset[j]
      subset.gr = samp[,,ss]
      subset.Y = Y[ss]
      iter <- (i)*(length(zset) + 1) + j
      Yhat <- sg.bern.subgraph_classifier(subset.gr, sg, sg_ests$p, sg_ests$pi, sg_ests$classes)
      lhat[iter] <- (sum(ss) - sum(Yhat$Yhat == Y[ss]))/sum(ss)
      dataset[iter] <- zset[j]
      count[iter] <- i
      if (i > 0) {
        count[iter] <- i
      } else if (i == 0) {
        count[iter] <- 'global sg'
      }
    }
    iter <- iter + 1
    dataset[iter] <- "global"
    if (i > 0) {
      count[iter] <- i
    } else if (i == 0) {
      count[iter] <- 'global sg'
    }
    Yhat <- sg.bern.subgraph_classifier(samp, sg, sg_ests$p, sg_ests$pi, sg_ests$classes)
    lhat[iter] <- (length(Y) - sum(Yhat$Yhat == Y))/length(Y)
  }
  result <- data.frame(count=count, dataset=dataset, lhat=lhat)
  return(result)
}


# givens
#=========================#
sig = 0.01  # global param for significance threshold
nrep = 100  # number of replicates for synthetic bootstrap to get null

