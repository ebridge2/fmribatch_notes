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
  jaccard <- array(NaN, dim=c(d))
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

model.compute <- function(samp, Y, Z, global.ssg, nedge, tstat=tstat, verb=FALSE, xval=FALSE) {
  zset = unique(Z)
  # alternate consists of the tstat computed between the subgraphs estimated per-z
  # with the global sub-graph
  lhat <- array(NaN, dim=c(length(zset)))
  ssgs <- array(0, dim=c(nroi, nroi, length(zset)))  # save the ssgs per dataset
  for (i in 1:length(zset)) {
    test.ss <- Z == zset[i]
    test.gr <- samp[,,test.ss]
    test.Y <- Y[test.ss]
    if (xval) {
      train.ss <- Z != zset[i]
    } else {
      train.ss <- test.ss
    }
    train.gr <- samp[,,train.ss]
    train.Y <- Y[train.ss]
    sg_ests <- sg.bern.subgraph_estimator(train.gr, train.Y)
    # compute test statistics per edge from the contingency matrix
    tstats <- sg.bern.edge_test(sg_ests$cont_matrix)
    sg_edge <- sg.bern.subgraph_edge_estimation(tstats, nedge)
    sg <- array(0, dim=c(nroi, nroi))
    sg[sg_edge] <- 1
    ssgs[,,i] <- sg
    Yhat <- sg.bern.subgraph_classifier(test.gr, sg_edge, sg_ests$p, sg_ests$pi, sg_ests$classes)
    lhat[i] <- (length(test.Y) - sum(Yhat$Yhat == test.Y))/length(test.Y)
  }
  tstat <- do.call(tstat, list(ssgs, global.ssg))
  ret <- list(tstats=tstat, ssgs=ssgs, lhat=lhat)
  return(ret)
}

model.null.permutation <- function(samp, Y, Z, global.ssg, global.est, ne, nrep=100, tstat=tstat.jaccard, verb=FALSE, xval=FALSE) {
  nroi = dim(samp)[2]
  zset = unique(Z)
  tstat.null <- array(NaN, dim=c(nrep))
  null.counts <- array(NaN, dim=c(length(zset) + 1))
  lhat.null <- array(NaN, dim=c(nrep, length(zset)))
  class.count <- list()
  for (k in 1:length(global.est$classes)) {
    class.count[[k]] <- which(Y == global.est$classes[k])  # indices of where each class is
  }
  for (i in 1:nrep) {
    yclass = class.count
    null.ssgs <- array(0, dim=c(nroi, nroi, length(zset)))
    counter <- 1
    dats <- list()  # synthetic null data
    y.labs <- list()  # labels for synthetic null data
    z.labs <- list()
    for (j in 1:length(zset)) {
      ss = Z == zset[j]  # compute which indices correspond to our current study
      subset.Y = Y[ss]  # partition the class abels for the study we want to construct synthetic data under the null for
      ns <- array(0, dim=c(length(global.est$classes)))
      for (k in 1:length(global.est$classes)) {
        ns[k] <- sum(subset.Y == global.est$classes[k])  # compute the number of each class we will sample for this synthetic study
        idx.samp <- sample(yclass[[k]], ns[k], replace=FALSE)  # sample ns[k] graphs from our graphs from all studies
        yclass[[k]] <- setdiff(yclass[[k]], idx.samp)  # remove already-chosen indices from our particular sample
        dats[[counter]] <- samp[,,idx.samp]  # and sample from the original collection of graphs
        y.labs[[counter]] <- Y[idx.samp]  # along w the appropriate labels
        z.labs[[counter]] <- j
        counter <- counter + 1
      }
    }
    perm.samp <- abind(dats, along=3)
    perm.Y <- abind(y.labs, along=1)
    perm.Z <- abind(z.labs, along=1)
    global.est <- sg.bern.subgraph_estimator(perm.samp, perm.Y)
    # compute test statistics per edge from the contingency matrix
    global.tstat <- sg.bern.edge_test(global.est$cont_matrix)
    global.edge <- sg.bern.subgraph_edge_estimation(global.tstat, ne)
    global.ssg <- array(0, dim=c(nroi, nroi))
    result <- model.compute(perm.samp, perm.Y, perm.Z, global.ssg, ne, tstat=tstat, verb=verb, xval=xval)
    tstat.null[i] <- result$tstat
    lhat.null[i,] <- result$lhat
  }
  rlist <- list(tstat.null=tstat.null, lhat.null=lhat.null)
  return(rlist)
}

model.case1 <- function(samp, Y, Z, sessions, nedges=xscale.log10(1, 3), tstat=tstat.jaccard, nrep=nrep, xval=FALSE) {
  nroi <- dim(samp)[2]
  zset = unique(Z)
  result <- list()
  for (j in 1:length(nedges)) {
    print(paste("nedges:", nedges[j], "iters:", j, "of:", length(nedges)))
    ne <- nedges[j]
    global.est <- sg.bern.subgraph_estimator(samp, Y)
    # compute test statistics per edge from the contingency matrix
    global.tstat <- sg.bern.edge_test(global.est$cont_matrix)
    global.edge <- sg.bern.subgraph_edge_estimation(global.tstat, ne)
    global.ssg <- array(0, dim=c(nroi, nroi))
    global.ssg[global.edge] <- 1
    edge_result <- list()
    for (i in 1:length(zset)) {
      print(paste("dataset:", zset[i], "iters:", i, "of:", length(zset)))
      ss <- Z == zset[i]
      ses_unique <- sample(unique(sessions[ss]))  # find the unique labels and resort
      ss1 <- ss & sessions %in% ses_unique[1:(length(ses_unique)/2)]  # first half of sessions randomly
      ss2 <- ss & sessions %in% ses_unique[(length(ses_unique)/2 + 1):length(ses_unique)]  # second half of sessions randomly
      res <- model.analyze(samp[,,ss1], samp[,,ss2], Y[ss1], Y[ss2], sessions[ss1], sessions[ss2], ne,
                           global.ssg, global.est, tstat=tstat, nrep=nrep, xval=xval)
      res$dataset <- zset[i]
      edge_result[[i]] <- res
    }
    edge_result$nedge <- ne
    result[[j]] <- edge_result
  }
  return(result)
}

model.case2 <- function(samp, Y, Z, subjects, nedges=xscale.log10(1, 3), tstat=tstat.jaccard, nrep=nrep, xval=FALSE) {
  nroi <- dim(samp)[2]
  zset = unique(Z)
  result <- list()
  for (j in 1:length(nedges)) {
    print(paste("nedges:", nedges[j], "iters:", j, "of:", length(nedges)))
    ne <- nedges[j]
    global.est <- sg.bern.subgraph_estimator(samp, Y)
    # compute test statistics per edge from the contingency matrix
    global.tstat <- sg.bern.edge_test(global.est$cont_matrix)
    global.edge <- sg.bern.subgraph_edge_estimation(global.tstat, ne)
    global.ssg <- array(0, dim=c(nroi, nroi))
    global.ssg[global.edge] <- 1
    edge_result <- list()
    for (i in 1:length(zset)) {
      print(paste("dataset:", zset[i], "iters:", i, "of:", length(zset)))
      ss <- Z == zset[i]
      sub_unique <- sample(unique(subjects[ss]))  # find the unique labels and resort
      ss1 <- ss & subjects %in% sub_unique[1:(length(sub_unique)/2)]  # first half of sessions randomly
      ss2 <- ss & subjects %in% sub_unique[(length(sub_unique)/2 + 1):length(sub_unique)]  # second half of sessions randomly
      res <- model.analyze(samp[,,ss1], samp[,,ss2], Y[ss1], Y[ss2], subjects[ss1], subjects[ss2], ne,
                           global.ssg, global.est, tstat=tstat, nrep=nrep, xval=xval)
      res$dataset <- zset[i]
      edge_result[[i]] <- res
    }
    edge_result$nedge <- ne
    result[[j]] <- edge_result
  }
  return(result)
}

model.case3 <- function(samp, Y, Z, sites, nedges=xscale.log10(1, 3), tstat=tstat.jaccard, nrep=nrep, xval=FALSE) {
  nroi <- dim(samp)[2]
  zset = unique(Z)
  sset <- unique(sites)
  result <- list()
  for (j in 1:length(nedges)) {
    print(paste("nedges:", nedges[j], "iters:", j, "of:", length(nedges)))
    ne <- nedges[j]
    global.est <- sg.bern.subgraph_estimator(samp, Y)
    # compute test statistics per edge from the contingency matrix
    global.tstat <- sg.bern.edge_test(global.est$cont_matrix)
    global.edge <- sg.bern.subgraph_edge_estimation(global.tstat, ne)
    global.ssg <- array(0, dim=c(nroi, nroi))
    global.ssg[global.edge] <- 1
    edge_result <- list()
    for (i in 1:length(sset)) {
      print(paste("site:", sset[i], "iters:", i, "of:", length(sset)))
      ss <- sites == sset[i]
      zset_unique <- unique(Z[ss])
      ss1 <- ss & (Z == zset_unique[1])  # first half of sessions randomly
      ss2 <- ss & (Z == zset_unique[2])  # second half of sessions randomly
      res <- model.analyze(samp[,,ss1], samp[,,ss2], Y[ss1], Y[ss2], sites[ss1], sites[ss2], ne,
                           global.ssg, global.est, tstat=tstat, nrep=nrep, xval=xval)
      res$site <- sset[i]
      res$datasets <- c(zset_unique[1], zset_unique[2])
      edge_result[[i]] <- res
    }
    edge_result$nedge <- ne
    result[[j]] <- edge_result
  }
  return(result)
}

model.case45 <- function(samp, Y, Z, nedges=xscale.log10(1, 3), tstat=tstat.jaccard, nrep=nrep, xval=FALSE) {
  nroi <- dim(samp)[2]
  zset = unique(Z)
  result <- list()
  for (j in 1:length(nedges)) {
    print(paste("nedges:", nedges[j], "iters:", j, "of:", length(nedges)))
    ne <- nedges[j]
    global.est <- sg.bern.subgraph_estimator(samp, Y)
    # compute test statistics per edge from the contingency matrix
    global.tstat <- sg.bern.edge_test(global.est$cont_matrix)
    global.edge <- sg.bern.subgraph_edge_estimation(global.tstat, ne)
    global.ssg <- array(0, dim=c(nroi, nroi))
    global.ssg[global.edge] <- 1
    edge_result <- list()
    counter <- 1
    for (i in 1:(length(zset)-1)) {
      print(paste("dataset:", zset[i], "iters:", i, "of:", length(zset)))
      ss1 <- Z == zset[i]
      for (k in (i+1):length(zset)) {
        print(paste("indataset:", zset[k], "iters:", k, "of:", length(zset)))
        ss2 <- Z == zset[k]
        res <- model.analyze(samp[,,ss1], samp[,,ss2], Y[ss1], Y[ss2], Z[ss1], Z[ss2], ne,
                             global.ssg, global.est, tstat=tstat, nrep=nrep, xval=xval)
        res$datasets <- c(zset[i], zset[k])
        edge_result[[counter]] <- res
        counter <- counter + 1
      }
    }
    edge_result$nedge <- ne
    result[[j]] <- edge_result
  }
  return(result)
}

model.analyze <- function(samp1, samp2, Y1, Y2, Z1, Z2, nedge, global.ssg, global.est, tstat=tstat.jaccard, nrep=100, xval=FALSE) {
  samp <- abind(samp1, samp2, along=3)
  Y <- c(Y1, Y2)
  Z <- c(array(0, dim=c(length(Z1))), array(1, dim=c(length(Z2))))
  # classifier built on first subset and applied to second
  sg_ests <- sg.bern.subgraph_estimator(samp1, Y1)
  # compute test statistics per edge from the contingency matrix
  tstats <- sg.bern.edge_test(sg_ests$cont_matrix)
  class1 <- sg.bern.subgraph_edge_estimation(tstats, nedge)
  
  Yhat2 <- sg.bern.subgraph_classifier(samp2, class1, sg_ests$p, sg_ests$pi, sg_ests$classes)$Yhat
  lhat2 <- (length(Y2) - sum(Yhat2 == Y2))/length(Y2)
  
  # classifier built on second subset and applied to first
  sg_ests <- sg.bern.subgraph_estimator(samp2, Y2)
  # compute test statistics per edge from the contingency matrix
  tstats <- sg.bern.edge_test(sg_ests$cont_matrix)
  class2 <- sg.bern.subgraph_edge_estimation(tstats, nedge)
  
  Yhat1 <- sg.bern.subgraph_classifier(samp1, class2, sg_ests$p, sg_ests$pi, sg_ests$classes)$Yhat
  lhat1 <- (length(Y1) - sum(Yhat1 == Y1))/length(Y1)
  
  # estimate null and alternative model
  alt <- model.compute(samp, Y, Z, global.ssg, nedge, tstat=tstat, xval=xval)
  null <- model.null.permutation(samp, Y, Z, global.ssg, global.est, nedge, nrep=nrep, tstat=tstat.jaccard, xval=xval)
  
  return(list(tstat.alt=alt$tstat, tstat.null=null$tstat.null, lhat.alt=alt$lhat,
              lhat.null=null$lhat.null, nedge=nedge, lhat1=lhat1, lhat2=lhat2, xval=xval))
}
