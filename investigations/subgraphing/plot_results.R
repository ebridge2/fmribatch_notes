require(ggplot2)
require(gridExtra)
require(Rmisc)
require(latex2exp)


g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

prepare_case1 <- function(result) {
  data <- data.frame(nedges=c(), dsets=c(), lhat=c(), dset_lhat=c(), tobs=c(), which=c(), p=c())
  for (i in 1:length(result)) {
    res_edge <- result[[i]]
    for (j in 1:(length(res_edge)-1)) {
      res_dsets <- res_edge[[j]]
      ids <- sprintf('%s', res_dsets$dataset)
      this_dat1 <- data.frame(nedges=res_dsets$nedge, dataset=ids, lhat=res_dsets$lhat1, tobs=res_dsets$tstat.alt,
                              test="ses-1", p=c(sum(res_dsets$tstat.null < res_dsets$tstat.alt)/length(res_dsets$tstat.null)))
      this_dat2 <- data.frame(nedges=res_dsets$nedge, dataset=ids, lhat=res_dsets$lhat2, tobs=res_dsets$tstat.alt,
                              test="ses-2", p=c(sum(res_dsets$tstat.null < res_dsets$tstat.alt)/length(res_dsets$tstat.null)))
      data <- rbind(data, rbind(this_dat1, this_dat2))
    }
  }
  data$nedges <- factor(data$nedges)
  data$test <- factor(data$test)
  data$dataset <- factor(data$dataset)
  return(data)
}

prepare_case2 <- function(result) {
  data <- data.frame(nedges=c(), dsets=c(), lhat=c(), dset_lhat=c(), tobs=c(), which=c(), p=c())
  for (i in 1:length(result)) {
    res_edge <- result[[i]]
    for (j in 1:(length(res_edge)-1)) {
      res_dsets <- res_edge[[j]]
      ids <- sprintf('%s', res_dsets$dataset)
      this_dat1 <- data.frame(nedges=res_dsets$nedge, dataset=ids, lhat=res_dsets$lhat1, tobs=res_dsets$tstat.alt,
                              test="half-1", p=c(sum(res_dsets$tstat.null < res_dsets$tstat.alt)/length(res_dsets$tstat.null)))
      this_dat2 <- data.frame(nedges=res_dsets$nedge, dataset=ids, lhat=res_dsets$lhat2, tobs=res_dsets$tstat.alt,
                              test="half-2", p=c(sum(res_dsets$tstat.null < res_dsets$tstat.alt)/length(res_dsets$tstat.null)))
      data <- rbind(data, rbind(this_dat1, this_dat2))
    }
  }
  return(data)
}

prepare_case3 <- function(result) {
  data <- data.frame(nedges=c(), dsets=c(), lhat=c(), dset_lhat=c(), tobs=c(), which=c(), p=c())
  for (i in 1:length(result)) {
    res_edge <- result[[i]]
    for (j in 1:(length(res_edge)-1)) {
      res_dsets <- res_edge[[j]]
      ids <- sprintf('%s', res_dsets$site)
      this_dat1 <- data.frame(nedges=res_dsets$nedge, dataset=ids, lhat=res_dsets$lhat1, tobs=res_dsets$tstat.alt,
                              test=res_dsets$datasets[1], p=c(sum(res_dsets$tstat.null < res_dsets$tstat.alt)/length(res_dsets$tstat.null)))
      this_dat2 <- data.frame(nedges=res_dsets$nedge, dataset=ids, lhat=res_dsets$lhat2, tobs=res_dsets$tstat.alt,
                              test=res_dsets$datasets[2], p=c(sum(res_dsets$tstat.null < res_dsets$tstat.alt)/length(res_dsets$tstat.null)))
      data <- rbind(data, rbind(this_dat1, this_dat2))
    }
  }
  return(data)
}

dset_labs <- c('BNU1', 'BNU3', 'HNU1', 'SWU4')
prepare_case45 <- function(result) {
  data <- data.frame(nedges=c(), dsets=c(), lhat=c(), dset_lhat=c(), tobs=c(), which=c(), p=c())
  for (i in 1:length(result)) {
    res_edge <- result[[i]]
    for (j in 1:(length(res_edge)-1)) {
      res_dsets <- res_edge[[j]]
      ids <- sprintf('%s-%s', res_dsets$datasets[1], res_dsets$datasets[2])
      this_dat1 <- data.frame(nedges=res_dsets$nedge, pair=ids, lhat=res_dsets$lhat1, tobs=res_dsets$tstat.alt,
                              test=res_dsets$datasets[1], p=c(sum(res_dsets$tstat.null < res_dsets$tstat.alt)/length(res_dsets$tstat.null)))
      this_dat2 <- data.frame(nedges=res_dsets$nedge, pair=ids, lhat=res_dsets$lhat2, tobs=res_dsets$tstat.alt,
                              test=res_dsets$datasets[2], p=c(sum(res_dsets$tstat.null < res_dsets$tstat.alt)/length(res_dsets$tstat.null)))
      data <- rbind(data, rbind(this_dat1, this_dat2))
    }
  }
  return(data)
}

plot_case1 <- function(data) {
  data$nedges <- factor(data$nedges)
  p1 <- ggplot(data=data, aes(x=nedges, y=lhat, color=dataset, shape=test)) +
    geom_point(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$\\hat{L}$")) +
    ylim(0, 1) +
    ggtitle("Training Data on half of sessions and performance on other half") +
    theme_bw()
  p2 <- ggplot(data=data, aes(x=nedges, y=tobs, group=dataset, color=dataset)) +
    geom_line(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$\\tau_{observed}$")) +
    ggtitle("Observed Test Statistic") +
    ylim(0, 1) +
    theme_bw()
  p3 <- ggplot(data=data, aes(x=nedges, y=p, group=dataset, color=dataset)) +
    geom_line(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$p$-value")) +
    ggtitle(TeX("$p$-Value of Observed Test Statistic")) +
    ylim(0, 1) +
    theme_bw()
  leg <- g_legend(p1)
  grid.arrange(arrangeGrob(p1 + theme(legend.position=NaN), p2 + theme(legend.position=NaN), p3 + theme(legend.position=NaN), nrow=3),
               leg, nrow=1, widths=c(.88, .12))
  
}

plot_case2 <- function(data) {
  data$nedges <- factor(data$nedges)
  p1 <- ggplot(data=data, aes(x=nedges, y=lhat, color=dataset, shape=test)) +
    geom_point(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$\\hat{L}$")) +
    ylim(0, 1) +
    ggtitle("Training Data on First DSet and performance on second") +
    theme_bw()
  p2 <- ggplot(data=data, aes(x=nedges, y=tobs, group=dataset, color=dataset)) +
    geom_line(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$\\tau_{observed}$")) +
    ggtitle("Observed Test Statistic") +
    ylim(0, 1) +
    theme_bw()
  p3 <- ggplot(data=data, aes(x=nedges, y=p, group=dataset, color=dataset)) +
    geom_line(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$p$-value")) +
    ggtitle(TeX("$p$-Value of Observed Test Statistic")) +
    ylim(0, 1) +
    theme_bw()
  leg <- g_legend(p1)
  grid.arrange(arrangeGrob(p1 + theme(legend.position=NaN), p2 + theme(legend.position=NaN), p3 + theme(legend.position=NaN), nrow=3),
               leg, nrow=1, widths=c(.88, .12))
  
}

plot_case3 <- function(data) {
  data$nedges <- factor(data$nedges)
  p1 <- ggplot(data=data, aes(x=nedges, y=lhat, color=dataset, shape=test)) +
    geom_point(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$\\hat{L}$")) +
    ylim(0, 1) +
    ggtitle("Training Data on First DSet and performance on second") +
    theme_bw()
  p2 <- ggplot(data=data, aes(x=nedges, y=tobs, group=dataset, color=dataset)) +
    geom_line(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$\\tau_{observed}$")) +
    ggtitle("Observed Test Statistic") +
    ylim(0, 1) +
    theme_bw()
  p3 <- ggplot(data=data, aes(x=nedges, y=p, group=dataset, color=dataset)) +
    geom_line(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$p$-value")) +
    ggtitle(TeX("$p$-Value of Observed Test Statistic")) +
    ylim(0, 1) +
    theme_bw()
  leg <- g_legend(p1)
  grid.arrange(arrangeGrob(p1 + theme(legend.position=NaN), p2 + theme(legend.position=NaN), p3 + theme(legend.position=NaN), nrow=3),
               leg, nrow=1, widths=c(.88, .12))
  
}

plot_case45 <- function(data) {
  data$nedges <- factor(data$nedges)
  p1 <- ggplot(data=data, aes(x=nedges, y=lhat, color=pair, shape=test)) +
    geom_point(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$\\hat{L}$")) +
    ylim(0, 1) +
    ggtitle("Training Data on First DSet and performance on second") +
    theme_bw()
  p2 <- ggplot(data=data, aes(x=nedges, y=tobs, group=pair, color=pair)) +
    geom_line(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$\\tau_{observed}$")) +
    ggtitle("Observed Test Statistic") +
    ylim(0, 1) +
    theme_bw()
  p3 <- ggplot(data=data, aes(x=nedges, y=p, group=pair, color=pair)) +
    geom_line(size=2) +
    xlab("Number of Edges in Subgraph") +
    ylab(TeX("$p$-value")) +
    ggtitle(TeX("$p$-Value of Observed Test Statistic")) +
    ylim(0, 1) +
    theme_bw()
  leg <- g_legend(p1)
  grid.arrange(arrangeGrob(p1 + theme(legend.position=NaN), p2 + theme(legend.position=NaN), p3 + theme(legend.position=NaN), nrow=3),
               leg, nrow=1, widths=c(.8, .2))
  
}