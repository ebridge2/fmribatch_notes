require(ggplot2)
require(gridExtra)
require(Rmisc)

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
    for (j in 1:length(res_edge)) {
      res_dsets <- res_edge[[j]]
      ids <- sprintf('%s', res_dsets$dataset)
      this_dat1 <- data.frame(nedges=res_dsets$nedge, dsets=ids, lhat=res_dsets$lhat1, tobs=res_dsets$tstat.alt,
                             which="ses-1", p=c(sum(res_dsets$tstat.alt < res_dsets$tstat.null)/length(res_dsets$tstat.null)))
      this_dat2 <- data.frame(nedges=res_dsets$nedge, dsets=ids, lhat=res_dsets$lhat2, tobs=res_dsets$tstat.alt,
                              which="ses-", p=c(sum(res_dsets$tstat.alt < res_dsets$tstat.null)/length(res_dsets$tstat.null)))
      data <- rbind(data, rbind(this_dat1, this_dat2))
    }
  }
  return(data)
}
plot_case1 <- function(result) {

  for (i in 1:length(result)) {

  }
}
