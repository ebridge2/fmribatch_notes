source('./givens.R')

# givens
#=========================#
basepath='/mnt/nfs2/MR/all_mr'
#basepath='/data/all_mr'
nrep = 100  # number of replicates for synthetic bootstrap to get null
nedge=10

# dMRI ------------------------------------
# ========================================#
nroi = 70
atlases = c("desikan")
dwi.dsets = c('BNU1', 'BNU3', 'HNU1', 'SWU4', 'KKI2009', 'Templeton114')

graphobj <- fmriu.io.collection.open_graphs(basepath=paste(basepath,'/dwi/edgelists/', sep=""), datasets = dwi.dsets,
                                            atlases = atlases, gname = "graphs", fmt='edgelist', rtype='array')
graphs = graphobj$graphs
gcpy = graphs
datasets = graphobj$dataset
subjects = graphobj$subjects
sessions = graphobj$sessions

sexpath = paste(basepath,'/phenotypic', sep="")
class = parse_class(sexpath, dwi.dsets, subjects)
sexs = class$sex
# dwi.diseases = class$disease
ages = class$age

dwi.graphs = graphs[!is.nan(ages) & !is.nan(sexs),,]
dwi.datasets = datasets[!is.nan(ages) & !is.nan(sexs)]
dwi.subjects = subjects[!is.nan(ages) & !is.nan(sexs)]
dwi.sessions = sessions[!is.nan(ages) & !is.nan(sexs)]
dwi.ages = ages[!is.nan(ages) & !is.nan(sexs)]
dwi.sexs = sexs[!is.nan(ages) & !is.nan(sexs)]

dwi.bin_graphs = apply(dwi.graphs, c(2,3), function(x) thresh_matrix(x, thresh=0))
dwi.bin_graphs <- aperm(dwi.bin_graphs, c(2,3,1))

result <- model.case45(dwi.bin_graphs, dwi.sexs, dwi.datasets, nedges=xscale.log10(1, 3, n=nedge),
                       tstat=tstat.jaccard, nrep=nrep, xval=FALSE)
saveRDS(result, 'case5dwi.rds')

# givens
#=========================#
basepath='/mnt/nfs2/MR/all_mr'
nrep = 100  # number of replicates for synthetic bootstrap to get null
nedge = 10

# fMRI ------------------------------------
# ========================================#
nroi = 70
atlases = c('desikan-2mm')

fmri.dsets = c('BNU1','BNU3','HNU1','SWU4','BNU2','IBATRT',
               'IPCAS1', 'IPCAS2', 'IPCAS8',
               'NYU1', 'SWU1', 'SWU2', 'SWU3', 'UWM', 'XHCUMS')

graphobj <- fmriu.io.collection.open_graphs(basepath=paste(basepath,'/fmri/ranked/edgelists/', sep=""), datasets = fmri.dsets,
                                            atlases = atlases, fmt='edgelist', rtype='array')
graphs = graphobj$graphs
gcpy = graphs
datasets = graphobj$dataset
subjects = graphobj$subjects
sessions = graphobj$sessions

sexpath = paste(basepath,'/phenotypic/', sep="")
class = parse_class(sexpath, fmri.dsets, subjects)
sexs = class$sex
# dwi.diseases = class$disease
ages = class$age

fmri.graphs = graphs[!is.nan(ages) & !is.nan(sexs),,]
fmri.datasets = datasets[!is.nan(ages) & !is.nan(sexs)]
fmri.subjects = subjects[!is.nan(ages) & !is.nan(sexs)]
fmri.sessions = sessions[!is.nan(ages) & !is.nan(sexs)]
fmri.ages = ages[!is.nan(ages) & !is.nan(sexs)]
fmri.sexs = sexs[!is.nan(ages) & !is.nan(sexs)]

fmri.bin_graphs = apply(fmri.graphs, c(2,3), function(x) thresh_matrix(x, thresh=0))
fmri.bin_graphs <- aperm(fmri.bin_graphs, c(2,3,1))

result <- model.case45(fmri.bin_graphs, fmri.sexs, fmri.datasets, nedges=xscale.log10(1, 3, n=nedge),
                       tstat=tstat.jaccard, nrep=nrep, xval=FALSE)
saveRDS(result, 'case5fmri.rds')