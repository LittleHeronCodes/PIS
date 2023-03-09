###=========================================================###
###  Mouse Lung Fibrosis dataset PIS run _(package version  ###
###=========================================================###

library(yaml)
library(data.table)
library(tidyverse)

DIR_RESOURCE='/spstorage/USERS/yeogha/resources'
library(Lazy2, lib.loc=DIR_RESOURCE)

setDTthreads(15)
theme_set(theme_bw())


out_dir <- 'output/fibrosis_2022-08-23'

library(PIS, lib.loc=DIR_RESOURCE) # 0.2.0.9000


##============================================================
## Mouse gene info

load(paste0(DIR_RESOURCE, '/other/mousegeneInfo.RData' ))

## Ortholog mapping
mmu2hsa <- fread(paste0(DIR_RESOURCE,'/other/mgi2human_genes.csv'))
mmu2hsa$HUMAN_ENTREZ = as.character(mmu2hsa$HUMAN_ENTREZ)
mmu2hsa$MGI_ENTREZ = as.character(mmu2hsa$MGI_ENTREZ)

##=====================================================================
## Bleomycin models
##=====================================================================
data_dir <- '/spstorage/USERS/yeogha/Projects/LungFib/'

## Mouse Bleomycin Fibrosis Models
yml_mmu <- read_yaml(paste0(data_dir,'datafiles_mouse.yaml'))

load(paste0(data_dir, yml_mmu$resultls_rdata))		#resultsLS

##==================================================================================

kpaths <- PIS::readGMT('/spstorage/DB/Kpath/K-Path_genesets.gmt')
kpinfo <- fread('/spstorage/DB/Kpath/K-Path_Reference.txt')
colnames(kpinfo)[1] <- 'Path_ID'

# c2cp <- Lazy2::readGMT('/spstorage/DB/MSigDB/v7.2/c2.cp.v7.4.entrez.gmt')

intgpath <- kpaths

intgpath <- lapply(intgpath, function(hgs) with(mmu2hsa, unique(MGI_ENTREZ[which(HUMAN_ENTREZ %in% hgs)])) )
intgpath <- lapply(intgpath, function(hgs) hgs[which(!is.na(hgs))] )
intgpath <- intgpath[which(sapply(intgpath, length) > 10)]

gspace.path <- Reduce(union, intgpath)
# 13969

##==================================================================================
## PIS parameters

ef_cut <- 2				# EF cut off for pathways to consider
min.overlap <- 0
ncore <- 50

fcos <- seq(1.2,3,0.1)
qcos <- seq(0.01,0.2, 0.01)


## RUN PIS

aid = 'GSE97546A'
resultDF <- resultsLS[[aid]]
resultDF3 <- resultDF %>% 
	# dplyr::rename( adj.P.Val='[[adjusted p value]]', logFC='[[log2 fold change]]', entGene='[[entrez gene id]]') %>% 
	filter(!is.na(entGene)) %>%
	mutate(entGene = as.character(entGene))

gspace <- intersect(unique(as.character(resultDF3$entGene)), gspace.path)

pathways <- lapply(intgpath, function(g) intersect(g,gspace))
pathways <- pathways[which(sapply(pathways, length) > 10)]


geneList.conv <- getGenesByCutoffs(resultDF3, fcos, qcos)

cat('==============\nRunning PIS...\n==============\n')

scoresMat.st1 <- calculatePathwayScores(geneList.conv$up, gspace, pathways, ef_cut, min.overlap, ncore=ncore)
pisres.up <- getPeakResults2(geneList.conv$up, scoresMat.st1,verbose=TRUE)

scoresMat.st2 <- calculatePathwayScores(geneList.conv$dn, gspace, pathways, ef_cut, min.overlap, ncore=ncore)
pisres.dn <- getPeakResults2(geneList.conv$dn, scoresMat.st2,verbose=TRUE)

pisres.up
pisres.dn

