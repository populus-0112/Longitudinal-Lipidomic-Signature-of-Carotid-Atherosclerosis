#https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html#44_Perform_DESeq2_normalization_and_transformation
# igraph package
#https://rstudio-pubs-static.s3.amazonaws.com/90587_39301ea4b9514d20880b4c9cda81ca4f.html
#https://rstudio-pubs-static.s3.amazonaws.com/337696_c6b008e0766e46bebf1401bea67f7b10.html
# network visualization package: ggnet2
#https://briatte.github.io/ggnet/
#https://mran.microsoft.com/snapshot/2017-12-11/web/packages/ggCompNet/vignettes/examples-from-paper.html
library(WGCNA)
require(flashClust)
library(magrittr)

# read data
plaque=read.table('plaque_progression.csv',sep=',',header=F)  
plaque=plaque[!is.na(plaque[,'S4BMI']),]
plaque=plaque[!is.na(plaque[,'S4LDL']),]
plaque=plaque[!is.na(plaque[,'S4GFR_CKD_EPI']),] 
plaque=plaque[plaque[,'S4HTN2']!='',]

plaque[,'progression']=plaque[,'AFFSEG5']-plaque[,'EFFSEG4']
table(plaque[,'progression'],useNA='always')

plaque[,'progression']=ifelse(plaque[,'progression']>0,1,0)
plaque=plaque[!is.na(plaque[,'progression']),] # delete subjects with missing progression
table(plaque[,'progression'],useNA='always')

nm='AC(10:0)'
if (nm%in%colnames(plaque)){
  start_lipid=which(colnames(plaque)=='AC(10:0)')
} else {
  start_lipid=which(colnames(plaque)=='AC.10.0.')
}

lipids_select=plaque[plaque[,'progression']==1,start_lipid:(start_lipid+1541)] 
dim(lipids_select)

# select model parameters
sft <- pickSoftThreshold(lipids_select,powerVector =seq(1,30,1),
                         dataIsExpr = TRUE,
                         corFnc = bicor,
                         networkType = "signed"
)
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)


TOM=TOMsimilarityFromExpr(lipids_select, corType='bicor',maxPOutliers=0.05, networkType='signed', power = 10)
diss=1-TOM
Tree= hclust(as.dist(diss),method="average")
diag(diss) = NA
#cutreeDynamic(dendro = Tree, distM = as.matrix(1-TOM),deepSplit = 4, pamRespectsDendro = FALSE,minClusterSize = 50)
tree = cutreeHybrid(dendro = Tree, pamStage=FALSE,minClusterSize = 50,# cutHeight = 0.9, ## the maximum joining heights for dendrogram
                    deepSplit = 4, distM = as.matrix(1-TOM))

mergedColors = labels2colors(tree$labels)
table(mergedColors)

# Calculate eigengenes
MEList = moduleEigengenes(lipids_select, colors = mergedColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Call an automatic merging function
merge = mergeCloseModules(lipids_select, mergedColors, cutHeight = 0.15, verbose = 3)
# The merged module colors
mergedColors = merge$colors
table(mergedColors)
# Eigengenes of the new merged modules    
mergedMEs = merge$newMEs
plotEigengeneNetworks(mergedMEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))

plotDendroAndColors(Tree, mergedColors, "", # "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)



