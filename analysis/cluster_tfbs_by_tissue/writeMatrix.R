## Get those clusters with more than 1 entry
tab <- read.table("graph.tree.L9.gz")
tab_fact <- as.factor(tab$V2)
tab_fact <- data.frame(summary(tab_fact)[summary(tab_fact)>1])
tab_fact = data.frame(cluster= rownames(tab_fact), number= tab_fact[1])
tab_fact <- tab_fact[-NROW(tab_fact),] ## Remove (other)
indx <- as.integer(as.character(tab_fact[,1]))

## Read in the sort data table.
#mat <- read.table("intersect_unique_sort.tsv.gz")
#save.image("intersect_unique.sort.Rdata")
load("intersect_unique.sort.Rdata")

## Build a matrix of elements in the same cluster.
k <- 50
plot_mat <- c()
plot_key <- c()
for(i in indx) {
  plot_mat <- c(plot_mat, sample(which(tab$V2 == i), k))
  plot_key <- c(plot_key, rep(i,k))
}
plot_mat <- mat[plot_mat,]

## Build a distance matrix and order rows by clustering.
d <- dist(t(plot_mat), method="manhattan")
h <- hclust(d)

## Write out a heatmap
require(lattice)
require(grid)
require(gridExtra)

M3 <- levelplot(as.matrix(plot_mat[,h$order]), cuts=1, labels=FALSE, axes=FALSE, scales=list(draw=FALSE), col.regions=c("white", "red"), 
	xlab=list(cex=4), ylab=list(cex=4), colorkey=list(cex.axis=4),
	panel= 
	function(...) {
		panel.levelplot(...) 
		for(i in seq(k, NROW(indx)*k, by= k)) {
			panel.abline(v=i)
		}
	})

png("~/transfer/matrix.png", width = 12000, height = 500)
  print(M3)
dev.off()


## Plot more of each cluster.
k <- 1000

png("~/transfer/rows-clusters-matrix.png", width = 6*k, height = 500*NROW(indx))

for(i in 1:NROW(indx)) {
  plot_mat <- mat[sample(which(tab$V2 == indx[i]), k),]
  plot_mat_l <- levelplot(as.matrix(plot_mat[,h$order]), cuts=1, labels=FALSE, axes=FALSE, scales=list(draw=FALSE), col.regions=c("white", "red"))
  print(plot_mat_l, split = c(1, i, 1, NROW(indx)), more=(i<NROW(indx)))
}

dev.off()

## Visualize as a UMAP.
require(umap)
require(DescTools)
library(RColorBrewer)

k <- 100
plot_mat <- c()
plot_key <- c()
for(i in indx) {
  plot_mat <- c(plot_mat, sample(which(tab$V2 == i), k))
  plot_key <- c(plot_key, rep(i,k))
}
plot_mat <- mat[plot_mat,]

## Produce a distance matrix with a better distance function.
dist.fun <- function(matr, origin, target) {
# sum(xor(matr[,origin], matr[,target])) / (sum((matr[,origin] & matr[,target])) + sum(xor(matr[,origin], matr[,target])))
 sum((matr[,origin] | matr[,target])) / sum((matr[,origin] & matr[,target]))-1
}
umap.settings <- umap.defaults
umap.settings$metric <- dist.fun #"manhattan"

um <- umap(plot_mat, config= umap.settings)

n <- NROW(indx)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cols <- SetAlpha(cols, alpha=0.5)


png("~/transfer/Tfbs-UMAP.png", width=750, height=750)
 cols <- sample(colors(), NROW(indx))
 plot(y= um$layout[,1], x= um$layout[,2], col=cols[as.integer(as.factor(plot_key))], pch=19, xlab="UMAP2", ylab="UMAP1", cex=1)
dev.off()



