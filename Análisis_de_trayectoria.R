#-------Presentación 1 scRNA seq--------
#.....Alejandro Alarcón del Carmen......
#---------------------------------------
#Análisis de datos de células hematopoyeticas de ratón
#Nestorowa et al. 2016
#

#Librerias necesarias
library(scRNAseq)
library(AnnotationHub)
#Importe de datos y generación de la matriz de cuentas.
sce.nest<- NestorowaHSCData()
#Anotación de cada gen conforme al genoma de ratón
ens.mm.v97<- AnnotationHub()[["AH73905"]]
anno <- select(ens.mm.v97, keys=rownames(sce.nest), 
               keytype="GENEID", columns=c("SYMBOL", "SEQNAME"))
rowData(sce.nest) <- anno[match(rownames(sce.nest), anno$GENEID),]
#Visualizar los datos del set
sce.nest
#Filtrado de datos
library(scater)
stats <- perCellQCMetrics(sce.nest)
qc<- quickPerCellQC(stats, percent_subsets = "altexps_ERCC_percent")
sce.nest<-sce.nest[,!qc$discard]
#Normalización
library(scran)
set.seed(101000110)
clusters <- quickCluster(sce.nest)
sce.nest <- computeSumFactors(sce.nest, clusters=clusters)
sce.nest <- logNormCounts(sce.nest)
sce.nest
#Modelado de varianza
set.seed(00010101)
dec.nest <- modelGeneVarWithSpikes(sce.nest, "ERCC")
top.nest <- getTopHVGs(dec.nest, prop=0.1)
#Reducción de dimensiones
set.seed(101010011)
sce.nest <- denoisePCA(sce.nest, technical=dec.nest, subset.row=top.nest)
sce.nest <- runTSNE(sce.nest, dimred="PCA")
#Clustering y gráfico TSNE
snn.gr <- buildSNNGraph(sce.nest, use.dimred="PCA")
colLabels(sce.nest) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
plotTSNE(sce.nest, colour_by="label")
#Detección de genes marcadores en los grupos o clusters
markers <- findMarkers(sce.nest, colLabels(sce.nest), 
                       test.type="wilcox", direction="up", lfc=0.5,
                       row.data=rowData(sce.nest)[,"SYMBOL",drop=FALSE])
#En este caso específico se pueden identificar genes relacionados con la
#identiad celular, por ejemplo Car2 y Hebp1 marcan que el cluster 8 es 
#de precursores de eritrocitos
markers$`8`
#Proceso de anotación de los clusters, asignación de identidad
#celular
library(SingleR)
mm.ref<- MouseRNAseqData()
#Cambiamos el número de cluster por el nombre de identidad
#celular asignado
renamed <- sce.nest
rownames(renamed)<- uniquifyFeatureNames(rownames(renamed), 
                                         rowData(sce.nest)$ SYMBOL)
label<- SingleR(renamed, mm.ref, labels = mm.ref$label.fine)
tab <- table(labels$labels, colLabels(sce.nest))
#Inicio de análisis de trayectoria: Dibujando la línea de trayectoria
library(scater)
by.cluster <- aggregateAcrossCells(sce.nest, ids=colLabels(sce.nest))
centroids <- reducedDim(by.cluster, "PCA")
library(TSCAN)
mst <- createClusterMST(centroids, clusters=NULL)
mst
#Dibujando el árbol en el TSNE
line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="TSNE")

plotTSNE(sce.nest, colour_by="label") + 
  geom_line(data=line.data, mapping=aes(x=dim1, y=dim2, group=edge))
#El siguiente paso es complicado ya que la selección del flujo de este
#pseudotiempo o trayectoria, se debe elegir el punto de partida o el punto más
#basal utilizando conocimiento biológico sobre el estado de pluripotencia de 
#cada grupo celular representado en los clusters. En este ejercicio se eligio
#uno de los puntos finales del arbol como punto de inicio.
map.tscan <- mapCellsToEdges(sce.nest, mst=mst, use.dimred="PCA")
tscan.pseudo <- orderCells(map.tscan, mst)
head(tscan.pseudo)
#Dibujando la trayectoria temporal sobre el TSNE
common.pseudo <- averagePseudotime(tscan.pseudo) 
plotTSNE(sce.nest, colour_by=I(common.pseudo), 
         text_by="label", text_colour="red") +
  geom_line(data=line.data, mapping=aes(x=dim1, y=dim2, group=edge))
#Caracterización de las trayectorias y lo cambios de los genes a traves del pseudotiempo
library(TSCAN)
pseudo<- testPseudotime(sce.nest, pseudotime = tscan.pseudo[,1])[[1]]
pseudo$SYMBOL <- rowData(sce.nest)$SYMBOL
sorted<-pseudo[order(pseudo$p.value),]
sce.nest$TSCAN.first <- pathStat(tscan.pseudo)[,1]
up.left <- sorted[sorted$logFC < 0,]
head(up.left, 10)
best<- head(up.left$SYMBOL,10)
plotExpression(sce.nest, features = best, swap_rownames = "SYMBOL",
               x ="TSCAN.first", colour_by = "label")