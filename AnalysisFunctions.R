detach("package:plyr", unload=TRUE); library(clusterProfiler); library(org.Hs.eg.db);library(mygene);library(ggplot2)
library(STRINGdb);library(igraph);library(dplyr); library(heatmap3); library(ggrepel); library(drc); library(ggExtra)

QueryKEGG <- function(genes,pvalue = 0.05,padjust = "bonferroni",keyType =  'uniprot'){
  "This function takes in a vector of gene names, converts them to Entrez gene symbols and queries KEGG"
  "It returns a dataframe containing all the pathways, and their pvalues/FDRs."
  EG_IDs <- list()
  mapping <- list()
  for(i in 1:length(genes)){EG_IDs[i] <- mygene::query(genes[i])$hits$entrezgene[1]}
  for(i in 1:length(genes)){mapping[as.character(EG_IDs[[i]])] <- as.character(genes[i])}
  pathways <- enrichKEGG(EG_IDs, pvalueCutoff = pvalue, pAdjustMethod = padjust)
  pathways <- as.data.frame(pathways)
  temp <- c()
  for(i in 1:length(pathways$geneID)){
    entrez <- strsplit(pathways$geneID[i],split="/")
    temp2 <- c()
    for(j in 1:length(entrez[[1]])){temp2 <- c(temp2,mapping[entrez[[1]][j]])}
    temp <- c(temp,paste(temp2,collapse="/"))}
  pathways$genes <- temp
  return(pathways)
}
QueryGO <- function(genes,ont = "bp",pvalue = 0.05,padjust = "bonferroni",keyType =  'uniprot'){
  "This function takes in a vector of gene names, converts them to Entrez gene symbols and queries KEGG"
  "It returns a dataframe containing all the pathways, and their pvalues/FDRs."
  EG_IDs <- list()
  mapping <- list()
  for(i in 1:length(genes)){EG_IDs[i] <- mygene::query(genes[i])$hits$entrezgene[1]}
  for(i in 1:length(genes)){mapping[as.character(EG_IDs[[i]])] <- as.character(genes[i])}
  pathways <- enrichGO(EG_IDs,OrgDb = 'org.Hs.eg.db',ont = ont, pvalueCutoff = pvalue, pAdjustMethod = padjust)
  pathways <- as.data.frame(pathways)
  temp <- c()
  for(i in 1:length(pathways$geneID)){
    entrez <- strsplit(pathways$geneID[i],split="/")
    temp2 <- c()
    for(j in 1:length(entrez[[1]])){temp2 <- c(temp2,mapping[entrez[[1]][j]])}
    temp <- c(temp,paste(temp2,collapse="/"))}
  pathways$genes <- temp
  return(pathways)
}
PlotKEGG <- function(KEGGpathways,pval = "p.adjust",pvalueCutoff,top10 = TRUE,col="black",theme = "classic"){
  "This function takes in the KEGG dataframe from QueryKEGG and plots the data."
  KEGGpathways$x <- factor(KEGGpathways$Description,levels=rev(KEGGpathways$Description))
  if(top10 == TRUE) {KEGGpathways <- KEGGpathways[1:10,]}
  if(pval=="pvalue"){
    p <- ggplot(data=KEGGpathways, aes(y=-log10(pvalue), x=x)) + 
      geom_bar(stat="identity",fill=col) + coord_flip()}
  if(pval=="p.adjust"){
    p <- ggplot(data=KEGGpathways, aes(y=-log10(p.adjust), x=x)) + 
      geom_bar(stat="identity",fill=col) + coord_flip()}
  
  if(theme== "classic") {p <- p + theme_classic()}
  if(theme== "b/w") {p <- p + theme_bw()}
  if(theme== "minimal") {p <- p + theme_minimal()}
  if(theme== "dark") {p <- p + theme_dark()}
  if(theme== "linedraw") {p <- p + theme_linedraw()}
  p <- p + xlab("")+theme(axis.title.y = element_text(size=rel(1.5),face="bold"),
                          axis.title.x = element_text(size=rel(1.5),face="bold"),
                          axis.text.x = element_text(size=rel(1.5),face="bold"),
                          axis.text.y = element_text(size=rel(1.5),face="bold"),
                          strip.text.x = element_text(size=rel(1.5),face="bold"),
                          legend.text = element_text(face="bold"),
                          legend.title = element_text(face="bold"))
  p
}
PlotOntology <- function(Ontology,pval = "p.adjust",pvalueCutoff,top10 = TRUE,col="black",theme = "classic"){
  "This function takes in the GO dataframe from QueryGO and plots the data."
  Ontology$x <- factor(Ontology$Description,levels=rev(Ontology$Description))
  if(top10 == TRUE) {Ontology <- Ontology[1:10,]}
  if(pval=="pvalue"){
    p <- ggplot(data=Ontology, aes(y=-log10(pvalue), x=x)) + 
      geom_bar(stat="identity",fill=col) + coord_flip()}
  if(pval=="p.adjust"){
    p <- ggplot(data=Ontology, aes(y=-log10(p.adjust), x=x)) + 
      geom_bar(stat="identity",fill=col) + coord_flip()}
  
  if(theme== "classic") {p <- p + theme_classic()}
  if(theme== "b/w") {p <- p + theme_bw()}
  if(theme== "minimal") {p <- p + theme_minimal()}
  if(theme== "dark") {p <- p + theme_dark()}
  if(theme== "linedraw") {p <- p + theme_linedraw()}
  p <- p + xlab("")+theme(axis.title.y = element_text(size=rel(1.5),face="bold"),
                          axis.title.x = element_text(size=rel(1.5),face="bold"),
                          axis.text.x = element_text(size=rel(1.5),face="bold"),
                          axis.text.y = element_text(size=rel(1.5),face="bold"),
                          strip.text.x = element_text(size=rel(1.5),face="bold"),
                          legend.text = element_text(face="bold"),
                          legend.title = element_text(face="bold"))
  p
}
QuerySTRING <- function(x,column,score_threshold = 800,header=TRUE){
  "This takes in a dataframe and column name (string) and queries STRING. It returns a SIF formatted dataframe."
  string_db <- STRINGdb$new( version="10", species=9606, score_threshold=score_threshold, input_directory="" )
  mapped <- string_db$map(x, column, removeUnmappedRows = TRUE )
  mapped <- unique(mapped)
  mapped$Gene <- as.character(mapped$Gene)
  map_dict <- list()
  for(i in 1:length(mapped$Gene)){
    map_dict[mapped$STRING_id[i]] <- mapped$Gene[i]
  }
  network <- string_db$get_interactions(mapped$STRING_id)
  network <- as.data.frame(network)
  to = c()
  from = c()
  for(i in 1:length(network[,1])){
    to <- c(to,unlist(map_dict[network$to[i]])[1])
    from <- c(from,unlist(map_dict[network$from[i]])[1])
  }
  sif <- data.frame(Source = to, Target = from)
  if(header==FALSE){colnames(sif) <- NULL}
  return(sif)
}
NetworkAttributes <- function(edges,centrality = "eigenvector",community="fast.greedy") {
  "This function takes in a SIF formatted dataframe and returns a node attribute dataframe."
  "It performs community optimization and calculates a centrality metric."
  nodes <- data.frame(V1 = unique(c(as.character(edges[,1]),as.character(edges[,2]))))
  graph <- graph.data.frame(edges, directed = FALSE, vertices = nodes)
  if(community == "optimal"){V(graph)$comm <- membership(optimal.community(graph))} # computationally intensive
  if(community == "fast.greedy"){V(graph)$comm <- membership(fastgreedy.community(graph))} # for larger datasets\
  if(community == "edge.betweenness"){V(graph)$comm <- membership(edge.betweenness.community(graph))}
  if(community == "walk.trap"){V(graph)$comm <- membership(walktrap.community(graph))}
  if(community == "spin.glass"){V(graph)$comm <- membership(spinglass.community(graph))}
  if(community == "leading.eigenvector"){V(graph)$comm <- membership(leading.eigenvector.community(graph))}
  if(community == "label.propagation"){V(graph)$comm <- membership(label.propagation.community(graph))}
  if(community == "multilevel"){V(graph)$comm <- membership(multilevel.community(graph))}
  
  if(centrality == "closeness"){V(graph)$closeness <- centralization.closeness(graph)$res}
  if(centrality == "betweenness"){V(graph)$betweenness <- centralization.betweenness(graph)$res}
  if(centrality == "eigenvector"){V(graph)$eigen <- centralization.evcent(graph)$vector}
  if(centrality == "PageRank"){V(graph)$page <- page_rank(graph)$vector}
  
  return(get.data.frame(graph, what = "vertices"))
}
CorrelationPlot <- function(x1, y1, xlab="", ylab="",pt_col = "black",xhist_col="gray",yhist_col = "gray",log2 = TRUE){
  "This function plots a scatterplot with histograms on the margins"
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  if(log2 == TRUE){
    xhist = hist(log2(x1), plot=FALSE)
    yhist = hist(log2(y1), plot=FALSE)
  }
  if(log2 == FALSE){
    xhist = hist(x1, plot=FALSE)
    yhist = hist(y1, plot=FALSE)
  }
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  if(log2 == TRUE){plot(log2(x1),log2(y1),col=pt_col)}
  if(log2 == FALSE){plot(x1,y1,col=pt_col)}
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0,col=xhist_col)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE,col=yhist_col)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=.8 * (mean(x1) - min(x1))/(max(x1)-min(x1)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=(.8 * (mean(y1) - min(y1))/(max(y1) - min(y1))))
  print(summary(lm(y1~x1)))
}
gg.addHistogram <- function(p,hist.col="gray"){
  ggMarginal(p, type = "histogram", xparams = list(binwidth = 1, fill = hist.col))
}
VolcanoPlot_andFilter <- function(x,y,xcutoffs = c(-2,2),ycutoff,cols = c("black","red"),xlab = "",ylab = "",xlim="default"){
  data <- data.frame(x=x,y=y)
  x1 <- subset(data,x >= xcutoffs[2]) ; x1 <- subset(x1, y <= ycutoff)
  x2 <- subset(data, x <= xcutoffs[1]); x2 <- subset(x2, y <= ycutoff)
  if(xlim != "default"){plot(data$x,-log10(data$y),xlab=xlab,ylab=ylab,col=cols[1],xlim=xlim)} else {plot(data$x,-log10(data$y),xlab=xlab,ylab=ylab,col=cols[1])}
  abline(v=xcutoffs[1],lty=2); abline(v=xcutoffs[2],lty=2);abline(h=-log10(ycutoff),lty=2)
  points(x1$x,-log10(x1$y),col=cols[2])
  points(x2$x,-log10(x2$y),col=cols[2])
}
VolcanoPlot_orFilter <- function(x,y,xcutoffs = c(-2,2),ycutoff,cols = c("black","red"),xlab = "",ylab = "",xlim="default"){
  data <- data.frame(x=x,y=y)
  x1 <- subset(data,x >= xcutoffs[2])
  x2 <- subset(data, x <= xcutoffs[1]); x3 <- subset(data, y <= ycutoff)
  if(xlim != "default"){plot(data$x,-log10(data$y),xlab=xlab,ylab=ylab,col=cols[1],xlim=xlim)} else {plot(data$x,-log10(data$y),xlab=xlab,ylab=ylab,col=cols[1])}
  abline(v=xcutoffs[1],lty=2); abline(v=xcutoffs[2],lty=2);abline(h=-log10(ycutoff),lty=2)
  points(x1$x,-log10(x1$y),col=cols[2])
  points(x2$x,-log10(x2$y),col=cols[2])
  points(x3$x,-log10(x3$y),col=cols[2])
}
ModulePlot <- function(Attributes,Module,pSTY_data, col="black",theme = "classic",plot="FoldChange",identifier = "ID",
                       showColDendro = FALSE,cexRow = 1, margins = c(5,5), Colv = "default"){
  "This function takes in the GO dataframe from QueryGO and plots the data."
  colnames(Attributes) <- c("Gene","comm","eigen")
  final <- merge(pSTY_data,Attributes,all=TRUE)
  final <- subset(final,comm == Module)
  if(plot=="FoldChange"){
    final <- arrange(final,desc(FoldChange))
    if(identifier == "ID"){
      p <- ggplot(data=final, aes(y=log2(FoldChange), x=ID)) + geom_bar(stat="identity",fill=col) + coord_flip() + geom_hline(yintercept=0)
    } 
    if(identifier == "ID2"){
      p <- ggplot(data=final, aes(y=log2(FoldChange), x=ID2)) + geom_bar(stat="identity",fill=col) + coord_flip() + geom_hline(yintercept=0)
    } else {
      p <- ggplot(data=final, aes(y=log2(FoldChange), x=Gene)) + geom_bar(stat="identity",fill=col) + coord_flip() + geom_hline(yintercept=0)
    }
  }
  if(plot=="eigen"){
    final <- arrange(final,desc(eigen))
    if(identifier == "ID"){
      p <- ggplot(data=final, aes(y=eigen, x=ID)) + geom_bar(stat="identity",fill=col) + coord_flip() + geom_hline(yintercept=0)
    } 
    if(identifier == "ID2"){
      p <- ggplot(data=final, aes(y=eigen, x=ID2)) + geom_bar(stat="identity",fill=col) + coord_flip() + geom_hline(yintercept=0)
    } else {
      p <- ggplot(data=final, aes(y=eigen, x=Gene)) + geom_bar(stat="identity",fill=col) + coord_flip() + geom_hline(yintercept=0)
    }
  }
  if(plot != "heatmap"){
    if(theme== "classic") {p <- p + theme_classic()}
    if(theme== "b/w") {p <- p + theme_bw()}
    if(theme== "minimal") {p <- p + theme_minimal()}
    if(theme== "dark") {p <- p + theme_dark()}
    if(theme== "linedraw") {p <- p + theme_linedraw()}
    p <- p + xlab("")+theme(axis.title.y = element_text(size=rel(1),face="bold"),
                            axis.title.x = element_text(size=rel(1),face="bold"),
                            axis.text.x = element_text(size=rel(1),face="bold"),
                            axis.text.y = element_text(size=rel(1),face="bold"),
                            strip.text.x = element_text(size=rel(1),face="bold"),
                            legend.text = element_text(face="bold"),
                            legend.title = element_text(face="bold"))
    return(p)
  }
  if(plot == "heatmap"){
    matrix <- as.matrix(final[,8:13])
    if(identifier == "ID2"){row.names(matrix) <- final$ID2}
    if(identifier == "Gene"){row.names(matrix) <- final$Gene}
    if(identifier == "ID"){row.names(matrix) <- final$ID}
    if(Colv == "NA"){heatmap3(matrix,showColDendro = showColDendro,cexRow = cexRow,Colv = NA, margins = margins)} else {
      heatmap3(matrix,showColDendro = showColDendro,cexRow = cexRow, margins = margins)}
  }
}
PlotNetwork <- function(edges,node_attributes,type = "ForceDirected", open = "viewer",saveHTML=FALSE){
  colnames(edges) <- NULL
  write.table(edges,file="tmpEdge1.txt",row.names = FALSE,quote = FALSE,sep="\t")
  write.table(node_attributes,file = "tmpNode1.txt",row.names = FALSE,quote = FALSE,sep="\t")
  if(type == "ForceDirected"){arg1 = "plotForce.py"}
  if(type == "EV-ForceDirected"){arg1 = "plotForce_EV.py"}
  if(type == "EdgeBundle"){arg1 = "plot_edge_bundling.py"}
  if(type == "Adjacency"){arg1 = "plot_adjacency.py"}
  system2(command = "python",args = c(arg1,"tmpNode1.txt","tmpEdge1.txt"),stdout=TRUE)
  file.remove("tmpNode1.txt");file.remove("tmpEdge1.txt")
  tempDir <- tempfile()
  dir.create(tempDir)
  viewer <- getOption("viewer")
  if(open == "viewer"){
    if(type == "EdgeBundle"){file.copy(from = "Edge_bundled_network.html", to = file.path(tempDir)); htmlFile <- file.path(tempDir,"Edge_bundled_network.html")}
    if(type == "Adjacency"){file.copy(from = "Adjacency.html", to = file.path(tempDir)); htmlFile <- file.path(tempDir,"Adjacency.html")}
    if(type == "ForceDirected"){file.copy(from = "Force_directed_network.html", to = file.path(tempDir)); htmlFile <- file.path(tempDir,"Force_directed_network.html")}
    if(type == "EV-ForceDirected"){file.copy(from = "EV_Force_directed_network.html", to = file.path(tempDir)); htmlFile <- file.path(tempDir,"EV_Force_directed_network.html")}
    viewer(htmlFile)}
  if(open == "browser"){
    if(type == "EdgeBundle"){file.copy(from = "Edge_bundled_network.html", to = file.path(tempDir)); htmlFile <- file.path(tempDir,"Edge_bundled_network.html")}
    if(type == "Adjacency"){file.copy(from = "Adjacency.html", to = file.path(tempDir)); htmlFile <- file.path(tempDir,"Adjacency.html")}
    if(type == "ForceDirected"){file.copy(from = "Force_directed_network.html", to = file.path(tempDir)); htmlFile <- file.path(tempDir,"Force_directed_network.html")}
    if(type == "EV-ForceDirected"){file.copy(from = "EV_Force_directed_network.html", to = file.path(tempDir)); htmlFile <- file.path(tempDir,"EV_Force_directed_network.html")}
    browseURL(htmlFile)}
  if(saveHTML == FALSE){
    if(type == "EdgeBundle"){file.remove("Edge_bundled_network.html")}
    if(type == "Adjacency"){file.remove("Adjacency.html")}
    if(type == "ForceDirected"){file.remove("Force_directed_network.html")}
    if(type == "EV-ForceDirected"){file.remove("EV_Force_directed_network.html")}}
  
}
FastWrite <- function(table,file){write.table(table,file,row.names = FALSE,quote = FALSE,sep="\t")}
FastRead <- function(file){read.table(file,header = TRUE, sep="\t",stringsAsFactors = FALSE)}

CompareModule <- function(edges,community = "fast.greedy",mode = "hclust",type="dendrogram") {
  nodes <- data.frame(V1 = unique(c(as.character(edges[,1]),as.character(edges[,2]))))
  graph <- graph.data.frame(edges, directed = FALSE, vertices = nodes)
  if(community == "optimal"){comm <- optimal.community(graph)} # computationally intensive
  if(community == "fast.greedy"){comm <- fastgreedy.community(graph)} # for larger datasets\
  if(community == "edge.betweenness"){comm <- edge.betweenness.community(graph)}
  if(community == "walk.trap"){comm <- walktrap.community(graph)}
  if(community == "spin.glass"){comm <- spinglass.community(graph)}
  if(community == "leading.eigenvector"){comm <- leading.eigenvector.community(graph)}
  if(community == "label.propagation"){comm <- label.propagation.community(graph)}
  if(community == "multilevel"){comm <- multilevel.community(graph)}
  if(type=="dendrogram"){plot_dendrogram(comm)}
  if(type=="graph"){modularity(comm);membership(comm);plot(comm,graph)}
}
Function2Network <- function(pathway_DF,data,KEGG_IDs,score_threshold = 800,header=TRUE){
  master_list <- c()
  for(i in 1:length(KEGG_IDs)){master_list <- c(master_list, subset(pathway_DF,ID==KEGG_IDs[i])$genes)}
  genes <- paste(master_list,collapse="/")
  genes <- unlist(strsplit(genes,split="/"))
  output_data <- subset(data, Gene %in% genes)
  network <- QuerySTRING(as.matrix(output_data),"Accession",score_threshold = score_threshold , header=header)
  return(network)
}
QueryGene <- function(DF,column = "Accession",colname = "Gene"){
  write.table(DF,file="tmpDF1.txt",row.names = FALSE,quote = FALSE,sep="\t")
  system2(command = "python",args = c("GetGeneName.py","tmpDF1.txt",column, colname),stdout=TRUE)
  file.remove("tmpDF1.txt")
  output <- read.table("tmpDF2.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)
  file.remove("tmpDF2.txt")
  return(output)
}
plot.IC50 <- function(DF,size=7,xlab="IC50",col){
  "This function plots IC50 values. Requires DF in following format:
  Drug | Protein1 | Protein2 | [...]
  ID1  |  value1  |  value2  | [...]
  "
  require(ggplot2);require(ggrepel);require(reshape2)
  p <- ggplot(melt(DF)) + geom_point(aes(x=value,y=Drug,color=Drug),size=size) + ylab("") + 
    log10_x_sci() + annotation_logticks(base=10,side="tb") + xlab(xlab) + geom_text_repel(aes(x=value,y=Drug,label=variable),color="black")
  if(!missing(col)){p <- p + scale_color_manual(values=col)}
  p
}
theme_matplotlib <- function(){theme(axis.title.y = element_text(size=rel(1.5),face="bold"),
                                     axis.title.x = element_text(size=rel(1.5),face="bold"),
                                     axis.text.x = element_text(size=rel(1.5),face="bold"),
                                     axis.text.y = element_text(size=rel(1.5),face="bold"),
                                     strip.text.x = element_text(size=rel(1.5),face="bold"),
                                     legend.text = element_text(face="bold"),
                                     legend.title = element_text(face="bold"),
                                     panel.border = element_rect(colour = "black", fill=NA, size=1),
                                     panel.background = element_blank(),
                                     panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank(),
                                     axis.ticks.x=element_blank(),
                                     axis.ticks.y=element_blank(),
                                     strip.background = element_blank(),
                                     legend.background = element_blank(),
                                     legend.key = element_blank())}
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
not_fancy <- function(l) {
  l <- format(l, scientific = FALSE)
  parse(text=l)
}
human_numbers <- function(x = NULL, smbl ="", signif = 1){
  humanity <- function(y){
    
    if (!is.na(y)){
      tn <- round(abs(y) / 1e12, signif)
      b <- round(abs(y) / 1e9, signif)
      m <- round(abs(y) / 1e6, signif)
      k <- round(abs(y) / 1e3, signif)
      
      if ( y >= 0 ){
        y_is_positive <- ""
      } else {
        y_is_positive <- "-"
      }
      
      if ( k < 1 ) {
        paste0( y_is_positive, smbl, round(abs(y), signif ))
      } else if ( m < 1){
        paste0 (y_is_positive, smbl,  k , "k")
      } else if (b < 1){
        paste0 (y_is_positive, smbl, m ,"m")
      }else if(tn < 1){
        paste0 (y_is_positive, smbl, b ,"bn")
      } else {
        paste0 (y_is_positive, smbl,  comma(tn), "tn")
      }
    } else if (is.na(y) | is.null(y)){
      "-"
    }
  }
  
  sapply(x,humanity)
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
log10_x_sci <- function(){scale_x_log10(labels = fancy_scientific)}
log10_y_sci <- function(){scale_y_log10(labels=fancy_scientific)}
log10_x_human <- function(){scale_x_log10(labels = human_numbers)}
log10_y_human <- function(){scale_y_log10(labels = human_numbers)}
log10_x_norm <- function(){scale_x_log10(labels = not_fancy)}
log10_y_norm <- function(){scale_y_log10(labels= not_fancy)}
human_gbp   <- function(x){human_numbers(x, smbl = "£")}
human_usd   <- function(x){human_numbers(x, smbl = "$")}
human_euro  <- function(x){human_numbers(x, smbl = "€")} 
human_num <- function(x){human_numbers(x, smbl = "")}
##### Base Plotting #####
log10.axis <- function(side, at, ...) {
  at.minor <- log10(outer(1:9, 10^(min(at):max(at))))
  lab <- sapply(at, function(i) as.expression(bquote(10^ .(i))))
  axis(side=side, at=at.minor, labels=NA, tcl=par("tcl")*0.5, ...)
  axis(side=side, at=at, labels=lab, ...)}
shaded.DRC.lines <- function(ResponseData,curve.col=c("#FF7878","#AEC6CF"),curve.pch=c(0,2),
                             shade.col="#d3d3d3",bar.col = c("#FF7878","#AEC6CF"),AUCplot=TRUE,
                             legendpos="default",bar.ylim=NULL,curve.ylim=c(0,120),curve.xlim=NULL,
                             conc.range=FALSE,conc.col="#ffb347",conc.line.col="orange",
                             conc.pos=110,curve.xlab="Dose",curve.ylab="Response",bar.ylab="AUC",
                             curve.main=NULL,bar.main=NULL,vehicle = "DMSO"){
  # Accepts DRC of 2 drug set and returns a plot and AUC values for the curve
  # Can annotate a concentration range on the plot
  # Can also plot AUCs alongside of the curve
  # AUCs are returned in case the user wants to plot separately
  # Input is dataframe with columns Dose - Response - Drug
  # Must only contain 2 drugs (single agent and combination)
  require(dplyr)
  require(shape)
  require(pracma) # calculate AUC of the curves
  data2 <- ResponseData %>%
    group_by(Drug, Dose) %>%
    summarise(std = sd(Response),Response = mean(Response))
  
  drugs = unique(data2$Drug)
  AUC <- c()
  cnt = 1
  for(i in 1:length(drugs)){
    temp = subset(data2, Drug == drugs[i])
    temp[temp$Dose== 0,"Dose"] <- min(temp$Dose[temp$Dose>0])/10
    AUC <- c(AUC, trapz(temp$Dose,temp$Response))
  }
  names(AUC) <- drugs
  AUC <- AUC[order(AUC,decreasing = TRUE)] # force into higher AUC first
  ord <- names(AUC)
  for(i in 1:length(ord)){
    temp = subset(data2, Drug == ord[i])
    temp[temp$Dose== 0,"Dose"] <- min(temp$Dose[temp$Dose>0])/10
    if(i==1){
      u <- par("usr")
      rect(u[1], u[3], u[2], u[4], col = "red")
      if(AUCplot==TRUE){par(mfrow=c(1,2))}
      graphics::plot(temp$Dose,temp$Response,col=curve.col[i],pch=curve.pch[i],lwd=2,
                     type="p",log="x",ylim=curve.ylim,xlab=curve.xlab,ylab=curve.ylab,
                     xaxt='n',panel.last=TRUE,yaxt="n",curve.xlim,main = curve.main)
      axis(side=2)
      ax <- as.character(sort(temp$Dose)); ax[1] <- vehicle
      axis(1, at=temp$Dose, labels=ax,tick = FALSE)
    } else {points(temp$Dose,temp$Response,col=curve.col[i],pch=curve.pch[i],lwd=2)}
    if(cnt == 1){polygon(c(temp$Dose,rev(temp$Dose)),c(rep(0,times=length(temp$Response)),rev(temp$Response)),col=shade.col,border=NA)} # upper curve}
    if(cnt == 2){polygon(c(temp$Dose,rev(temp$Dose)),c(rep(0,times=length(temp$Response)),rev(temp$Response)),col='white',border=NA)} # lower curve}
    cnt <- cnt+1
    if(conc.range[1] != FALSE){
      polygon(x=c(c(conc.range[1],conc.range[1]),rev(c(conc.range[2],conc.range[2]))),
              y=c(curve.ylim,rev(curve.ylim)),border=NA,col=adjustcolor(conc.col,alpha.f=0.3))
      abline(v=conc.range[1],lwd=2,col=conc.line.col)
      abline(v=conc.range[2],lwd=2,col=conc.line.col)
      Arrows(x0=conc.range[1],x1=conc.range[2],y0=conc.pos,y1=conc.pos,col="black",code=3)
    }
    points(temp$Dose,temp$Response,col=curve.col[i],pch=curve.pch[i],lwd=2)
    arrows(temp$Dose, temp$Response-temp$std, temp$Dose, temp$Response+temp$std, length=0.05, angle=90, code=3,col=curve.col[i],lwd=2)
    lines(temp$Dose,temp$Response,col=curve.col[i],lwd=2)
  }
  axis(side=1, at=c(pretty(c(0,0.001),n=10), # need to automate this with a function
                    pretty(c(0.001,0.01),n=10),
                    pretty(c(0.01,0.1),n=10),
                    pretty(c(0.1,1),n=10),
                    pretty(c(1,10),n=10),
                    pretty(c(10,100),n=10)
  ),
  labels=FALSE)
  abline(h=0,lty=2)
  abline(v=min(temp$Dose),lty=2)
  if(legendpos[1] == "default"){
    legend(x=min(temp$Dose)*1.05,y=40,
           legend=ord,col=curve.col[1:length(ord)],
           pch=curve.pch[1:length(ord)],lty=rep(1,times=length(ord)),
           bty="n",lwd=rep(2,times=length(ord)),bg="transparent")}
  else {legend(x=legendpos[1],y=legendpos[2],
               legend=ord,col=curve.col[1:length(ord)],
               pch=curve.pch[1:length(ord)],lty=rep(1,times=length(ord)),
               bty="n",lwd=rep(2,times=length(ord)),bg="transparent")}
  if(AUCplot==TRUE){
    barplot(height=AUC,names.arg=ord,col = bar.col,ylab="AUC",ylim=bar.ylim,main=bar.main)
    par(mfrow=c(1,1))}
  return(AUC)
}
baseDRC <- function(data,fct=LL.3(),vehicle = "DMSO",lines=FALSE, estimates = c(50), ylab = "Viability",
                      xlab = "Conc", xlim = c(0,20), ylim = c(0,120),
                      col = c("red","blue","black","dark green","red","blue","black","dark green"),
                      pch = c(0,2,1,6,6,1,2,0),
                      legendPos = c(10,100)){
  require(drc);require(multcomp);require(lmtest);require(sandwich)
  colors=col
  data2 <- data %>%
    group_by(Drug, Dose) %>%
    summarise(std = sd(Response),Response = mean(Response))
  curve <- drm(Response ~ Dose,curveid = Drug, data = data2, fct = fct)
  #summary(curve)
  #coeftest(curve, vcov = sandwich)
  #summary(glht(curve))
  drugs = unique(data2$Drug)
  pl <- plot(curve, broken = FALSE, type = "obs",log="x",xlim=xlim,ylim=ylim,
             xlab = xlab,xttrim=FALSE,conName=vehicle,
             ylab = ylab,col=colors[1:length(drugs)],lty=rep(1,length(drugs)),
             lwd=rep(2,length(drugs)),legendPos = legendPos,
             pch=pch[1:length(drugs)])
  if(lines == TRUE){
    plot(curve, broken = FALSE, type = "obs",log="x",xlim=xlim,ylim=ylim,
         xlab = xlab,xttrim=FALSE,conName=vehicle,
         ylab = ylab,col=colors[1:length(drugs)],lty=rep(1,length(drugs)),
         lwd=rep(2,length(drugs)),legendPos = legendPos,
         pch=pch[1:length(drugs)])
  } else {
    plot(curve, broken = FALSE, type = "all",log="x",xlim=xlim,ylim=ylim,
         xlab = xlab,xttrim=FALSE,conName=vehicle,
         ylab = ylab,col=colors[1:length(drugs)],lty=rep(1,length(drugs)),
         lwd=rep(2,length(drugs)),legendPos = legendPos,
         pch=pch[1:length(drugs)])
  }
  cnt = 1
  for(i in 1:length(drugs)){
    temp = subset(data2, Drug == drugs[i])
    arrows(temp$Dose, temp$Response-temp$std, temp$Dose, temp$Response+temp$std, length=0.05, angle=90, code=3,col=colors[i],lwd=2)
    if(lines==TRUE){
      lines(temp$Dose,temp$Response,col=colors[i],lwd=2)
    }}
  axis(side=1, at=c(pretty(c(0,0.001),n=10), # need to automate this with a function
                    pretty(c(0.001,0.01),n=10),
                    pretty(c(0.01,0.1),n=10),
                    pretty(c(0.1,1),n=10),
                    pretty(c(0.1,10),n=10),
                    pretty(c(10,100),n=10)),labels=FALSE)
  ED(curve, estimates, interval = "delta")
  invisible(pl)
}
ggDRC <- function(data,fct=LL.3(),col=NULL,size=3,xlab="Dose",ylab="Response",estimates = c(50),lines = FALSE,xlim = c(0,10),ylim=c(0,120)){
  data2 <- data %>%
    group_by(Drug, Dose) %>%
    summarise(std = sd(Response),Response = mean(Response))
  curve <- drm(Response ~ Dose,curveid = Drug, data = data2, fct = fct)
  pl <- baseDRC(data,fct=fct,estimates=estimates)
  dev.off()
  drugs <- c()
  response <- c()
  for(i in 1:length(unique(data2$Drug))){
    drugs <- c(drugs,rep(unique(as.character(data2$Drug))[i], times=length(pl$Dose)))
    response <- c(response,pl[,i+1])
  }
  pl <- data.frame(Dose = rep(pl$Dose,times=length(unique(data2$Drug))), Drug = drugs, Response = response)
  data2 <- data.frame(data2)
  if(lines == FALSE){
  p <- ggplot(data=data2, aes(x=Dose,y=Response))+ geom_hline(yintercept = 0,lty=2) + geom_point(data=data2, aes(x=Dose,y=Response,color=Drug,shape=Drug),size=size) + 
    geom_errorbar(data=data2,aes(x=Dose,ymin=Response-std,ymax=Response+std,color=Drug),width=.1) + xlab(xlab) + ylab(ylab) + 
    theme_matplotlib() + log10_x_sci() + geom_line(data=pl,aes(x=Dose,y=Response,color=Drug)) + annotation_logticks(base=10,sides="b")
  } else {
    p <- ggplot(data=data2, aes(x=Dose,y=Response))+ geom_hline(yintercept = 0,lty=2) + geom_point(data=data2, aes(x=Dose,y=Response,color=Drug,shape=Drug),size=size) + 
      geom_errorbar(data=data2,aes(x=Dose,ymin=Response-std,ymax=Response+std,color=Drug),width=.05) + xlab(xlab) + ylab(ylab) + 
      theme_matplotlib() + log10_x_sci() + geom_line(data=data2,aes(x=Dose,y=Response,color=Drug)) + annotation_logticks(base=10,sides="b")
  }
  p <- p + xlim(xlim[1],xlim[2]) + ylim(ylim[1],ylim[2])
  if(!is.null(col)){
    p <- p + scale_color_manual(values = col)
  }
  p
}
# Example data for shaded.DRC.lines()
Combination.sample <- read.table("ComboExample.txt",sep="\t",header=TRUE)
# Example data for ggDRC() and baseDRC()
DRC.sample <- read.table("drcExample.txt",sep="\t",header=TRUE)
# Example data for Network analysis
example.network <- read.table("edges.txt",sep="\t",header=FALSE)
# Example data for plot.IC50()
SRC.IC50 = data.frame(Drug = c("Dasatinib", "Bosutinib"),MAP2K1 = c(1000,19), MAP2K2 = c(1400,9.9),PTK2 = c(10000, 570))






