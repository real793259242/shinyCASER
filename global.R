library(shinythemes)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)
library(shinycssloaders)
library(dashboardthemes)
library(rintrojs)
library(reactable)
library(markdown)
library(ggpubr)
library(ggplot2)
library(ggExtra)
library(plotly)
library(corrplot)
library(Rgraphviz)
library(visNetwork)
library(dplyr)
library(sparkline)
library(foreach)
library(doParallel)
library(parallel)
library(tidyverse)
library(Seurat)
library(scran)
library(SingleCellExperiment)
library(fst)
library(olsrr)
library(tippy)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(echarty)

source("./dsepTest/function.R")

options(scipen = 999)
options(shinyApp.data = NULL)
options(shiny.maxRequestSize = 100*1024^2) # max file size: 100MB

print("Loading data ...")

load("./extdata/hsa_and_mmu_pathway_name.RData") # hsaAndMmuPathwayName
load("./extdata/EN-hsa-amat.RData") # hsa_EN_DAG_amatList
# load("./extdata/allGeneName.RData") # allGeneName


print("Loading complete!")

# numCores <- parallel::detectCores()

.citAlgo <- c("dcc.gamma", "dcc.perm", "hsic.gamma", 
              "hsic.perm", "hsic.clust", "cmiKnn", 
              "KCIT", "RCIT", "RCoT", "Gauss")

myColors <- c("#098ebb", "#fdc23a", "#e96449", "#818286", "#4E84C4")
plotly_fontsize <- 16

.CMIKNN <- reticulate::import("tigramite.independence_tests")

# .fs <- reticulate::import("pyFeatureSelection")
# .fsObj <- .fs$FeatureSelection()

# .fsAlgo <- c("Random Forest", "Extra Trees", "XGBoost", 
#              "BAHSIC", "SHS", 
#              "Lasso", "Ridge", "Elastic Net", "OLS regression")
# .cfsAlgo <- c("MMPC", "pc_simple", "HITON", "semi_HITON", "getPC", "MBtoPC")

# Close the window and shut down app----
jscode <- "shinyjs.closeWindow = function() { window.close(); }"


# Reset app
jsResetCode <- "shinyjs.reset = function() {history.go(0)}" # Define the js method that resets the page



#' Create a Seurat object with performing linearly and nonlinearly dimension reduction.
#' 
#' @param data A matrix-like object with cells as columns and features as rows.
#' @param meta.data Additional cell-level metadata to add to the Seurat object. Should 
#' be a data.frame where the rows are cell names and the columns are additional metadata 
#' fields. Row names in the metadata need to match the column names of the data matrix.
#' @param normalizeQ A logical type which asks whether to perform data normalization.
#' @param project Project name for the Seurat object.
#' @param min.cells Include features detected in at least this many cells.
#' @param find.HVGs.method The method of choosing top variable features.
#' @param num.HVGs Number of features to select as top variable features; only used when find.HVGs.method is set to 'dispersion' 
#' or 'vst'.
#' @param num.pc Total Number of PCs to compute and store (50 by default).
#' @param dims Which dimensions to use as input features.
#' 
#' @return A list to store Seurat object, sample names and feature names.
makeSeuratObject <- function(data, 
                             meta.data = NULL, 
                             normalizeQ, 
                             project, 
                             min.cells = 0, 
                             find.HVGs.method = "vst", 
                             num.HVGs = 2000, 
                             num.pc = 50, 
                             perplexity = 5, # for tsne, generally 5~50
                             dims = 1:5) 
{
    seuratObj <- Seurat::CreateSeuratObject(counts = data, 
                                            meta.data = meta.data, 
                                            project = project, 
                                            min.cells = min.cells)

    if (normalizeQ) # normalize
    {
        seuratObj <- NormalizeData(object = seuratObj)
    } else # do not normalize
    {
        seuratObj[["RNA"]]@data <- seuratObj[["RNA"]]@counts
    }

    if (ncol(data) >= num.HVGs)
    {
        seuratObj <- Seurat::FindVariableFeatures(seuratObj, 
                                                  selection.method = find.HVGs.method, 
                                                  nfeatures = num.HVGs)      
    } else
    {
        seuratObj <- Seurat::FindVariableFeatures(seuratObj, 
                                          selection.method = find.HVGs.method, 
                                          nfeatures = ncol(data))
    }


    # PCA
    all.genes <- rownames(seuratObj)
    seuratObj <- Seurat::ScaleData(seuratObj, features = all.genes)
    seuratObj <- Seurat::RunPCA(seuratObj, features = SeuratObject::VariableFeatures(object = seuratObj), npcs = num.pc)

    ### tSNE
    seuratObj <- Seurat::RunTSNE(seuratObj, dims = dims, perplexity = perplexity)

    ### UMAP ###
    seuratObj <- Seurat::RunUMAP(seuratObj, dims = dims)

    return(list(seuratObj = seuratObj, sampleName = colnames(data), featureName = rownames(data)))
}


myTooltip <- function(title, note, class = "help-tip")
{
  htmlCode <- paste0("<table><tr><td><strong>", title, "</strong></td>", 
                     "<td><div class=", class, "><p>", note, "</p></div></td></tr></table>")
  return(shiny::HTML(htmlCode))
}


genObserver_menus <- function(pat="btn_results_", n=1, updateVal) {
  res <- paste0('observeEvent(input$',pat,n,', {
    curid <- "',pat,n,'"
    nn <- names(input)
    nn <- nn[grep("',pat,'",nn)]
    nn <- setdiff(nn, curid)
    for (btnid in nn) {
      updateButton(session, btnid, style="default")
    }
    obj$',updateVal,' <- "',pat,n,'"
    updateButton(session, curid, style="primary")
  });
  ')
  res
}


strongTitleHTML <- function(title)
{
    return(HTML(paste0('<h2 style="text-align: center; margin-top: 0"><strong>', 
                       title, 
                       '</strong></h2>')))
}


myActionButton <- function(inputId, label, btn.style="", css.class="") {
  if ( btn.style %in% c("primary","info","success","warning","danger","inverse","link")) {
    btn.css.class <- paste("btn", btn.style, sep="-")
  } else {
    btn.css.class <- btn.style
  }
  tags$button(id=inputId, 
                type="button", 
                class=paste("btn action-button", btn.css.class, css.class, collapse=" "), 
                label)
}


myActionButton2 <- function(inputId, label, class="", btn.style = "default", style = "vertical-align:middle", type = "button")
{
    if ( btn.style %in% c("primary","info","success","warning","danger","inverse","link")) {
      btn.style <- paste("btn", btn.style, sep="-")
    }
    
    # id <- paste0('id="', inputId, '"')
    # class <- paste0('class="btn action-button ', btn.style, '"')
    # Style <- paste0('style="', style, '"')
    # Type <- paste0('type="', type, '"')
    # spanLabel <- paste0('<span>', label, '</span>')
    # HTML(paste0('<button ', id, ' ', class, ' ', Style, ' ', Type, '>', spanLabel, '</button>'))

    shinyBS::bsButton(inputId = inputId, 
                      label = shiny::span(label), 
                      class = class, 
                      style = btn.style)
}


miniButtonInBoxTitle <- function(inputId, label, title, icon = NULL)
{
    shiny::actionButton(inputId = inputId,
                        label = label,
                        title = title,
                        icon = icon,
                        class = "btn-xs",
                        style = "margin-right: 3px"
    )
}


warningAlerts <- function(txt)
{
    # https://getbootstrap.com/docs/3.4/components/#panels

    htmlCode <- paste0('<div class="alert alert-warning alert-dismissible" role="alert">
                          <button type="button" class="close" data-dismiss="alert" aria-label="Close"><span aria-hidden="true">&times;</span></button>
                          <strong>Warning!</strong> ', 
                        txt, 
                       '</div>')
    return(htmlCode)
}


# feature selection, Render a bar chart with positive and negative values
bar_chart_pos_neg <- function(label, value, max_value = 1, height = "16px",
                              pos_fill = "#005ab5", neg_fill = "#dc3220") {
  neg_chart <- div(style = list(flex = "1 1 0"))
  pos_chart <- div(style = list(flex = "1 1 0"))
  width <- paste0(abs(value / max_value) * 100, "%")

  if (value < 0) {
    bar <- div(style = list(marginLeft = "8px", background = neg_fill, width = width, height = height))
    chart <- div(style = list(display = "flex", alignItems = "center", justifyContent = "flex-end"), label, bar)
    neg_chart <- tagAppendChild(neg_chart, chart)
  } else {
    bar <- div(style = list(marginRight = "8px", background = pos_fill, width = width, height = height))
    chart <- div(style = list(display = "flex", alignItems = "center"), bar, label)
    pos_chart <- tagAppendChild(pos_chart, chart)
  }

  div(style = list(display = "flex"), neg_chart, pos_chart)
}


# Render a bar chart in the background of the cell
bar_style <- function(width = 1, fill = "#e6e6e6", height = "75%", align = c("left", "right"), color = NULL) {
  align <- match.arg(align)
  if (align == "left") {
    position <- paste0(width * 100, "%")
    image <- sprintf("linear-gradient(90deg, %1$s %2$s, transparent %2$s)", fill, position)
  } else {
    position <- paste0(100 - width * 100, "%")
    image <- sprintf("linear-gradient(90deg, transparent %1$s, %2$s %1$s)", position, fill)
  }
  list(
    backgroundImage = image,
    backgroundSize = paste("100%", height),
    backgroundRepeat = "no-repeat",
    backgroundPosition = "center",
    color = color
  )
}


# See the ?tippy documentation to learn how to customize tooltips
with_tooltip <- function(value, tooltip, ...) {
  div(style = "text-decoration: underline; text-decoration-style: dotted; cursor: help",
      tippy::tippy(value, tooltip, ...))
}


findNeighbors <- function(g, node)
{
  if(class(g) == "graphNEL")
  {
    g <- as(g, "matrix")
  }
  
  if(node %in% colnames(g))
  {
    parents <- colnames(g)[g[, node] == TRUE]
    children <- colnames(g)[g[node, ] == TRUE]
  }
  
  return(c(parents, children))
}


randomSubGraph <- function(amat, n_node) 
{
  v_names <- colnames(amat)

  while (1)
  {
    node <- sample(v_names, 1)
    neighbor <- findNeighbors(amat, node) # maybe has no neighborhood!
    if(length(neighbor) > 0)
    {
      subGraph_nodes <- c(node)
      for(i in 1:(n_node-1))
      {
        node <- sample(neighbor, 1)
        subGraph_nodes <- c(subGraph_nodes, node)
        tmp <- findNeighbors(amat, node)
        neighbor <- union(neighbor, tmp)
        neighbor <- setdiff(neighbor, subGraph_nodes)
      }
      break
    }
  }

  sub_amat <- amat[subGraph_nodes, subGraph_nodes]
  return(sub_amat)
}


# getAllSubgraphs <- function(amat, n_node)
# {
#   if(ggm::isAcyclic(amat) == FALSE){stop("The graph is not acyclic!")}
#   v_names <- colnames(amat)
#   cb <- combn(v_names, n_node)
#   c <- 0
#   l <- list()
#   for(i in 1:ncol(cb))
#   {
#     amat_sub <- amat[cb[, i], cb[, i]]
#     if(!(0 %in% (rowSums(amat_sub) + colSums(amat_sub))))
#     {
#       c <- c + 1
#       l[[c]] <- amat_sub
#     }
#   }
#   return(l)
# }


getAllConnectedSubgraphs <- function(amat, n_node)
{
    allPossibleNumberLimit <- 50000
    if(ggm::isAcyclic(amat) == FALSE){stop("The graph is not acyclic!")}
    v_names <- colnames(amat)
    cb <- combn(v_names, n_node)
    print(dim(cb))
    if (ncol(cb) > allPossibleNumberLimit)
    {
        cb <- cb[, sample(ncol(cb), allPossibleNumberLimit)]
    }

    l <- list()
    for(i in 1:ncol(cb))
    {
        cat(i, " ")
        amat_sub <- amat[cb[, i], cb[, i]]
        l[[i]] <- amat_sub
    }

    connectedQ <- unlist(lapply(l, function(x){graph::isConnected(as(x, "graphNEL"))}))
    connectedIndex <- which(connectedQ == TRUE)

    return(l[connectedIndex])
}


compoundAmatQ <- function(amat)
{
    # KEGG
    v <- rownames(amat)
    return(any(sapply(v, function(x){return(grepl(" ", x))})))
}


extractEdges <- function(AdjMat)
{
  # library(tidyverse)
  genes <- rownames(AdjMat)
  ind <- which(AdjMat == 1, arr.ind = TRUE)
  row_ind <- ind[, "row"]
  col_ind <- ind[, "col"]
  
  edges <- matrix(nrow=nrow(ind), ncol=3)
  colnames(edges) <- c("start", "->", "end")
  
  for(i in 1:nrow(ind))
  {
    edges[i, ] <- c(genes[row_ind[i]], "->", genes[col_ind[i]])
  }
  edges <- as.data.frame(edges)
  
  
  if (nrow(edges) > 1)
  {
    # check undirected edges !
    del <- c()
    for(i in 1:(nrow(edges)-1)) 
    {
      for(j in (i+1):nrow(edges)) 
      {
        if(edges$start[j] == edges$end[i]) 
        {
          if(edges$end[j] == edges$start[i]) 
          {
            edges$`->`[i] <- "<->"
            #edges$`->`[j] <- "<->"
            del <- c(del, j)
          }
        }
      }
    }
    
    edges <- edges[setdiff(rownames(edges), del), ]
    rownames(edges) <- 1:nrow(edges)
    edges2 <- matrix(nrow(edges), nrow(edges), 1)
    
    for(j in 1:nrow(edges))
    {
      edges2[j, ] <- paste(edges$start[j], edges$`->`[j], edges$end[j], sep = "")
    }
    
    edges2 <- as.data.frame(edges)
    
    edges_vec <- unlist(pmap(edges2, paste, sep=" "))
  } else
  {
    edges_vec <- unlist(pmap(edges, paste, sep=" "))
  }
  
  return(edges_vec)
}


edgeToAmat <- function(edge_vec)
{
  l <- lapply(strsplit(edge_vec, split = "->"), function(x) {unname(sapply(x, remove_punct))})
  node_vec <- Reduce(union, l)
  amat <- matrix(0, nrow = length(node_vec), ncol = length(node_vec))
  rownames(amat) <- colnames(amat) <- node_vec
  for (i in 1:length(l))
  {
    amat[l[[i]][1], l[[i]][2]] <- 1
  }
  return(amat)
}


generateDOTcode <- function(amat) 
{
    # node
    label <- letters[1:ncol(amat)]
    tmpTxt <- ""
    for (i in 1:length(label))
    {
        tmpTxt <- paste(tmpTxt, paste0(label[i], " [label = '@@", i, "']"), sep = "\n")
    }
    nodeTxt <- tmpTxt

    # edge
    edgePosition <- which(amat == 1, arr.ind = TRUE)
    tmpTxt <- ""
    for (j in 1:nrow(edgePosition))
    {
        index <- unname(edgePosition[j, ])
        tmpTxt <- paste(tmpTxt, paste0(label[index[1]], " -> ", label[index[2]]), sep = "\n")
    }
    edgeTxt <- tmpTxt
    # cat(edgeTxt)


    # footer
    tmpTxt <- ""
    for (k in 1:length(label))
    {
        tmpTxt <- paste(tmpTxt, paste0("[", k, "]: '", colnames(amat)[k], "'"), sep = "\n")
    }
    footerTxt <- tmpTxt
    # cat(footerTxt)

    txt <- paste('digraph a_nice_graph {', 
                 'node [footname = Helvetica]', 
                 nodeTxt, 
                 edgeTxt, 
                 '}', 
                 footerTxt, 
                 sep = "\n")
    # cat(txt)
    return(txt)
}


amatToDgrGraph <- function(amat, 
                           node_type = "lower", 
                           node_fontname = "Helvetica", 
                           node_fontsize = 20, 
                           node_fontcolor = "gray50", 
                           node_fillcolor = "aliceblue", 
                           node_fixedsize = TRUE, 
                           node_width = 0.5, 
                           node_style = "filled",
                           node_color = "gray70",
                           node_shape = "circle", 
                           edge_rel = "leading_to", 
                           edge_fontname = "Helvetica", 
                           edge_fontsize = 8, 
                           edge_len = 1.5, 
                           edge_color = "gray80", 
                           edge_arrowsize = 1, 
                           edge_arrowhead = "vee")
{
	ndf <- create_node_df(
			n = ncol(amat),
			label = colnames(amat),
			type = node_type, 
			fontname = node_fontname, 
			fontsize = node_fontsize, 
			fontcolor = node_fontcolor, 
			fillcolor = node_fillcolor, 
			fixedsize = node_fixedsize, 
			width = node_width, 
			style = node_style,
			color = node_color,
			shape = rep(node_shape, ncol(amat)))      

    tmpEdges <- as.data.frame(which(amat == 1, arr.ind = TRUE))      

    edf <- create_edge_df(
			from = tmpEdges$row,
			to = tmpEdges$col,
			rel = edge_rel, 
			fontname = edge_fontname, 
			fontsize = edge_fontsize, 
			len = edge_len, 
			color = edge_color, 
			arrowsize = edge_arrowsize, 
			arrowhead = edge_arrowhead)      

    graph <- create_graph(
				nodes_df = ndf,
				edges_df = edf)

    return(graph)
}



parseDsepClaim <- function(strvec)
{
    paste(paste(strvec[1:2], collapse = " "), " | { ", paste(strvec[-(1:2)], collapse = " "), " }", sep="")
}


makeDAGCandidateDF <- function(l)
{
    if (is.null(l)) {return(NULL)}
    # make candidate DAG list as dataframe
    res <- matrix(NA, nrow = length(l), ncol = 4)
    colnames(res) <- c("DAG", "#nodes", "#edges", "edges list")

    if (is.null(names(l)))
    {
        for (i in 1:length(l)) 
        {
            amat <- l[[i]]
            g <- as(amat, "graphNEL")
            res[i, ] <- c(paste0("DAG ", i), numNodes(g), numEdges(g), paste(str_replace(extractEdges(amat), "->", "→"), collapse = ", "))
        }
    } else 
    {
        for (i in 1:length(l)) 
        {
            amat <- l[[i]]
            g <- as(amat, "graphNEL")
            if (identical(class(as.numeric(names(l)[i])), "numeric") && !is.na(as.numeric(names(l)[i])))
            {
              res[i, ] <- c(paste0("DAG ", names(l)[i]), numNodes(g), numEdges(g), paste(str_replace(extractEdges(amat), "->", "→"), collapse = ", "))
            } else
            {
              res[i, ] <- c(paste0("", names(l)[i]), numNodes(g), numEdges(g), paste(str_replace(extractEdges(amat), "->", "→"), collapse = ", "))
            }
        }
    }
    res <- as.data.frame(res)
    return(res)
}


boxTitle <- function(title) {
    p(title, style = "padding-right: 5px; display: inline")
}


getMethodsForDataTransform <- function()
{
    return(c("rescale to a specific interval", "log-transform", "z-score", 
             "l2 norm", 
             "pooling normalization and log (for scRNA-seq count data)"))
}


zscore <- function(vec)
{
    if (length(unique(vec)) == 1)
    {
        return(vec)
    } else
    {
        return((vec - mean(vec)) / sd(vec))
    }
    
}


modified_zscore <- function(vec)
{   
    # https://www.statology.org/modified-z-score/
    if (length(unique(vec)) == 1)
    {
        return(vec)
    } else
    {
        mad <- stats::mad(vec, constant = 1)
        return(0.6745 * (vec - median(vec)) / mad)
    }
    
}


l2norm <- function(vec)
{
    return(vec / base::norm(vec, type = "2"))
}


normForSingleCell <- function(mat, base = 2, pseudoCount = 1)
{
    # mat: gene x cell counts matrix

    # Data will be normalized by the pooling normalization method implemented in the scran R package
    # and log2-transformed.

    # ref: https://github.com/Winnie09/imputationBenchmark/blob/master/data/code/process/02_make10xcellline.R

    sce <- SingleCellExperiment::SingleCellExperiment(list(counts=mat))
    if (ncol(mat) < 21){
        sce <- scran::computeSumFactors(sce, BPPARAM=BiocParallel::MulticoreParam(workers = 5), sizes=c(5,10,15,20))
    } else {
        sce <- scran::computeSumFactors(sce, BPPARAM=BiocParallel::MulticoreParam(workers = 5))  
    }
    sf <- SingleCellExperiment::sizeFactors(sce)
    normmat <- sweep(mat, 2, sf, '/')
    normmat <- log(normmat + pseudoCount, base = base)
    return(normmat)
}


calPercentageExpressedCells <- function(vec)
{
    if (any(is.na(vec)))
    {
        return(NA)
    }

    if (min(vec) >= 0)
    {
        return(round(sum(vec > 0) / length(vec), 4))
    } else
    {
        stop("Have negative values!")
    }
  
}


runFeatureSelection <- function(fs_obj,
                                filteredData, 
                                responseData, 
                                method,
                                selectedFeatureNumber = 50) 
{
  if (method == "Random Forest"){
      fs_results <- fs_obj$RandomForest(filteredData, responseData, 
                                        as.integer(selectedFeatureNumber))
  } else if (method == "Extra Trees"){
      fs_results <- fs_obj$ExtraTrees(filteredData, responseData, 
                                      as.integer(selectedFeatureNumber))
  } else if (method == "XGBoost"){
      fs_results <- fs_obj$XGBoost(filteredData, responseData, 
                                   as.integer(selectedFeatureNumber))
  } else if (method == "BAHSIC"){
      fs_results <- fs_obj$BAHSIC(filteredData, responseData, 
                                  as.integer(selectedFeatureNumber))
  } else if (method == "SHS"){
      fs_results <- fs_obj$SHS(filteredData, responseData, 
                               as.integer(selectedFeatureNumber))
  } else if (method == "HSICLasso"){
      fs_results <- fs_obj$HSICLassoRegression(filteredData, responseData, 
                                               as.integer(selectedFeatureNumber))
  } else if (method == "Lasso"){
      fs_results <- fs_obj$LassoRegression(filteredData, responseData, 
                                           as.integer(selectedFeatureNumber))
  } else if (method == "Ridge"){
      fs_results <- fs_obj$RidgeRegression(filteredData, responseData, 
                                           as.integer(selectedFeatureNumber))
  } else if (method == "Elastic Net"){
      fs_results <- fs_obj$ElasticNetRegression(filteredData, responseData, 
                                                as.integer(selectedFeatureNumber))
  }
  
  return(fs_results)
}


makeVisNetwork <- function(amat, 
                           main = NULL, 
                           submain = NULL, 
                           navigationButtons = TRUE)
{
    # submain = list(text = "Custom subtitle", style = "font-family:Comic Sans MS;color:#ff0000;font-size:15px;text-align:center;")

    nodes <- data.frame(id = 1:ncol(amat), 
                    # add labels on nodes
                    label = colnames(amat), 
                    
                    title = paste0("<p>", colnames(amat), "</p>"), 
                    stringsAsFactors = TRUE
    )

    edges <- data.frame(from = unname(which(amat == 1, arr.ind = TRUE)[, "row"]), 
                        to = unname(which(amat == 1, arr.ind = TRUE)[, "col"])
    )

    visnetworkObj <- visNetwork(nodes, edges, 
                                main = main, 
                                submain = submain, 
                                # footer = "Fig.1 minimal example",
                                width = "100%") %>%
                       visNodes(shape = "box", 
                                color = list(background = "lightblue", 
                                             border = "darkblue", 
                                             highlight = "yellow")) %>%
                       visEdges(arrows = list(to = list(enabled = TRUE)), 
                                color = list(color = "lightblue", highlight = "red")) %>%
                       visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T), 
                                  nodesIdSelection = TRUE)

    # if (ggm::isAcyclic(amat))
    # {
    #     visnetworkObj <- visnetworkObj %>% visHierarchicalLayout()
    # }

    visnetworkObj <- visnetworkObj %>% visInteraction(navigationButtons = navigationButtons)

    return(visnetworkObj)
}


getAvailableDataset <- function(obj)
{
    vec <- c()
    if (!is.null(obj$importedData))
    {
        vec <- c(vec, "imported dataset")
    }
    if (!is.null(obj$subsetData))
    {
        vec <- c(vec, "subset dataset")
    }
    if (!is.null(obj$transformedData))
    {
        vec <- c(vec, "transformed dataset")
    }
    return(vec)
}


getGenesOfKeggAmat <- function(amat)
{
  return(unique(Reduce(c, 
                       sapply(rownames(amat), function(x){strsplit(x, split = " ")[[1]]}))
         )
  )
}


findIndexInKeggAmat <- function(gene, keggAmat)
{
  for (i in 1:ncol(keggAmat))
  {
    node <- colnames(keggAmat)[i]
    node <- strsplit(node, " ")[[1]]
    if (gene %in% node)
    {
      return(i)
    }
  }
  return(NA)
}


sorted_by_group <- function(df, group_name, specific_order)
{
  new_order <- c()
  for (i in 1:length(specific_order))
  {
    new_order <- c(new_order, which(df[, group_name] == specific_order[i]))
  }
  sortedData <- df[new_order, ]
  rownames(sortedData) <- 1:nrow(df)
  return(sortedData)
}


changeStyleExpressedCellsForReactable <- function()
{
    jscode <- "function(cellInfo) {
                 let value = '' 
                 let pct = '';
                 
                 if (!cellInfo || isNaN(cellInfo.value)) {
                  value = 'NA'
                  value = value.padStart(5)
                  pct = '0%';
                 } else {  

                   if (cellInfo.value == 1)
                   {
                      // Format as percentage
                      pct = (cellInfo.value * 100) + '%'
                   } else 
                   {
                     // Format as percentage
                     pct = (cellInfo.value * 100).toFixed(1) + '%';
                   }
                     

                   // Pad single-digit numbers
                   value = pct.padStart(5)
                 }
                  
                 // Render bar chart
                 return (
                   '<div class=\"bar-cell\">' +
                     '<span class=\"number\">' + value + '</span>' +
                     '<div class=\"bar-chart\" style=\"background-color: #e1e1e1\">' +
                       '<div class=\"bar\" style=\"width: ' + pct + '; background-color: #fc5185\"></div>' +
                     '</div>' +
                   '</div>'
                 )
              }"
    return(jscode)
}


cor.mtest2 <- function(mat, ...) 
{
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


dccRelationMatrix <- function(df, digits = 3)
{
  # df: row is sample, column is feature.
  
  if (nrow(df) < 100)
  {
    k <- 10
  } else
  {
    k <- 100
  }
  
  dccRes <- pRes <- matrix(0, nrow = ncol(df), ncol = ncol(df))
  for (i in 1:(ncol(df)-1))
  {
    for (j in (i+1):ncol(df))
    {
      dccRes[i, j] <- dccRes[j, i] <- round(as.numeric(energy::dcor.test(df[, i], df[, j])$statistic), digits)
      pRes[i, j] <- pRes[j, i] <- round(kpcalg::dcov.gamma(df[, i], df[, j], numCol = k)$p.value, digits)
    }
  }
  
  # compute dcc on principle diagonal
  for (i in 1:ncol(df))
  {
    res <- dcor.test(df[, i], df[, i])
    dccRes[i, i] <- round(as.numeric(energy::dcor.test(df[, i], df[, i])$statistic), digits)
    pRes[i, i] <- round(kpcalg::dcov.gamma(df[, i], df[, i], numCol = k)$p.value, digits)
  }
  
  rownames(dccRes) <- colnames(dccRes) <- rownames(pRes) <- colnames(pRes) <- colnames(df)
  
  return(list(dccMatrix = dccRes, pMatrix = pRes))
}


mi <- function(x, y)
{
    # 2019 Nature Methods, association benchmarking
    d <- data.frame(x=x, y=y)
    return(WGCNA::mutualInfoAdjacency(d)$AdjacencyUniversalVersion1[1,2])
}


mi.test <- function(x, y, n_perm = 100)
{
    estimate <- mi(x, y)
        
    # permutation test
    vec <- c()
    for (i in 1:n_perm)
    {
        y_shuffle <- sample(y, length(y))
        vec[i] <- mi(x, y_shuffle)
    }

    p <- length(vec[abs(vec) >= abs(estimate)]) / (n_perm + 1) + runif(1, 1e-6, 1) * 1e-8

    return(list(estimate = estimate, pvalue = p))
}


miRelationMatrix <- function(df, n_perm = 100, digits = 3)
{
  # df: row is sample, column is feature.
  
  miRes <- pRes <- matrix(0, nrow = ncol(df), ncol = ncol(df))
  for (i in 1:(ncol(df)-1))
  {
    for (j in (i+1):ncol(df))
    {
      res <- mi.test(df[, i], df[, j], n_perm)
      miRes[i, j] <- miRes[j, i] <- round(res$estimate, digits)
      pRes[i, j] <- pRes[j, i] <- round(res$pvalue, digits)
    }
  }
  
  # compute dcc on principle diagonal
  for (i in 1:ncol(df))
  {
    res <- mi.test(df[, i], df[, i], n_perm)
    miRes[i, i] <- round(res$estimate, digits)
    pRes[i, i] <- round(res$pvalue, digits)
  }
  
  rownames(miRes) <- colnames(miRes) <- rownames(pRes) <- colnames(pRes) <- colnames(df)
  
  return(list(miMatrix = miRes, pMatrix = pRes))
}


similarity.test <- function(x, y, method, n_perm = 100)
{
    # method:
    #  - pearson
    #  - spearman
    #  - kendall
    #  - MI
    #  - dcor

    if (method == "distance correlation")
    {
        estimate <- as.numeric(energy::dcor.test(x, y)$statistic)

        # permutation test
        # vec <- c()
        # for (i in 1:n_perm)
        # {
        #     y_shuffle <- sample(y, length(y))
        #     vec[i] <- as.numeric(energy::dcor.test(x, y_shuffle)$statistic)
        # }

        # p <- length(vec[abs(vec) >= abs(estimate)]) / (n_perm + 1) + runif(1, 1e-6, 1) * 1e-8

        # calculate p by dcov.gamma
        if (length(x) < 100)
        {
            k <- 10
        } else
        {
            k <- 100
        }
        p <- kpcalg::dcov.gamma(x, y, numCol = k)$p.value
        
    } else if (method == "mutual information")
    {
        res <- mi.test(x, y, n_perm = n_perm)
        estimate <- res$estimate
        p <- res$pvalue
    } else
    {
        tmp_result <- cor.test(x, y, method = method)
        estimate <- tmp_result$estimate
        p <- tmp_result$p.value
    }

    return(list(estimate = estimate, pvalue = p, method = method))
}


getSimilarityGTEx <- function(gene1, gene2, group_data, similarity.method)
{
    # oneData <- split(data, data$primary.disease.or.tissue)
    do.call(rbind,lapply(group_data, function(x){
        if (nrow(x) >= 10)
        {
            dd <- similarity.test(as.numeric(x[, gene1]), as.numeric(x[, gene2]), method=similarity.method)
            data.frame(similarity=dd$estimate, pvalue=dd$pvalue )
        }
        
    }))
}


getSimilarityTCGA <- function(gene1, gene2, group_data, similarity.method)
{
    # oneData <- split(data, data$primary.disease.or.tissue)
    do.call(rbind,lapply(group_data, function(x){
        if (nrow(x) >= 10)
        {
            dd  <- similarity.test(as.numeric(x[, gene1]), as.numeric(x[, gene2]), method=similarity.method)
            data.frame(similarity=dd$estimate, pvalue=dd$pvalue)
        }
    }))
}


getColorName <- function(category)
{
    return(rownames(RColorBrewer::brewer.pal.info %>% dplyr::filter(category == category)))
}


plotlyGenePairsCCLE <- function(gene1, gene2, similarity.method)
{
    if (similarity.method == "dcor")
    {
        method <- "Distance correlation"
    } else if (similarity.method == "MI")
    {
        method <- "Mutual information"
    } else if (similarity.method == "pearson")
    {
        method <- "Pearson"
    } else if (similarity.method == "spearman")
    {
        method <- "Spearman"
    } else if (similarity.method == "kendall")
    {
        method <- "Kendall"
    }

    res <- similarity.test(CCLE_data[, gene1], CCLE_data[, gene2], similarity.method, n_perm = 100)
    df <- data.frame(x = CCLE_data[, gene1], y = CCLE_data[, gene2])
    colnames(df) <- c(gene1, gene2)

    g <- eval(parse(text = paste0("ggplot(df, aes(x = ", gene1, ", y = ", gene2, "))"))) +
           geom_point(size = 1, fill = myColors[5], shape = 21, colour = myColors[5], stroke = 1) +
           geom_smooth(method = 'lm', colour = "black") + 
           ggtitle(paste0("Corr (", method, ") = ", round(res$estimate, 3), "; P-value = ", round(res$pvalue, 3))) + 
           theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), 
                 axis.line = element_line(color="black"), axis.line.x = element_line(color="black")) +
           theme_bw()
    g_plotly <- ggplotly(g)
    g_plotly <- g_plotly %>%
                  layout(title = paste0("Corr (", method, ") = ", round(res$estimate, 3), "; P-value = ", round(res$pvalue, 3)), 
                         font = list(size = plotly_fontsize))
    return(g_plotly)
}


remove_punct <- function(string) 
{
  ### Remove punct in the beginning and the end of a string.
  while (grepl("[[:punct:]]", substr(string, 1, 1))) 
  {
    string <- sub("^.", "", string)
  }

  while (grepl("[[:punct:]]", substr(string, nchar(string), nchar(string)))) 
  {
    string <- sub(".$", "", string)
  }

  # remove whitespace
  string <- trimws(string, which = c("both"))

  return(string)
}


ggcorrplot2 <- function (corr, 
                         method = c("square", "circle"), 
                         type = c("full", "lower", "upper"), 
                         ggtheme = ggplot2::theme_minimal, 
                         title = "", 
                         show.legend = TRUE, 
                         legend.title = "Corr", 
                         show.diag = FALSE, 
                         colors = c("blue", "white", "red"), ## TODO
                         cor.limit = c(-1, 1), ## TODO
                         outline.color = "gray", 
                         hc.order = FALSE, 
                         hc.method = "complete", 
                         lab = FALSE, 
                         lab_col = "black", 
                         lab_size = 4, 
                         p.mat = NULL, 
                         sig.level = 0.05, 
                         insig = c("pch", "blank"), 
                         pch = 4, 
                         pch.col = "black", 
                         pch.cex = 5, 
                         tl.cex = 12, 
                         tl.col = "black", 
                         tl.srt = 45, 
                         digits = 2) 
{
    type <- match.arg(type)
    method <- match.arg(method)
    insig <- match.arg(insig)
    if (inherits(corr, "cor_mat")) {
        cor.mat <- corr
        corr <- .tibble_to_matrix(cor.mat)
        p.mat <- .tibble_to_matrix(attr(cor.mat, "pvalue"))
    }
    if (!is.matrix(corr) & !is.data.frame(corr)) {
        stop("Need a matrix or data frame!")
    }
    corr <- as.matrix(corr)
    corr <- base::round(x = corr, digits = digits)
    if (hc.order) {
        ord <- ggcorrplot:::.hc_cormat_order(corr)
        corr <- corr[ord, ord]
        if (!is.null(p.mat)) {
            p.mat <- p.mat[ord, ord]
            p.mat <- base::round(x = p.mat, digits = digits)
        }
    }
    if (type == "lower") {
        corr <- .get_lower_tri(corr, show.diag)
        p.mat <- .get_lower_tri(p.mat, show.diag)
    }
    else if (type == "upper") {
        corr <- .get_upper_tri(corr, show.diag)
        p.mat <- .get_upper_tri(p.mat, show.diag)
    }
    corr <- reshape2::melt(corr, na.rm = TRUE)
    colnames(corr) <- c("Var1", "Var2", "value")
    corr$pvalue <- rep(NA, nrow(corr))
    corr$signif <- rep(NA, nrow(corr))
    if (!is.null(p.mat)) {
        p.mat <- reshape2::melt(p.mat, na.rm = TRUE)
        corr$coef <- corr$value
        corr$pvalue <- p.mat$value
        corr$signif <- as.numeric(p.mat$value <= sig.level)
        p.mat <- subset(p.mat, p.mat$value > sig.level)
        if (insig == "blank") {
            corr$value <- corr$value * corr$signif
        }
    }
    corr$abs_corr <- abs(corr$value) * 10
    p <- ggplot2::ggplot(data = corr, mapping = ggplot2::aes_string(x = "Var1", 
        y = "Var2", fill = "value"))
    if (method == "square") {
        p <- p + ggplot2::geom_tile(color = outline.color)
    }
    else if (method == "circle") {
        p <- p + ggplot2::geom_point(color = outline.color, shape = 21, 
            ggplot2::aes_string(size = "abs_corr")) + ggplot2::scale_size(range = c(4, 
            10)) + ggplot2::guides(size = FALSE)
    }

    ## Modify here!
    if (identical(cor.limit, c(-1, 1)))
    {
        p <- p + ggplot2::scale_fill_gradient2(low = colors[1], 
                                               high = colors[3], 
                                               mid = colors[2], 
                                               midpoint = 0, 
                                               limit = cor.limit, 
                                               space = "Lab", 
                                               name = legend.title)
    } else if (identical(cor.limit, c(0, 1)))
    {
        # colours = viridis::viridis(256, option = "D")
        # colours = RColorBrewer::brewer.pal(10, "YlOrRd")
        p <- p + ggplot2::scale_fill_gradientn(colours = viridis::viridis(256, option = "D"), 
                                               limit = cor.limit, 
                                               space = "Lab", 
                                               name = legend.title)
    }
    ## End modify!
    
    if (class(ggtheme)[[1]] == "function") {
        p <- p + ggtheme()
    }
    else if (class(ggtheme)[[1]] == "theme") {
        p <- p + ggtheme
    }
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = tl.srt, 
        vjust = 1, size = tl.cex, hjust = 1), axis.text.y = ggplot2::element_text(size = tl.cex)) + 
        ggplot2::coord_fixed()
    label <- round(x = corr[, "value"], digits = digits)
    if (!is.null(p.mat) & insig == "blank") {
        ns <- corr$pvalue > sig.level
        if (sum(ns) > 0) 
            label[ns] <- " "
    }
    if (lab) {
        p <- p + ggplot2::geom_text(mapping = ggplot2::aes_string(x = "Var1", 
            y = "Var2"), label = label, color = lab_col, 
            size = lab_size)
    }
    if (!is.null(p.mat) & insig == "pch") {
        p <- p + ggplot2::geom_point(data = p.mat, mapping = ggplot2::aes_string(x = "Var1", 
            y = "Var2"), shape = pch, size = pch.cex, color = pch.col)
    }
    if (title != "") {
        p <- p + ggplot2::ggtitle(title)
    }
    if (!show.legend) {
        p <- p + ggplot2::theme(legend.position = "none")
    }
    p <- p + ggcorrplot:::.no_panel()
    p
}






































# table <- read.table(paste0("./extdata/hg38_gene_ID_name.tsv.gz"), 
#                     sep = "\t", 
#                     header = TRUE, 
#                     stringsAsFactors = FALSE)
# human_url1 <- "https://www.genecards.org/Search/Keyword?queryString="
# human_url2 <- ""
# humanReactable <- reactable(table, 
#                             resizable = FALSE, 
#                             wrap = FALSE, 
#                             bordered = TRUE, 
#                             showSortable = TRUE, 
#                             filterable = TRUE, 
#                             searchable = TRUE, 
#                             highlight = TRUE, 
                            
#                             showPageSizeOptions = TRUE, 
#                             pageSizeOptions = seq(5, 25, 10), 
#                             defaultPageSize = 25, 
#                             paginationType = "jump", 
#                             pagination = TRUE, 
                            
#                             columns = list(
#                               ensembl = colDef(
#                                 cell = function(value) {
#                                   tags$a(href = paste0(human_url1, value, human_url2), target = "_blank", value)
#                                 }
#                               ), 
#                               symbol = colDef(
#                                 cell = function(value) {
#                                   tags$a(href = paste0(human_url1, value, human_url2), target = "_blank", value)
#                                 }
#                               )
#                             )
# )

# table <- read.table(paste0("./extdata/mm10_gene_ID_name.tsv.gz"), 
#                     sep = "\t", 
#                     header = TRUE, 
#                     stringsAsFactors = FALSE)
# res <- clusterProfiler::bitr(table$symbol, 
#                              fromType = "SYMBOL", 
#                              toType = c("ENTREZID"), 
#                              OrgDb = "org.Mm.eg.db")
# newColumn <- data.frame(matrix(NA, nrow = nrow(table), ncol = 1))
# colnames(newColumn) <- "entrez"
# newColumn[match(res$SYMBOL, table$symbol), ] <- res$ENTREZID
# newTable <- cbind(table[, 1:3], newColumn, table[, 4:5])
# write.csv(newTable, file = "./extdata/mm10_gene_ID_name.csv", row.names = FALSE)
# mouse_url1 <- "http://www.informatics.jax.org/quicksearch/summary?queryType=exactPhrase&query="
# mouse_url2 <- "&submit=Quick+Search"
# mouseReactable <- reactable(newTable, 
#                             resizable = FALSE, 
#                             wrap = FALSE, 
#                             bordered = TRUE, 
#                             showSortable = TRUE, 
#                             filterable = TRUE, 
#                             searchable = TRUE, 
#                             highlight = TRUE, 
                            
#                             showPageSizeOptions = TRUE, 
#                             pageSizeOptions = seq(5, 25, 10), 
#                             defaultPageSize = 25, 
#                             paginationType = "jump", 
#                             pagination = TRUE, 
                            
#                             columns = list(
#                               ensembl = colDef(
#                                 cell = function(value) {
#                                   tags$a(href = paste0(mouse_url1, value, mouse_url2), target = "_blank", value)
#                                 }
#                               ), 
#                               symbol = colDef(
#                                 cell = function(value) {
#                                   tags$a(href = paste0(mouse_url1, value, mouse_url2), target = "_blank", value)
#                                 }
#                               ), 
#                               entrez = colDef(
#                                 cell = function(value) {
#                                   tags$a(href = paste0(mouse_url1, value, mouse_url2), target = "_blank", value)
#                                 }
#                               )
#                             )
# )
# save(mouseReactable, file = "./extdata/mm10_geneInformation_reactable.RData")