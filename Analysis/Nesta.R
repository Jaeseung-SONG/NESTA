# Required packages

if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if(!require("Seurat", quietly = TRUE))
  BiocManager::install("Seurat")

if(!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if(!require("igraph", quietly = TRUE))
  BiocManager("igraph")

if(!require("devtools", quietly = TRUE))
  install.packages("devtools")

if(!require("hdWGCNA", quietly = TRUE))
  BiocManager::install("hdWGCNA")

if(!require("diffuStats", quietly = TRUE))
  BiocManager::install("diffuStats")

if(!require("optparse", quietly = TRUE))
  install.packages("optparse")

if(!require("data.table", quietly = TRUE))
  install.packages("data.table")


# Loading packages

suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library(devtools))
suppressMessages(library(hdWGCNA))
suppressMessages(library(diffuStats))
suppressMessages(library(optparse))
suppressMessages(library(data.table))

# Required option list

option_list = list(
  make_option("--TWAS_res", action = 'store', default = NA, type = 'character',
              help = 'Rdata file for TWAS result'),
  make_option("--Reference_net", action = 'store', default = NA, type = 'character',
              help = 'Reference network file: .rds file or table format (.tsv/.csv)'),
  make_option("--Is_expression_network", action = 'store', default = 'TRUE', type = 'character',
              help = 'Please specify the type of the reference network. \n
              Will try to adjust for cell/tissue-specific expression pattern if it is TRUE'),
  make_option("--Expression_slot", action = 'store', default = NA, type = 'character',
              help = 'If this field is not provided, @assays$SCT@dat slot will be used as default.'),
  make_option("--Diffuse_grid", action = 'store', default = FALSE, type = 'logical',
              help = 'If TRUE, propagation will run for all available methods for quantitative initial heats. \n
              See diffuStats package for the detailed descriptions for each method.'),
  make_option("--Diffuse_method", action = 'store', default = 'ber_p', type = 'character',
              help = "The name of the method passed to the diffuse() function."),
  make_option("--Diffuse_nperm", action = 'store', default = 300, type = 'integer',
              help = "Number of permutations passed to diffuStats mc/ber_p diffusion. Default 300 matches the submitted script intent."),
  make_option("--check_bias", action = 'store', default = TRUE, type= 'logical',
              help = "Checking the presence bias of final heats toward the initial heat using Pearson's R.\n
              Default value is TRUE"),
  make_option("--out_dir", action = 'store', default = './', type= 'character',
              help = "Directory for storing output."),
  make_option("--prefix", action = 'store', default = 'Diffuse_out', type= 'character',
              help = "Prefix for output file"),
  make_option("--TWAS_cutoff", action = 'store', default = 1, type = 'double',
              help = 'TWAS.P value cutoff for selecting the genes for node weight. \n
              The default value is 1 (Using all available TWAS.Z for the initial weight)'),
  make_option("--edge_cutoff", action = 'store', default = 0.1, type = 'double',
              help = 'Edge weight cutoff for selecting the edges. \n
              The default value is 0.2'),
  make_option("--Analysis_name", action = 'store', default = NULL, type = 'character',
              help = 'Please provide the name of the analysis. It can be cell type name or others'),
  make_option("--Initial_weight_mode", action = 'store', default = 'auto', type = 'character',
              help = 'Initial weight mode: auto, twas_only, or nesta_expression_weighted. Default auto preserves existing behavior.'),
  make_option("--Expression_reference_net", action = 'store', default = NA, type = 'character',
              help = 'Optional expression network RDS used only to compute cell-type mean expression when Is_expression_network is FALSE and Initial_weight_mode is nesta_expression_weighted.'),
  make_option("--Kernel_rds", action = 'store', default = NA, type = 'character',
              help = 'Optional precomputed diffuStats kernel RDS for the supplied reference network.'),
  make_option("--Save_kernel_rds", action = 'store', default = NA, type = 'character',
              help = 'Optional path for saving the computed diffuStats kernel RDS.')

)

opt = parse_args(OptionParser(option_list = option_list))

parse_logical_flag <- function(x, option_name){
  if(is.logical(x)){
    if(length(x) == 1 && !is.na(x)) return(x)
  }
  x_chr <- tolower(trimws(as.character(x)))
  if(x_chr %in% c('true', 't', 'yes', 'y', '1')) return(TRUE)
  if(x_chr %in% c('false', 'f', 'no', 'n', '0')) return(FALSE)
  stop('Error: ', option_name, ' must be TRUE/FALSE or YES/NO.')
}

opt$Is_expression_network <- parse_logical_flag(opt$Is_expression_network, '--Is_expression_network')



# Main function

NESTA <- function(TWAS_res, Refer.net, Is.expression, Expression.slot,
                          grid.opt, out.dir, prefix, diffuse.method, Check.redundancy, TWAS.cut, edge_cutoff, Analysis_name,
                          Initial.weight.mode = 'auto', Expression.reference.net = NA, Diffuse.nperm = 300,
                          Kernel.rds = NA, Save.kernel.rds = NA){

  # Checking for input file availability

  if (is.na(TWAS_res) || !file.exists(TWAS_res)) {
    stop("Error: TWAS result file is required and must exist.")
  }
  if (is.na(Refer.net) || !file.exists(Refer.net)) {
    stop("Error: Reference network file is required and must exist.")
  }

  # Creating the output directory (if needed)

  if (!dir.exists(out.dir)) {
    cat(paste0('Creating output directory: ', out.dir, '\n'))
    dir.create(out.dir, recursive = TRUE)
  }

  # Loading reference network to the working environment.

  Is.expression <- parse_logical_flag(Is.expression, '--Is_expression_network')
  Initial.weight.mode <- tolower(Initial.weight.mode)
  if(Initial.weight.mode == 'auto'){
    Initial.weight.mode <- ifelse(Is.expression, 'nesta_expression_weighted', 'twas_only')
  }
  if(!Initial.weight.mode %in% c('twas_only', 'nesta_expression_weighted')){
    stop("Error: Initial_weight_mode must be auto, twas_only, or nesta_expression_weighted.")
  }

  read_mean_expression <- function(expr.path, expr.slot){
    Expr.dat.local <- readRDS(file = expr.path)
    if(is.na(expr.slot)){
      cat('Using default expression slot \n')
      Mean.expr.local <- Expr.dat.local@assays$SCT@data %>%
        rowMeans() %>% data.frame
    }else{
      cat('Using specified expression slot: ', expr.slot, ' \n')
      Mean.expr.local <- Expr.dat.local[[expr.slot]]$data %>%
        rowMeans() %>% data.frame
    }
    Mean.expr.local$SYMBOL <- rownames(Mean.expr.local)
    colnames(Mean.expr.local)[1] <- 'Mean_expression'
    Mean.expr.local
  }

  if(Is.expression){

    cat('Loading Gene Expression-based Network\n')


    Expr.dat <- readRDS(file = Refer.net)

    Ref.net <- Expr.dat %>%
      GetTOM %>%
      graph.adjacency(., mode = 'undirected', weighted = T)


    Mean.expr <- read_mean_expression(Refer.net, Expression.slot)

  }else{


    cat('Loading Topology-only Network\n')

    file.arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
    if(length(file.arg) > 0){
      script.dir <- dirname(normalizePath(sub("^--file=", "", file.arg[1])))
    }else{
      script.dir <- "/home/js/NESTA/Analysis"
    }
    string.helper <- file.path(script.dir, 'string_ppi_utils.R')
    if(!file.exists(string.helper)){
      string.helper <- '/home/js/NESTA/Analysis/string_ppi_utils.R'
    }
    if(file.exists(string.helper)){
      source(string.helper)
      Ref.net <- read_ppi_graph(Refer.net, edge_cutoff = edge_cutoff)
    }else{
      Ref.net <- fread(Refer.net) %>%
        graph.adjacency(., mode = 'undirected', weighted = T)
    }

    if(Initial.weight.mode == 'nesta_expression_weighted'){
      if(is.na(Expression.reference.net) || !file.exists(Expression.reference.net)){
        stop("Error: Expression_reference_net is required for nesta_expression_weighted mode with Is_expression_network FALSE.")
      }
      Mean.expr <- read_mean_expression(Expression.reference.net, Expression.slot)
    }

  }

  E(Ref.net)$weight <- ifelse(E(Ref.net)$weight > edge_cutoff, E(Ref.net)$weight, 0)

  # Generating the initial score vector

  cat('Matching available genes between TWAS and reference network.\n')

  Net.gene.list <- V(Ref.net) %>%
    names()

  tmp.score <- data.frame(SYMBOL = Net.gene.list,
                          weight = rep(0, length(Net.gene.list)))

  test.cutoff <- TWAS.cut

  TWAS.res <- readRDS(file = TWAS_res) %>%
    dplyr::filter(!is.na(SYMBOL),
                  SYMBOL != '',
                  is.finite(TWAS.Z),
                  is.finite(TWAS.P),
                  TWAS.P < test.cutoff) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::slice_min(order_by = TWAS.P, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(SYMBOL, TWAS.Z)

  tmp.score <- tmp.score %>%
    dplyr::left_join(TWAS.res, by = 'SYMBOL') %>%
    dplyr::mutate(weight = if_else(is.na(TWAS.Z), 0, TWAS.Z)) %>%
    dplyr::select(SYMBOL, weight)

  cat('Reference network nodes: ', nrow(tmp.score), '\n')
  cat('Nodes with TWAS-derived nonzero initial weight: ', sum(tmp.score$weight != 0), '\n')
  cat('Zero-imputed network nodes without nonzero TWAS weight: ', sum(tmp.score$weight == 0), '\n')

  if(Initial.weight.mode == 'nesta_expression_weighted'){

    cat('Generating cell type-specific initial weight vector.\n')


    Init.score <- tmp.score %>%
      dplyr::left_join(., Mean.expr, by = 'SYMBOL') %>%
      dplyr::mutate(Mean_expression = if_else(is.na(Mean_expression), 0, Mean_expression))

    if(sd(Init.score$Mean_expression) > 0 & sd(Init.score$weight) > 0){
      Init.score$comb.weight <- (Init.score$Mean_expression/sd(Init.score$Mean_expression)) *
        (Init.score$weight/sd(Init.score$weight))
    }else{
      Init.score$comb.weight <- 0
    }

  }else{

    cat('Generating initial weight vector.\n')

    Init.score <- tmp.score %>%
      mutate(., comb.weight = weight)
  }


  # Diffusing the scores across the provided network.

  Init.score <- Init.score %>%
    dplyr::mutate(., Pos.weight = if_else(comb.weight >= 0, comb.weight, 0),
                  Neg.weight = if_else(comb.weight <= 0, abs(comb.weight), 0))

  score.vec.Pos <- Init.score$Pos.weight
  names(score.vec.Pos) <- Init.score$SYMBOL
  score.vec.Neg <- Init.score$Neg.weight
  names(score.vec.Neg) <- Init.score$SYMBOL

  if(!is.na(Kernel.rds) && file.exists(Kernel.rds)){
    cat('Loading precomputed diffusion kernel: ', Kernel.rds, '\n')
    Diffuse.kernel <- readRDS(Kernel.rds)
  }else{
    cat('Computing diffusion kernel once for this network.\n')
    Diffuse.kernel <- regularisedLaplacianKernel(graph = Ref.net)
    if(!is.na(Save.kernel.rds)){
      cat('Saving diffusion kernel: ', Save.kernel.rds, '\n')
      saveRDS(Diffuse.kernel, file = Save.kernel.rds)
    }
  }

  if(grid.opt){

    cat('Performing network-wise diffusion with multiple methods. \n')

    use.methods <- c('raw', 'ber_s', 'ber_p', 'mc', 'z')


    cat('Diffusing for positive weights \n')


    F.score.Pos <- diffuse_grid(graph = Ref.net,
                            scores = score.vec.Pos,
                            grid_param = expand.grid(method = use.methods),
                            n.perm = Diffuse.nperm,
                            K = Diffuse.kernel,
                            seed = 9703)

    F.score.Pos <- F.score.Pos %>%
      group_by(method) %>%
      mutate(sd_node_score = sd(node_score)) %>%
      ungroup() %>%
      mutate(norm_node_score = node_score/sd_node_score) %>%
      dplyr::select(., -sd_node_score)


    cat('Done \n')


    cat('Diffusing for negative weights \n')


    F.score.Neg <- diffuse_grid(graph = Ref.net,
                            scores = score.vec.Neg,
                            grid_param = expand.grid(method = use.methods),
                            n.perm = Diffuse.nperm,
                            K = Diffuse.kernel,
                            seed = 9703)

    F.score.Neg <- F.score.Neg %>%
      group_by(method) %>%
      mutate(sd_node_score = sd(node_score)) %>%
      ungroup() %>%
      mutate(norm_node_score = node_score/sd_node_score) %>%
      dplyr::select(., -sd_node_score)

    cat('Done \n')


    cat('Computing final scores \n')


    F.score <- data.frame(method = F.score.Pos$method,
                          node_id = F.score.Pos$node_id,
                          Pre_score = F.score.Pos$norm_node_score - F.score.Neg$norm_node_score)


    F.score <- F.score %>%
      group_by(method) %>%
      mutate(sd_Pre_score = sd(Pre_score)) %>%
      ungroup() %>%
      mutate(F.score = Pre_score/sd_Pre_score) %>%
      dplyr::select(., -sd_Pre_score)

    F.score$Initial.Heat <- Init.score$comb.weight %>%
      rep(., times = length(use.methods))



  }else{

    cat('Performing network-wise diffusion with ', diffuse.method,' method. \n')


    cat('Diffusing for positive weights \n')

    F.score.Pos <- diffuse(graph = Ref.net,
                            scores = score.vec.Pos,
                            method = diffuse.method,
                            n.perm = Diffuse.nperm,
                            K = Diffuse.kernel,
                            seed = 9703)

    cat('Done \n')


    cat('Diffusing for negative weights \n')

    F.score.Neg <- diffuse(graph = Ref.net,
                           scores = score.vec.Neg,
                           method = diffuse.method,
                           n.perm = Diffuse.nperm,
                           K = Diffuse.kernel,
                           seed = 9703)

    cat('Done \n')


    cat('Computing final scores \n')


    F.score <- F.score.Pos/sd(F.score.Pos) - F.score.Neg/sd(F.score.Neg)

    F.score <- F.score/sd(F.score)


    F.score <- F.score %>%
      data.frame(Final.Heat = .)

    F.score$SYMBOL <- rownames(F.score)

    F.score <- F.score %>%
      dplyr::select(SYMBOL, Final.Heat)

    F.score <- Init.score %>%
      inner_join(F.score, by = 'SYMBOL')

    F.score$Initial.Heat <- Init.score$comb.weight

  }


  F.score$Analysis_name <- rep(Analysis_name, nrow(F.score))


  write_rds(F.score, file = paste0(out.dir, prefix, '_scores.rds'))


  if(Check.redundancy){
    if(grid.opt){

      corr_vec <- c()

      for(j in use.methods){

        tmp.score <- F.score %>%
          dplyr::filter(., method == j)

        corr.val <- cor(x = tmp.score$F.score, y = Init.score$comb.weight)

        corr_vec <- corr_vec %>% append(., corr.val)

      }


      corr_res <- data.frame(method = use.methods,
                             cor = corr_vec)


      for(k in 1:nrow(corr_res)){

        cat('Correlation for ', corr_res[k,1], ': ', corr_res[k,2], '\n')

      }

      save(corr_res, file = paste0(out.dir, prefix, '_corr.rds'))


    }else{
      corr_res <- cor(x = F.score$Final.Heat, y = Init.score$comb.weight)

      cat('Correlation: ', corr_res, '\n')

    }


  }else{
    cat('Skip for checking redundancies')
  }


}


NESTA(opt$TWAS_res, opt$Reference_net, opt$Is_expression_network,
              opt$Expression_slot,  opt$Diffuse_grid, opt$out_dir,
              opt$prefix, opt$Diffuse_method, opt$check_bias,
              opt$TWAS_cutoff, opt$edge_cutoff, opt$Analysis_name,
              opt$Initial_weight_mode, opt$Expression_reference_net, opt$Diffuse_nperm,
              opt$Kernel_rds, opt$Save_kernel_rds)
