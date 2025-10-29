rm(list = ls())
library(Seurat)
library(hdWGCNA)
library(diffuStats)
library(tidyverse)
library(igraph)
library(devtools)
# devtools::install_github("cysouw/qlcMatrix")
library(qlcMatrix)
library(ggpubr)
library(data.table)
gc()

###############################################################################################################################################

## Loading all diffusion results
diff_dir <- './Grid_test_0715/'
# Graves_dir <- './Graves_Grid/'

list.res <- list.files(diff_dir, pattern = '.scores') 


tmp <- NULL
anal.data <- NULL


for(i in list.res){
  
  # colnames(tmp) <- 'Final.Heat'
  if(length(grep('hyper', i, value = F)) != 0){
    
    tmp <- readRDS(file = paste0(diff_dir, i)) %>% data.frame()
    
    tmp$Phenotype <- rep('hyper', nrow(tmp))
    
  }else{
    
    if(length(grep('hypo', i, value = F)) != 0){
      
      tmp <- readRDS(file = paste0(diff_dir, i)) %>% data.frame()
      
      tmp$Phenotype <- rep('hypo', nrow(tmp))
      
    }else{
      tmp <- readRDS(file = paste0(diff_dir, i)) %>% data.frame()
      
      tmp$Phenotype <- rep('Graves', nrow(tmp))
      
    }
  }
  
  if(i == list.res[1]){
    anal.data <- tmp
  }else{
    anal.data <- rbind(anal.data, tmp)
  }
  
}

anal.data <- anal.data %>%
  mutate(.,Cell_type = Analysis_name, gene = node_id) %>%
  dplyr::select(., -c('node_id', 'Analysis_name'))


list.cell <- list.files('./Coex_Net_Thyr/', pattern = 'rds') %>%
  gsub(' cells|.rds', '', .) %>%
  gsub(' or ', '.', .)

list.net <- list.files('./Coex_Net_Thyr/', pattern = 'rds',
                       full.names = T)

corr_vec <- c(rep(0, length(list.net)))
names(corr_vec) <- list.cell

anal.data$Cell_type <- anal.data$Cell_type %>%
  gsub('Myeloids', 'Myeloid_cells', .) %>%
  gsub('Fibroblasts', 'Fibroblasts_cells', .) %>%
  gsub('SMCs or Pericytes', 'SMCs.Pericytes_cells', .) %>%
  gsub('T cells or NK cells', 'T.NK_cells', .) %>%
  gsub(' cells', '_cells', .) 


head(anal.data)

rownames(anal.data) <- seq(1:nrow(anal.data))

Module_frame <- c()


###################################### Above this line  ###############################################
gc()
library(future)
options(future.globals.maxSize = 5 * 1024^4)
plan(multisession, workers = 64)

i <- NULL

for(i in 1:length(list.net)){
  
  
  whole.net <- readRDS(list.net[i]) %>%
    ModuleEigengenes() %>%
    ModuleConnectivity(., group_name = list.cell[i]) 
  
  tmp.net <- whole.net %>%
    GetModules() %>%
    dplyr::select(gene_name, module)
  
  assign(paste0(list.cell[i], '_network'), whole.net)
  
  tmp.net$Cell_type <- rep(paste0(list.cell[i], '_cells'), nrow(tmp.net))
  
  colnames(tmp.net) <- c('gene', 'Module', 'Cell_type')
  
  
  Module_frame <- rbind(Module_frame, tmp.net)
  
}

rm(whole.net)

new_anal <- left_join(anal.data, Module_frame, by = c('gene', 'Cell_type'))

load('./All_Thyroid_TWAS.Rdata')
head(hyper_Data)

hyper_Data <- hyper_Data %>%
  dplyr::filter(., Tissue == 'Thyroid') 

hypo_Data <- hypo_Data %>%
  dplyr::filter(., Tissue == 'Thyroid') 

Graves_data <- read_rds(file = './TWAS_res/Graves_res.rds')


head(new_anal)

new_anal$Module %>% is.na %>% sum

new_anal$Cell_type <- new_anal$Cell_type %>%
  gsub('B_cells', 'B Cells', .) %>%
  gsub('Endothelial_cells', 'Endothelial Cells', .) %>%
  gsub('Epithelial_cells', 'Epithelial Cells', .) %>%
  gsub('Fibroblasts_cells', 'Fibroblasts', .) %>%
  gsub('Myeloid_cells', 'Myeloid Cells', .) %>%
  gsub('Proliferating_cells', 'Proliferating Cells', .) %>%
  gsub('SMCs.Pericytes_cells', 'SMCs and Pericytes', .) %>%
  gsub('T.NK_cells', 'T/NK cells', .)


new_anal$TWAS.Z <- rep(0,nrow(new_anal))
new_anal$TWAS.P <- rep(0,nrow(new_anal))

head(new_anal)

for(i in 1:nrow(new_anal)){
  tmp.pheno <- new_anal$Phenotype[i]
  tmp.gene <- new_anal$gene[i]
 
  
  if(tmp.pheno == 'hyper'){
    tmp.TWAS.ind <- hyper_Data[hyper_Data$SYMBOL == tmp.gene, ]
    
  }else{
    if(tmp.pheno == 'hypo'){
      tmp.TWAS.ind <- hypo_Data[hypo_Data$SYMBOL == tmp.gene, ]
    }else{
      tmp.TWAS.ind <- Graves_data[Graves_data$SYMBOL == tmp.gene, ]
      
    }
    
  }
  
  if(nrow(tmp.TWAS.ind) == 0){
    new_anal$TWAS.Z[i] <- 0
    new_anal$TWAS.P[i] <- 1
  }else{
    if(nrow(tmp.TWAS.ind) > 1){
      
      cat(paste0(tmp.gene, ' in ', tmp.pheno, '\n'))
      
      signif.P <- tmp.TWAS.ind$TWAS.P %>%
        min
      
      tmp.Z <- tmp.TWAS.ind[tmp.TWAS.ind$TWAS.P == signif.P, ] %>% unique()
      
      new_anal$TWAS.Z[i] <- tmp.Z$TWAS.Z
      new_anal$TWAS.P[i] <- tmp.Z$TWAS.P
      
    }else{
      new_anal$TWAS.Z[i] <- tmp.TWAS.ind$TWAS.Z
      new_anal$TWAS.P[i] <- tmp.TWAS.ind$TWAS.P
      
    }
  }
}

new_anal <- new_anal %>%
  dplyr::filter(., Module != 'grey')

colnames(new_anal)[3] <- 'Final.Heat'

new_anal_2 <- new_anal %>%
  dplyr::mutate(., delta = Final.Heat - TWAS.Z) %>%
  na.omit

rm(list.net, tmp.net, tmp.Z, tmp.TWAS.ind, hyper_Data, hypo_Data, Module_frame, Graves_data)


## Evaluating diffusion algorithms used in diffuStats package.
## Trying to investigate the algorithm that captured the effects of the Initial.heat and network topology evenly.

list.net <- ls() %>%
  grep('_network', ., value = T)

i <- NULL
Net_stat <- NULL

for(i in list.net){
  tmp <- get(i)
  
  tmp.net <- tmp %>% 
    GetTOM() %>%
    as.matrix() %>%
    graph.adjacency(adjmatrix = ., mode = 'undirected', weighted = TRUE) 
  
  E(tmp.net)$weight <- ifelse(E(tmp.net)$weight > 0.1, E(tmp.net)$weight, 0)
  
  
  tmp.degree <- degree(tmp.net)
  tmp.strength <- strength(tmp.net, weights = edge_attr(tmp.net)$weight)
  tmp.eig <- eigen_centrality(tmp.net, weights = edge_attr(tmp.net)$weight)$vector
  tmp.Hub <- hub_score(tmp.net, weights = edge_attr(tmp.net)$weight)$vector
  tmp.auth <- authority_score(tmp.net, weights = edge_attr(tmp.net)$weight)$vector

  
  res.tmp <- data.frame(gene = names(tmp.degree),
                        Degree = tmp.degree,
                        Strength = tmp.strength, 
                        Eigenvector_centrality = tmp.eig,
                        Hub_score = tmp.Hub,
                        Authority_score = tmp.auth)
  
  res.tmp$Network.type <- rep(i, nrow(res.tmp))
  
  Net_stat <- Net_stat %>% rbind(., res.tmp)
  
  rm(tmp, tmp.net, tmp.degree, tmp.auth, tmp.eig, tmp.Hub, res.tmp)
  
  gc()
}

Net_stat <- Net_stat %>%
  dplyr::mutate(., Cell_type = Network.type)

Net_stat$Cell_type <- Net_stat$Cell_type %>%
  gsub('B_network', 'B Cells', .) %>%
  gsub('Endothelial_network', 'Endothelial Cells', .) %>%
  gsub('Epithelial_network', 'Epithelial Cells', .) %>%
  gsub('Fibroblasts_network', 'Fibroblasts', .) %>%
  gsub('Myeloid_network', 'Myeloid Cells', . ) %>%
  gsub('Proliferating_network', 'Proliferating Cells', .) %>%
  gsub('SMCs.Pericytes_network', 'SMCs and Pericytes', .) %>%
  gsub('T.NK_network', 'T/NK cells', .)

setdiff(unique(Net_stat$gene), unique(new_anal_2$gene)) %>% length

setdiff(unique(new_anal_2$gene), unique(Net_stat$gene)) %>% length




Net_check <- left_join(new_anal_2, Net_stat, by = c('Cell_type', 'gene'), 
                       relationship = 'many-to-many') %>%
  dplyr::filter(., Phenotype == 'Graves') 


Net_check <- Net_check %>%
  dplyr::mutate(., Module_per_Cell = paste0(Module, '_', Cell_type))







Check_effect <- Net_check %>%
  group_by(method, Module_per_Cell, .add = TRUE) %>%
  dplyr::summarise(`Cor(Initial Heat, Strength)` = cor(Initial.Heat, Strength, use = 'pairwise.complete.obs', method = 'pearson'),
                   `Cor(Initial Heat, Eigenvector Centrality)` = cor(Initial.Heat, Eigenvector_centrality, use = 'pairwise.complete.obs', method = 'pearson'),
                   `Cor(Eigenvector Centrality, Strength)` = cor(Eigenvector_centrality, Strength, use = 'pairwise.complete.obs', method = 'pearson'),
                   `Cor(Final Heat, Eigenvector Centrality)` = cor(Final.Heat, Eigenvector_centrality, use = 'pairwise.complete.obs', method = 'pearson'),
                   `Cor(Final Heat, Strength)` = cor(Final.Heat, Strength, use = 'pairwise.complete.obs', method = 'pearson'),
                   `Cor(Final Heat, Initial Heat)` = cor(Final.Heat, Initial.Heat, use = 'pairwise.complete.obs', method = 'pearson'))

Method.dir <- './Method_selection'

if(!dir.exists(Method.dir)){
  dir.create(Method.dir)
}


i <- NULL

for(i in 1:3){
  x_col <- 2+i
  y_col <- 3+i
  
  use.frame <- Check_effect %>%
    dplyr::select(., 1:2, colnames(Check_effect)[x_col], colnames(Check_effect)[y_col])
  
  x_name <- colnames(Check_effect)[x_col]
  y_name <- colnames(Check_effect)[y_col]
  
  colnames(use.frame)[c(3,4)] <- c('X', 'Y')
  
  ggplot(use.frame) +
    geom_jitter(aes(x = X, y = Y,
                    colour = method, size = X*Y)) +
    geom_smooth(aes(x = X, y = Y, colour = method), alpha = 0.6, se = F) +
    theme_bw() +
    xlab(as.character(x_name)) +
    ylab(as.character(y_name))
  
}


long.frame <- Check_effect %>%
  pivot_longer(cols = `Cor(Initial Heat, Strength)`:`Cor(Final Heat, Initial Heat)`)

long.frame <- long.frame %>%
  dplyr::filter(., name %in% c('Cor(Final Heat, Initial Heat)', 'Cor(Final Heat, Eigenvector Centrality)', 'Cor(Final Heat, Strength)'))

ggplot(long.frame) +
  geom_jitter(aes(x = name, y = value, colour = method), alpha = 0.3) +
  geom_boxplot(aes(x = name, y = value, fill = method), alpha = 0.65) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5))+
  xlab('Correlation between network and diffusion statistics') +
  ylab('Correlation')


long.frame_2 <- long.frame %>%
  dplyr::mutate(., name_2 = value^2)

ggplot(long.frame_2) +
  geom_jitter(aes(x = name, y = name_2, colour = method), alpha = 0.3) +
  geom_boxplot(aes(x = name, y = name_2, fill = method), alpha = 0.65) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5))+
  xlab('Correlation between network and diffusion statistics') +
  ylab('Correlation')


## Evaluating the model performance using the correlations were deemly failed....
## Trying to evaluate the best method with known Graves markers.
gc()
library(readxl)
# install.packages('UpSetR')

## Importing results form Barrio-Hernandez et al

PPR_res <- read_xlsx('./Benchmark_opps/Barrio-Hernandez et al Supple.xlsx',
                     sheet = 1)

thyr.vec <- PPR_res$EFO_name %>% grep('Graves', .)

PPR_res <- PPR_res[thyr.vec, ]

PPR_thresh <- PPR_res$page.rank %>%
  quantile(., probs = seq(0, 1, 0.25)) %>%
  .[4]

PPR_signi <- PPR_res %>%
  dplyr::filter(., page.rank > PPR_thresh) %>%
  .$gene %>%
  unique()
  


Target.markers <- fread('./Opentarget_Graves_Marker.tsv') %>%
  dplyr::filter(., symbol %in% Net_check$gene)

marker.thresh <- 0.0

Target.markers <- Target.markers %>%
  dplyr::filter(., globalScore > marker.thresh) %>%
  .$symbol

library(EnsDb.Hsapiens.v86)

ensemble.genes <- genes(EnsDb.Hsapiens.v86)


ent2sym <- ensembldb::select(EnsDb.Hsapiens.v86, 
                              keys= ensemble.genes$gene_id, 
                              keytype = "GENEID", 
                              columns = c("SYMBOL","GENEID", "ENTREZID"))
colnames(ent2sym)[3] <- 'Entrez Gene ID'

OMIM_markers <- fread('./OMIM-Entry-Retrieval.tsv')

OMIM_markers <- OMIM_markers %>%
  dplyr::inner_join(., ent2sym, by = 'Entrez Gene ID', relationship = 'many-to-many') %>%
  na.omit() %>%
  .$SYMBOL %>%
  unique()


Malacard_markers <- fread('./MalaCards - Genes associated with Graves Disease.csv') %>%
  .$Symbol %>%
  unique

All.markers <- Target.markers %>%
  union(., OMIM_markers) %>%
  union(., Malacard_markers) %>%
  unique


Net_check <- Net_check %>%
  dplyr::mutate(., is.known.target = if_else(gene %in% All.markers,
                                             'Known_Marker', 'Random'),
                TWAS.signif = if_else(TWAS.P < 0.05/11337, 'Sig', 'non-sig'))

# ggplot(Net_check) +
#   geom_boxplot(aes(x = method, y = Strength, fill = is.known.target)) +
#   theme_bw()
# 
# 
# ggplot(Net_check) +
#   geom_boxplot(aes(x = method, y = abs(Final.Heat), fill = is.known.target)) +
#   theme_bw()
# 
# ggplot(Net_check) +
#   geom_boxplot(aes(x = method, y = abs(TWAS.Z), fill = is.known.target)) +
#   theme_bw()
# 
# ggplot(Net_check) +
#   geom_boxplot(aes(x = method, y = abs(TWAS.Z), fill = TWAS.signif)) +
#   theme_bw()
# 
# 
# ggplot(subset(Net_check, Final.Heat > 0.9997|Final.Heat < -0.14935)) +
#   geom_boxplot(aes(x = method, y = abs(Final.Heat), fill = is.known.target)) +
#   theme_bw()
# 
# 
# ggplot(Net_check) +
#   geom_boxplot(aes(x = method, y = abs(Initial.Heat), fill = is.known.target)) +
#   theme_bw()
# 
# ggplot(subset(Net_check, abs(delta))) +
#   geom_boxplot(aes(x = method, y = abs(delta), fill = is.known.target)) +
#   theme_bw()

## 071824 Box plot cannot present meaningful results.
## Trying to calculate the overlaps with Fisher's test

sig_tab_dir <- './Sig_gene_per_method'

if(!dir.exists(sig_tab_dir)){
  dir.create(sig_tab_dir)
}


Raw_frame <- Net_check %>%
  dplyr::filter(., method == 'raw') %>%
  dplyr::mutate(., Diffuse.signif = if_else(abs(Final.Heat) > quantile(abs(Final.Heat), probs = seq(0, 1, 0.01))[100],
                'Sig', 'non-sig'),
                Delta.signif = if_else(abs(delta) > quantile(abs(delta), probs = seq(0, 1, 0.01))[100],
                                         'Sig', 'non-sig'))

Raw_TWAS <- Raw_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & TWAS.signif == 'Sig') %>%
  .$gene %>% unique() %>% length()

Raw_Diffuse <- Raw_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & (Diffuse.signif == 'Sig' | Delta.signif == 'Sig')) %>%
  .$gene %>% unique() %>% length()


Ber_S_frame <- Net_check %>%
  dplyr::filter(., method == 'ber_s') %>%
  dplyr::mutate(., Diffuse.signif = if_else(abs(Final.Heat) > quantile(abs(Final.Heat), probs = seq(0, 1, 0.01))[100],
                                            'Sig', 'non-sig'),
                Delta.signif = if_else(abs(delta) > quantile(abs(delta), probs = seq(0, 1, 0.01))[100],
                                       'Sig', 'non-sig'))

Ber_S_TWAS <- Ber_S_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & TWAS.signif == 'Sig') %>%
  .$gene %>% unique() %>% length()

Ber_S_Diffuse <- Ber_S_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & (Diffuse.signif == 'Sig' | Delta.signif == 'Sig')) %>%
  .$gene %>% unique() %>% length()


Ber_P_frame <- Net_check %>%
  dplyr::filter(., method == 'ber_p') %>%
  dplyr::mutate(., Diffuse.signif = if_else(abs(Final.Heat) > quantile(abs(Final.Heat), probs = seq(0, 1, 0.01))[100],
                                            'Sig', 'non-sig'),
                Delta.signif = if_else(abs(delta) > quantile(abs(delta), probs = seq(0, 1, 0.01))[100],
                                       'Sig', 'non-sig'))

Ber_P_TWAS <- Ber_P_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & TWAS.signif == 'Sig') %>%
  .$gene %>% unique() %>% length()

Ber_P_Diffuse <- Ber_P_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & (Diffuse.signif == 'Sig' | Delta.signif == 'Sig')) %>%
  .$gene %>% unique() %>% length()

MC_frame <- Net_check %>%
  dplyr::filter(., method == 'mc') %>%
  dplyr::mutate(., Diffuse.signif = if_else(abs(Final.Heat) > quantile(abs(Final.Heat), probs = seq(0, 1, 0.01))[100],
                                            'Sig', 'non-sig'),
                Delta.signif = if_else(abs(delta) > quantile(abs(delta), probs = seq(0, 1, 0.01))[100],
                                       'Sig', 'non-sig'))

MC_TWAS <- MC_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & TWAS.signif == 'Sig') %>%
  .$gene %>% unique() %>% length()

MC_Diffuse <- MC_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & (Diffuse.signif == 'Sig' | Delta.signif == 'Sig')) %>%
  .$gene %>% unique() %>% length()

Z_frame <- Net_check %>%
  dplyr::filter(., method == 'z') %>%
  dplyr::mutate(., Diffuse.signif = if_else(abs(Final.Heat) > quantile(abs(Final.Heat), probs = seq(0, 1, 0.01))[100],
                                            'Sig', 'non-sig'),
                Delta.signif = if_else(abs(delta) > quantile(abs(delta), probs = seq(0, 1, 0.01))[100],
                                       'Sig', 'non-sig'))

Z_TWAS <- Z_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & TWAS.signif == 'Sig') %>%
  .$gene %>% unique() %>% length()

Z_Diffuse <- Z_frame %>%
  dplyr::filter(., is.known.target == 'Known_Marker' & (Diffuse.signif == 'Sig' | Delta.signif == 'Sig')) %>%
  .$gene %>% unique() %>% length()


signif.tab <- rbind(MC_frame, Ber_P_frame, Z_frame, Ber_S_frame, Raw_frame)

save(signif.tab, file = paste0(sig_tab_dir, '/Dat_with_significance.Rdata'))

# load('./Sig_gene_per_method/Dat_with_significance.Rdata')
New_Methods <- c('NESTA (raw)', 'NESTA (ber_s)', 'NESTA (ber_p)', 'NESTA (mc)', 'NESTA (z)')

Marker.rep <- data.frame(method = c(New_Methods, 'Barrio-Hernandez et al', 'TWAS'),
                         Replication = c(Raw_Diffuse, Ber_S_Diffuse, Ber_P_Diffuse, MC_Diffuse, Z_Diffuse, length(intersect(PPR_signi, All.markers)), Z_TWAS),
                         no.of.genes = c(length(unique(Raw_frame$gene)), length(unique(Ber_S_frame$gene)), length(unique(Ber_P_frame$gene)), 
                                         length(unique(MC_frame$gene)), length(unique(Z_frame$gene)), length(unique(PPR_res$gene)),length(unique(Raw_frame$gene))),
                         no.of.markers = c(rep(length(All.markers), 7)))

Marker.rep$method <- Marker.rep$method %>%
  factor(., levels = c('NESTA (mc)', 'NESTA (ber_s)', 'NESTA (ber_p)', 'NESTA (raw)', 'Barrio-Hernandez et al', 'NESTA (z)', 'TWAS'))

i <- NULL

fish.mat <- c()

for(i in 1:nrow(Marker.rep)){
  tmp.mat <- data.frame(Markers = c(Marker.rep[i,2], Marker.rep[i,4] - Marker.rep[i,2]),
                        non.Markers = c(Marker.rep[i,3] - Marker.rep[i,2], Marker.rep[i, 3] - Marker.rep[i, 4]))
  
  tmp.1 <- fisher.test(tmp.mat, 'greater') %>%
    .$p.value
  
  tmp.2 <- data.frame(method = Marker.rep[i,1],
                      p.value = tmp.1)
  
  fish.mat <- rbind(fish.mat, tmp.2)
  
  
  rm(tmp.mat, tmp.1, tmp.2)
}


ggplot(Marker.rep) +
  geom_text(data = Marker.rep, aes(label = Replication, y = Replication, x = method), 
            vjust = -0.8) +
  geom_col(aes(x = method, y = Replication)) +
  theme_bw() +
  xlab('Methods') +
  ylab('Number of replicated known markers') +
  ylim(c(0,50)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.8, hjust = 0.8))

ggsave(paste0(sig_tab_dir, '/No_of_known_markers.png'), width = 5, height = 5, units = 'in', dpi = 300)

MC_res <- new_anal_2 %>%
  dplyr::filter(method == 'mc')

save(Net_stat, MC_res, file = './073124_data_to_anal.Rdata')




Net_stat_2 <- Net_stat
MC_res_2 <- MC_res

load('./073124_data_to_anal.Rdata')
