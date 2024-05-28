suppressPackageStartupMessages(library(BoolNet))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(tictoc))

##############
# LOAD FILES #
##############

bnet_file <- 'M_model_for_PBN.txt'
patient_file <- 'mutation_subject_gene_complete.tsv'

bnet <- BoolNet::loadNetwork(bnet_file)
sample_infor <- read_delim(patient_file, delim = '\t', show_col_types = FALSE)

##################################################
# FUNCTIONS FOR GENERATING PATIENT SPECIFIC PBNS #
##################################################

convert_to_PBM <- function(bnet, genes, dummy_string = 'DUMMY_') {
  # Generate dummy genes
  dummy_genes <- paste0(dummy_string, genes)
  
  # indices
  dummy_indices <- seq_along(dummy_genes) + length(bnet$genes)
  
  # new_genes
  new_genes <- c(bnet$genes, dummy_genes)
  
  # new_fixed
  dummy_fixed <- rep_along(-1, along = dummy_genes)
  names(dummy_fixed) <- dummy_genes
  new_fixed <- c(bnet$fixed, dummy_fixed)
  
  # new_prob_interactions
  dummy_interactions <- lapply(seq_along(dummy_genes), function(x) {
    list(input = dummy_indices[x], func = c(0,1), expression = dummy_genes[x])
  }) %>% set_names(dummy_genes)
  
  new_interactions <- c(bnet$interactions, dummy_interactions)
  
  new_prob_interactions <- new_interactions %>% lapply(function(x) {
    x$probability <- 1.0
    list(x)
  })
  
  new_prob_interactions_for_given_genes <- lapply(seq_along(genes), function(x) {
    gene <- genes[x]
    interaction_org <- new_prob_interactions[[gene]][[1]]
    
    interaction_AND <- interaction_org
    interaction_AND$input <- c(dummy_indices[x], interaction_org$input)
    interaction_AND$func <- c(rep_along(0, along = interaction_org$func), interaction_org$func)
    interaction_AND$expression <- paste0(dummy_genes[x], ' & (', interaction_org$expression, ')')
    interaction_AND$probability <- 0
    
    interaction_OR <- interaction_org
    interaction_OR$input <- c(dummy_indices[x], interaction_org$input)
    interaction_OR$func <- c(interaction_org$func, rep_along(1, along = interaction_org$func))
    interaction_OR$expression <- paste0(dummy_genes[x], ' | (', interaction_org$expression, ')')
    interaction_OR$probability <- 0
    
    list(interaction_org, interaction_AND, interaction_OR)
  }) %>% set_names(genes)
  
  new_prob_interactions[genes] <- new_prob_interactions_for_given_genes
  
  # Generate new BoolNet object
  bnet_new <- list(genes = new_genes, interactions = new_prob_interactions, fixed = new_fixed, dummy_string = dummy_string)
  
  # set class
  class(bnet_new) <- c('ProbabilisticBooleanNetwork', 'BooleanNetworkCollection')
  
  bnet_new
}



generate_patient_PBN <- function(bnet, sample_infor) {
  names(sample_infor)[1] <- 'gene'
  
  if(!all(sample_infor$gene %in% bnet$genes)) {
    err_str <- paste(setdiff(sample_infor$gene, bnet$genes), sep = ',')
    stop(paste0('invalid gene names: ', err_str))
  }
  
  sample_infor <- sample_infor %>% as.data.frame %>% column_to_rownames('gene')
  
  sample_prob_list <- sample_infor %>% lapply(function(x) {
    names(x) <- rownames(.)
    x[!is.na(x) & x != 0]
  })
  
  prob_genes <- sort(unique(names(unlist(unname(sample_prob_list)))))
  
  
  # Convert Deterministic Boolean network into PBN
  m_prob <- convert_to_PBM(bnet = bnet, genes = prob_genes)
  
  # Generate patient-specific BoolNet object
  bnet_list <- lapply(sample_prob_list, function(x) {
    
    m_prob_new <- m_prob
    
    genes <- names(x)
    gene_probs <- x
    gene_interaction_indices <- ifelse(gene_probs < 0, 2, 3)
    gene_interaction_indices_1 <- ifelse(gene_probs < 0, 3, 2)
    
    #
    m_prob_new$interactions[genes] <-
      lapply(seq_along(genes), function(y) {
        new_interactions <- m_prob_new$interactions[[genes[y]]]
        new_interactions[[gene_interaction_indices[y]]]$probability <- abs(gene_probs[y])
        new_interactions[[gene_interaction_indices_1[y]]]$probability <- 0
        new_interactions[[1]]$probability <- 1 - abs(gene_probs[y])
        new_interactions
      }) %>% set_names(genes)
    
    
    # Fix model
    fix_values <- ifelse(gene_probs < 0, 0, 1)
    dummy_genes <- paste0(m_prob_new$dummy_str, genes)
    fixGenes(m_prob_new, fixIndices = dummy_genes, values = fix_values)
  }) %>% set_names(names(sample_prob_list))
  
  bnet_list
}


##########################################
#     GENERATE patient specific PBNS     #
##########################################

m_list <- generate_patient_PBN(bnet = bnet, sample_infor = sample_infor)


############################################################
# Functions for generating patient-specific initial states #
############################################################
generate_random_states <- function(bnet, spec = NULL, n = 100) {
  if(!is.null(spec) & !all(names(spec) %in% bnet$genes)) {
    err_genes <- paste(setdiff(names(spec), bnet$genes), sep = ',')
    stop(paste0('invalid gene names: ', err_genes))
  }
  
  n_var <- length(bnet$genes) - length(spec)
  
  gen_mat <- function(n_rec) {
    ifelse(sapply(seq(n_rec), function(x) { runif(n_var) }) > 0.5, 1, 0) %>% t
  }
  
  n_rec <- n
  rand_mat <- NULL
  while(n_rec > 0) {
    rand_mat <- rbind(rand_mat, gen_mat(n_rec))
    rand_mat <- rand_mat[!duplicated(rand_mat),,drop=FALSE]
    n_rec <- n - nrow(rand_mat)
  }
  
  fix_indices <- which(bnet$genes %in% names(spec))
  rand_indices <- setdiff(seq(length(bnet$genes)), fix_indices)
  
  res_mat <- matrix(0, nrow = length(bnet$genes), ncol = n, dimnames = list(bnet$genes, seq(n)))
  res_mat[rand_indices, ] <- rand_mat
  res_mat[fix_indices, ] <- spec
  
  res_mat
}


generate_patient_inits <- function(bnet, n = 100) {
  genes_ON <- names(bnet$fixed)[bnet$fixed == 1]
  genes_OFF <- names(bnet$fixed)[bnet$fixed == 0]
  genes_dummy <- setdiff(bnet$genes[str_starts(bnet$genes, bnet$dummy_string)], c(genes_ON, genes_OFF))
  
  gene_fixed <- c(rep(1, length(genes_ON)), rep(0, length(genes_OFF)), rep(0, length(genes_dummy)))
  names(gene_fixed) <- c(genes_ON, genes_OFF, genes_dummy)
  gene_fixed
  
  init_mat <- generate_random_states(bnet = bnet, spec = gene_fixed, n = n)
  init_mat
}








##########################################
#             RUN SIMULATION             #
##########################################

# Elapsed time
# - [N_samples: 1, N_INITS: 1e1, N_TIME_TRANSIENT: 1000, N_TIME_STEADY: 300] : 4.5s
# - [N_samples: 1, N_INITS: 1e2, N_TIME_TRANSIENT: 1000, N_TIME_STEADY: 300] : 14.485s
# - [N_samples: 1, N_INITS: 1e3, N_TIME_TRANSIENT: 1000, N_TIME_STEADY: 300] : 101.442s
# - [N_samples: 1, N_INITS: 1e4, N_TIME_TRANSIENT: 1000, N_TIME_STEADY: 300] : 1047s (~18m)
# - [N_samples: 1, N_INITS: 1e4, N_TIME_TRANSIENT: 500, N_TIME_STEADY: 100] : 476s (~8m)
# - [N_samples: 10, N_INITS: 1e4, N_TIME_TRANSIENT: 500, N_TIME_STEADY: 100] : 5075s (~8.4m ~ 1.4h)

N_INITS <- 1e4
N_TIME_TRANSIENT <- 500
N_TIME_STEADY <- 100
N_CORES <- 50
SUBSET_INDICES <- 1:length(m_list)

M_LIST <- m_list[SUBSET_INDICES]

res_save_name <- sprintf('res_obj_Microglia_PBN_code_20230509_%d_%d.RDS', head(SUBSET_INDICES,1), tail(SUBSET_INDICES, 1))

tic()
res_sim_list <- list()
for(i_pat in seq_along(M_LIST)) {
  m_pat <- M_LIST[[i_pat]]
  
  # Generate matrix for initial states
  init_mat <- generate_patient_inits(bnet = m_pat, n = N_INITS)
  
  res_list <- seq(ncol(init_mat)) %>% mclapply(mc.cores = N_CORES, mc.cleanup = TRUE, function(x) {
    # Run transient dynamics
    c_state <- init_mat[, x]
    for(i in seq(N_TIME_TRANSIENT)) {
      c_state <- stateTransition(m_pat, state = c_state, type="probabilistic")
    }
    
    # Run steady state dynamics
    c_state_aggr <- 0
    for(i in seq(N_TIME_STEADY)) {
      c_state <- stateTransition(m_pat, state = c_state, type="probabilistic")
      c_state_aggr <- c_state_aggr + c_state
    }
    averaged_state <- c_state_aggr / N_TIME_STEADY
    
    # Return result
    list(final_state = c_state, averaged_state = averaged_state)
  })
  
  sprintf('%dth loop completed', i_pat)
  
  res_sim_list[[i_pat]] <- res_list
}
toc()

# Set names to list
res_sim_list <- res_sim_list %>% set_names(names(M_LIST))

# single final state
res_sim_fstate_list <- res_sim_list %>% lapply(function(x_pat) {
  x_pat %>% lapply(function(x_init) { x_init$final_state[!str_starts(names(x_init$final_state), 'DUMMY_')] }) %>% do.call(cbind, .)
})
res_sim_fstate_mean_mtx <- res_sim_fstate_list %>% lapply(function(x) {apply(x, 1, mean)}) %>% do.call(cbind, .) %>% as.data.frame

# (pseudo-)steady state
res_sim_astate_list <- res_sim_list %>% lapply(function(x_pat) {
  x_pat %>% lapply(function(x_init) { x_init$averaged_state[!str_starts(names(x_init$averaged_state), 'DUMMY_')] }) %>% do.call(cbind, .)
})
res_sim_astate_mean_mtx <- res_sim_astate_list %>% lapply(function(x) {apply(x, 1, mean)}) %>% do.call(cbind, .) %>% as.data.frame

# Save results
list(M_LIST = M_LIST,
     res_sim_fstate_list = res_sim_fstate_list,
     res_sim_fstate_mean_mtx = res_sim_fstate_mean_mtx,
     res_sim_astate_list = res_sim_astate_list,
     res_sim_astate_mean_mtx = res_sim_astate_mean_mtx) %>%
  saveRDS(res_save_name)

