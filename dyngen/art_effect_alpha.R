# Check the controllability of dyngen + GENIE3 system    # 2021/7/24
library(tidyverse)
library(dyngen)
## used in change kinetics
library(assertthat)
library(rlang)
library(stats)
library(GillespieSSA2)
library(dynutils)

dataset_list = list('linear_simple','cycle_simple','bifurcating','converging')

alphas = 10**seq(-2,1,1)
n_para = length(alphas)
n_net = 10
dir = 'dyngen/check_controllability_change_para_error_bar/'
dir.create(dir)
dir.create(paste(dir,"evaluation",sep = ""))
dir.create(paste(dir,"Network_inference",sep = ""))


# control backbone
# data_num = 1
for(data_num in 1:4){
  for(net_id in 1:n_net){
    for(para_id in 1:n_para){
      set.seed(0) # set backbone
      backbone = list_backbones()[[dataset_list[[data_num]]]]()
      # other parameters
      num_cells = 100
      num_simu = 9
      num_tfs = nrow(backbone$module_info)
      TF_prop = 0.1
      num_genes = num_tfs/TF_prop
      num_targets = num_genes - num_tfs
      # change kinetic paramters
      if(T){
        # change splicing time: 2h -> 10/60 h
        splice_time = 10/60
        # change unit
        # alpha_min = 0.1*60
        # alpha_max = 1*60
        alpha = alphas[para_id]*60
        
        kinetics_change1 <- function() {
          # satisfy r cmd check
          transcription_rate <- translation_rate <- mrna_halflife <- protein_halflife <- 
            independence <- splicing_rate <- effect <- strength <- hill <- NULL
          
          sampler_tfs <-
            function(feature_info, feature_network, cache_dir = NULL, verbose = FALSE) {
              feature_info %>% mutate(
                transcription_rate = transcription_rate %|% runif(n(), 10, 20),
                translation_rate = translation_rate %|% runif(n(), 100, 150),
                mrna_halflife = mrna_halflife %|% runif(n(), 1, 10),
                protein_halflife = protein_halflife %|% runif(n(), 1, 10),
                independence = independence %|% 1,
                splicing_rate = splicing_rate %|% (log(2) / splice_time)
              )
            }
          
          sampler_nontfs <-
            function(feature_info, feature_network, cache_dir = NULL, verbose = FALSE) {
              feature_info %>% mutate(
                transcription_rate = transcription_rate %|% runif(n(), alpha, alpha),
                translation_rate = translation_rate %|% runif(n(), 100, 150),
                mrna_halflife = mrna_halflife %|% runif(n(), 1, 10),
                protein_halflife = protein_halflife %|% runif(n(), 1, 10),
                independence = independence %|% 1,
                splicing_rate = splicing_rate %|% (log(2) / splice_time)
              )
            }
          
          sampler_interactions <-
            function(feature_info, feature_network, cache_dir = NULL, verbose = FALSE) {
              feature_network %>% 
                mutate(
                  effect = effect %|% sample(c(-1L, 1L), n(), replace = TRUE, prob = c(.25, .75)),
                  strength = strength %|% 10 ^ runif(n(), log10(1), log10(100)),
                  hill = hill %|% rnorm_bounded(n(), 2, 2, min = 1, max = 10)
                )
            }
          
          lst(sampler_tfs, sampler_nontfs, sampler_interactions)
        }
        environment(kinetics_change1) = environment(kinetics_default)
      }
      
      set.seed(net_id)
      config <- 
        initialise_model(
          num_cells = num_cells,
          num_tfs = num_tfs,
          num_targets = num_targets,
          num_hks = 0,
          backbone = backbone,
          verbose = F,
          download_cache_dir = "./.cache/dyngen",
          num_cores = 20,
          simulation_params = simulation_default(experiment_params = bind_rows(simulation_type_wild_type(num_simulations = num_simu, seed = seq(1,num_simu,1)),
                                                                               simulation_type_knockdown(num_simulations = 0)),
                                                 census_interval = 1
          ),
          kinetics_params  = kinetics_change1()
        )
      # dyngen simulation
      model = config
      model <- generate_tf_network(model)
      model <- generate_feature_network(model)   ### 2 min for download
      # change kinetics calculation: incorporate nonTF parameters in function '.kinetics_generate_gene_kinetics'
      if(T){
        generate_kinetics1 <- function(model) {
          # satisfy r cmd check
          burn <- mol_premrna <- mol_mrna <- mol_protein <- val <- NULL
          
          assert_that(
            !is.null(model$feature_info),
            !is.null(model$feature_network)
          )
          
          # generate kinetics params
          model <- .kinetics_generate_gene_kinetics(model)
          
          # generate formulae
          formulae <- .kinetics_generate_formulae(model)
          
          # create variables
          fid <- model$feature_info$feature_id
          model$feature_info$mol_premrna <- paste0("mol_premrna_", fid)
          model$feature_info$mol_mrna <- paste0("mol_mrna_", fid)
          model$feature_info$mol_protein <- paste0("mol_protein_", fid)
          
          molecule_ids <- c(
            model$feature_info$mol_premrna, 
            model$feature_info$mol_mrna,
            model$feature_info$mol_protein
          )
          
          initial_state <- set_names(
            rep(0, length(molecule_ids)),
            molecule_ids
          )
          
          # extract params
          parameters <- .kinetics_extract_parameters(
            model$feature_info, 
            model$feature_network
          )
          
          # determine variables to be used during burn in
          burn_variables <- 
            model$feature_info %>% 
            filter(burn) %>% 
            select(mol_premrna, mol_mrna, mol_protein) %>% 
            gather(col, val) %>% 
            pull(val)
          
          # return system
          model$simulation_system <- lst(
            reactions = formulae, 
            molecule_ids,
            initial_state,
            parameters,
            burn_variables
          )
          
          model
        }
        
        #' @export
        #' @rdname generate_kinetics
        #' @importFrom stats runif
        kinetics_default <- function() {
          # satisfy r cmd check
          transcription_rate <- translation_rate <- mrna_halflife <- protein_halflife <- 
            independence <- splicing_rate <- effect <- strength <- hill <- NULL
          
          sampler_tfs <-
            function(feature_info, feature_network, cache_dir = NULL, verbose = FALSE) {
              feature_info %>% mutate(
                transcription_rate = transcription_rate %|% runif(n(), 10, 20),
                translation_rate = translation_rate %|% runif(n(), 100, 150),
                mrna_halflife = mrna_halflife %|% runif(n(), 2.5, 5),
                protein_halflife = protein_halflife %|% runif(n(), 5, 10),
                independence = independence %|% 1,
                splicing_rate = splicing_rate %|% (log(2) / 2)
              )
            }
          
          sampler_interactions <-
            function(feature_info, feature_network, cache_dir = NULL, verbose = FALSE) {
              feature_network %>% 
                mutate(
                  effect = effect %|% sample(c(-1L, 1L), n(), replace = TRUE, prob = c(.25, .75)),
                  strength = strength %|% 10 ^ runif(n(), log10(1), log10(100)),
                  hill = hill %|% rnorm_bounded(n(), 2, 2, min = 1, max = 10)
                )
            }
          
          lst(sampler_tfs, sampler_nontfs = sampler_tfs, sampler_interactions)
        }
        
        
        ############# modify this function: sampler_tfs -> sampler_nontfs in 
        ############# 'sample kinetics from dataset for non-tfs' section
        .kinetics_generate_gene_kinetics <- function(model) {
          # satisfy r cmd check
          is_tf <- mrna_halflife <- protein_halflife <- to <- feature_id <- effect <- 
            basal <- basal_2 <- num_molecules <- mult <- id <- NULL
          
          if (model$verbose) cat("Generating kinetics for ", nrow(model$feature_info), " features\n", sep = "")
          params <- model$kinetics_params
          
          # fetch feature info and network
          feature_info <- model$feature_info %>%
            .kinetics_add_columns(c("transcription_rate", "splicing_rate", "translation_rate", "mrna_halflife", "protein_halflife", "independence"), NA_real_)
          feature_network <- model$feature_network %>%
            .kinetics_add_columns(c("effect", "strength", "hill"), NA_real_)
          
          # generate relatively stable kinetics for TFs
          feature_info_tf <- params$sampler_tfs(
            feature_info %>% filter(is_tf),
            feature_network
          )
          
          # sample kinetics from dataset for non-tfs
          feature_info_nontf <- params$sampler_nontfs(
            feature_info %>% filter(!is_tf),
            feature_network, 
            cache_dir = model$download_cache_dir,
            verbose = model$verbose
          )
          
          # combine feature info
          feature_info <- 
            bind_rows(feature_info_tf, feature_info_nontf) %>% 
            mutate(
              mrna_decay_rate = log(2) / mrna_halflife,
              protein_decay_rate = log(2) / protein_halflife
            )
          
          # sample network kinetics from distributions
          feature_network <- params$sampler_interactions(feature_info, feature_network)
          
          # calculate k
          dis_out <- .kinetics_calculate_dissociation(feature_info, feature_network)
          feature_info <- dis_out$feature_info
          feature_network <- dis_out$feature_network
          
          # calculate ba and a
          feature_info <- 
            left_join(
              feature_info,
              feature_network %>% 
                rename(feature_id = to) %>% 
                group_by(feature_id) %>% 
                summarise(
                  basal_2 = .kinetics_calculate_basal(effect)
                ),
              by = "feature_id"
            ) %>% 
            mutate(
              # 1 for genes that are not being regulated by any other genes,
              # yet did not already have a value for 'basal' defined
              basal = basal %|% basal_2 %|% 1 
            ) %>% 
            select(-basal_2)
          
          model$feature_info <- feature_info
          model$feature_network <- feature_network
          
          model
        }
        
        #' @importFrom GillespieSSA2 reaction
        .kinetics_generate_formulae <- function(model) {
          # satisfy r cmd check
          from <- to <- `.` <- NULL
          
          
          if (model$verbose) cat("Generating formulae\n")
          
          # add helper information to feature info
          feature_info <- 
            model$feature_info %>% 
            left_join(
              model$feature_network %>% 
                group_by(feature_id = to) %>% 
                do({tibble(regulators = list(.))}) %>% 
                ungroup(),
              by = "feature_id"
            ) %>% 
            left_join(
              model$feature_network %>% 
                group_by(feature_id = from) %>% 
                summarise(num_targets = n()),
              by = "feature_id"
            )
          
          # generate formula per feature
          out <- furrr::future_map(
            seq_len(nrow(feature_info)),
            .progress = model$verbose,
            function(i) {
              info <- feature_info %>% extract_row_to_list(i)
              
              fid <- info$feature_id
              
              w <- paste0("mol_premrna_", fid)
              x <- paste0("mol_mrna_", fid)
              y <- paste0("mol_protein_", fid)
              
              transcription_rate <- paste0("transcription_rate_", fid)
              splicing_rate <- paste0("splicing_rate_", fid)
              translation_rate <- paste0("translation_rate_", fid)
              mrna_decay_rate <- paste0("mrna_decay_rate_", fid)
              protein_decay_rate <- paste0("protein_decay_rate_", fid)
              
              basal <- paste0("bas_", fid)
              independence <- paste0("ind_", fid)
              
              if (!is.null(info$regulators)) {
                rid <- info$regulators$from
                eff <- info$regulators$effect
                str <- info$regulators$strength
                reg_ys <- paste0("mol_protein_", rid)
                reg_diss <- paste0("dis_", rid, "_", fid)
                reg_hills <- paste0("hill_", rid, "_", fid)
                reg_strs <- paste0("str_", rid, "_", fid)
                regulation_var <- paste0("chi_", rid, "_", fid)
                
                reg_affinity_calc <- paste(paste0(regulation_var, " = ", reg_strs, " * pow(", reg_ys, "/", reg_diss, ", ", reg_hills, "); "), collapse = "")
                
                # Several optimisations have been applied.
                #
                # original:
                #   [ba + x0 + x0x1 + x1] / [x0 + x0x1 + x1 + 1],
                #   with xi = (yi / ki) ^ ci
                #
                # factorise:
                #   [ba + (x0 + 1) * (x1 + 1) - 1] / (x0 + 1) / (x1 + 1) / (x2 + 1)
                #
                # use buffer to remember calculations:
                #   [ba - 1 + buf0 * buf1] / buf0 / buf1 / buf2,
                # with buf0 = x0 + 1, buf1 = x1 + 1, buf2 = x2 + 1
                
                numerator <-
                  if (sum(eff > 0) > 0) {
                    paste0(basal, " - pow(", independence, ",", sum(eff > 0), ") + ", paste("(", regulation_var[eff > 0], " + ", independence, ")", collapse = " * ", sep = ""))
                  } else {
                    basal
                  }
                denominator <- paste("(", regulation_var, " + 1)", collapse = " * ", sep = "")
                
                act_function <- paste0(reg_affinity_calc, transcription_rate, " * (", numerator, ")/(", denominator, ")")
              } else {
                act_function <- paste0(transcription_rate, " * ", basal)
                regulation_var <- character()
              }
              
              formulae <- list(
                # pre-mRNA production
                reaction(
                  name = paste0("transcription_", fid),
                  effect = set_names(1, w),
                  propensity = paste0(act_function)
                ),
                # splicing
                reaction(
                  name = paste0("splicing_", fid),
                  effect = set_names(c(1, -1), c(x, w)),
                  propensity = paste0(splicing_rate, " * ", w)
                ),
                # protein production
                reaction(
                  name = paste0("translation_", fid), 
                  effect = set_names(1, y),
                  propensity = paste0(translation_rate, " * ", x)
                ),
                # pre-mRNA degradation
                reaction(
                  name = paste0("premrna_degradation_", fid),
                  effect = set_names(-1, w),
                  propensity = paste0(mrna_decay_rate, " * ", w)
                ),
                # mRNA degradation
                reaction(
                  name = paste0("mrna_degradation_", fid),
                  effect = set_names(-1, x),
                  propensity = paste0(mrna_decay_rate, " * ", x)
                ),
                # protein degradation
                reaction(
                  name = paste0("protein_degradation_", fid),
                  effect = set_names(-1, y),
                  propensity = paste0(protein_decay_rate, " * ", y)
                )
              )
              
              formulae[[1]]$buffer_ids <- regulation_var
              
              formulae
            }
          )
          
          unlist(out, recursive = FALSE)
        }
        
        .kinetics_extract_parameters <- function(feature_info, feature_network) {
          # satisfy r cmd check
          feature_id <- transcription_rate <- splicing_rate <- translation_rate <- mrna_decay_rate <- protein_decay_rate <-
            basal <- independence <- param <- value <- id <- from <- to <- dissociation <- hill <- strength <- `.` <- from <- NULL
          
          # extract production / degradation rates, ind and bas
          feature_params <- 
            feature_info %>% 
            select(feature_id, transcription_rate, splicing_rate, translation_rate, mrna_decay_rate, protein_decay_rate, bas = basal, ind = independence) %>% 
            gather(param, value, -feature_id) %>% 
            mutate(id = paste0(param, "_", feature_id)) %>% 
            select(id, value) %>% 
            deframe()
          
          # extract dis, hill, str
          edge_params <- 
            feature_network %>% 
            select(from, to, dis = dissociation, hill, str = strength) %>% 
            gather(param, value, -from, -to) %>% 
            mutate(id = paste0(param, "_", from, "_", to)) %>% 
            select(id, value) %>% 
            deframe()
          
          c(feature_params, edge_params)
        }
        
        .kinetics_calculate_basal <- function(effects) {
          case_when(
            all(effects == -1) ~ 1,
            all(effects == 1) ~ 0.0001,
            TRUE ~ 0.5
          )
        }
        
        .kinetics_calculate_dissociation <- function(feature_info, feature_network) {
          # satisfy r cmd check
          transcription_rate <- mrna_decay_rate <- splicing_rate <- max_premrna <- translation_rate <- 
            protein_decay_rate <- max_mrna <- feature_id <- max_protein <- NULL
          
          remove <- c("max_premrna", "max_mrna", "max_protein", "dissociation", "k", "max_protein")
          
          feature_info <- feature_info[, !colnames(feature_info) %in% remove]
          feature_network <- feature_network[, !colnames(feature_network) %in% remove]
          
          feature_info <- 
            feature_info %>%
            mutate(
              max_premrna = transcription_rate / (mrna_decay_rate + splicing_rate),
              max_mrna = splicing_rate / mrna_decay_rate * max_premrna,
              max_protein = translation_rate / protein_decay_rate * max_mrna
            )
          
          feature_network <- 
            feature_network %>% 
            left_join(feature_info %>% select(from = feature_id, max_protein), by = "from") %>% 
            mutate(
              dissociation = max_protein / 2
            )
          
          lst(feature_info, feature_network)
        }
        
        .kinetics_add_columns <- function(df, colnames, fill = NA) {
          for (colname in colnames) {
            if (!colname %in% colnames(df)) {
              df[[colname]] <- fill
            }
          }
          df
        }
        
        
        
      }
      model <- generate_kinetics1(model)
      # model[["feature_info"]][["mrna_halflife"]]
      # model[["feature_info"]][["transcription_rate"]]
      # model[["feature_info"]][["splicing_rate"]]
      model <- generate_gold_standard(model)
      model <- generate_cells(model)  ## 2 min ## 12 min for cell-specific network
      
      # check identity in 'model' level
      # identical(model0_1,model0_2)
      # identical(model0_1,model1_1)
      # find_diff(model0_1,model0_2)
      # find_diff(model0_1,model1_1)
      
      
      save(model,file = paste(dir,'dataset',data_num,'_net',net_id,"_para",para_id,'.RData',sep = ''))
      print(paste('dataset',data_num,'_net',net_id,"_para",para_id,sep = ''))
    }
    
  }
}





# Network inference with GENEIE3   # 2021/2/23
library(GENIE3)
#scale data
scale.matrix = function(a){
  for(i in 1:nrow(a)){
    if(!identical(as.numeric(a[i,]),numeric(length = ncol(a)))){      #considering special case: zero vector
      a[i,] = scale(a[i,])
    }
  }
  return(a)
}
for(data_num in 1:4){
  for(net_id in 1:n_net){
    for(para_id in 1:n_para){
      for(rep_id in 1:3){
        # load input matrix
        load(file = paste(dir,'dataset',data_num,'_net',net_id,"_para",para_id,'.RData',sep = ''))
        simu_counts = model[["simulations"]][["counts"]]
        n_genes = ncol(simu_counts)/3
        genes = model[["feature_info"]]$feature_id
        t = model[["simulations"]][["meta"]][["sim_time"]]
        
        sub.t = which(t>=0)   # select cells with t>=0
        
        if(rep_id!=0){  # if rep_id ==0 -> use all cells to infer network
          # select one trajectory of simulation
          sub.t = sub.t[which(sub.t>((rep_id-1)*3)*length(t)/num_simu & sub.t<=(rep_id*3)*length(t)/num_simu)]
        }
        subcells = sub.t         # select sub cells
        simu_counts = simu_counts[subcells,]
        n_genes = ncol(simu_counts)/3
        genes = model[["feature_info"]]$feature_id
        
        s = t(as.matrix(simu_counts[,(1:n_genes)+n_genes]))
        u = t(as.matrix(simu_counts[,1:n_genes]))
        rownames(s) = rownames(u) = genes
        
        # remove constant mRNA
        non_const_id = NULL
        for(i in 1:nrow(s)){
          if(sd(s[i,])!=0&sd(u[i,])!=0){
            non_const_id = c(non_const_id,i)
          }
        }
        s = s[non_const_id,]
        u = u[non_const_id,]
        
        #log transform
        u = log2(u+1)
        s = log2(s+1)
        #scale
        scaleintron = scale.matrix(u)
        scaleexon = scale.matrix(s)
        #intron names:  add '-intron'
        rownames(scaleintron) = paste(rownames(scaleintron),'-intron',sep = '')
        #combine exon and intron
        scaledata = rbind(scaleexon,scaleintron)
        
        # TF information
        feature_info = model[["feature_info"]]
        TFs = as.character(as.matrix(feature_info[which(feature_info$is_tf==TRUE),1]))
        TFs = intersect(TFs, rownames(s))
        
        for(seed_id in 1:1){
          # intron and exon input
          input = scaledata
          #run genie3
          exprMatrix = input
          set.seed(seed_id)
          weightMatrix <- GENIE3(exprMatrix,regulators = TFs, nCores = 20) 
          print(paste('dataset',data_num,'_net',net_id,"_para",para_id,"_rep",rep_id,"_seed",seed_id," inex done!",sep = ''))
          save(weightMatrix, file = paste(dir,'Network_inference/genie3_output_inex_',data_num,'_net',net_id,"_para",para_id,"_",
                                          rep_id,"_",seed_id,'.RData',sep = ""))
          
          # exon as input
          input = scaleexon
          #run genie3
          exprMatrix = input
          set.seed(seed_id)
          weightMatrix <- GENIE3(exprMatrix,regulators = TFs, nCores = 20)  #10 min
          print(paste('dataset',data_num,'_net',net_id,"_para",para_id,"_rep",rep_id,"_seed",seed_id," exon done!",sep = ''))
          save(weightMatrix, file = paste(dir,'Network_inference/genie3_output_exon_',data_num,'_net',net_id,"_para",para_id,"_",
                                          rep_id,"_",seed_id,'.RData',sep = ""))
        }
        
      } 
    }
  }
}


# Evaluation for network inference   # 2021/2/23
list_to_mat = function(a){
  TFs = unique(a[,1])
  targets = unique(a[,2])
  mat = matrix(0,nrow = length(TFs),ncol = length(targets))
  rownames(mat) = sort(TFs)
  colnames(mat) = sort(targets)
  if(ncol(a)==2){
    for(i in 1:nrow(a)){
      mat[a[i,1],a[i,2]] = 1
    }
  }else{
    for(i in 1:nrow(a)){
      mat[a[i,1],a[i,2]] = as.numeric(a[i,3])
    }
  }
  return(mat)
}
Evaluation <- function(A1,Agr1,thr){    #calculate the results of A1 with threshold, and compare to Agr1 to calculate precision, recall and FPR
  B = matrix(0,nrow=nrow(A1),ncol = ncol(A1)) 
  for(i in 1:nrow(A1)){
    for(j in 1:ncol(A1)){
      if(A1[i,j]>=thr) B[i,j]<- 1
    }
  } 
  a=b=c=d=0 
  for(i in 1:nrow(B)){
    for(j in 1:ncol(B)){
      if(i!=j){				#attention that here ignored the autoregulation
        if(Agr1[i,j]==1&&B[i,j]==1){a=a+1
        }else if(Agr1[i,j]==1&&B[i,j]==0) { b=b+1
        }else if(Agr1[i,j]==0&&B[i,j]==1) { c=c+1
        }else if(Agr1[i,j]==0&&B[i,j]==0) { d=d+1}
      }
    }
  }	
  if((a+c)==0){
    return(c(0,0,c/(c+d)))
  }else{
    return(c(a/(a+c),a/(a+b),c/(c+d))) ##c(precision,recall,FPR)
  }
}
df = NULL
for(data_num in 1:4){
  for(net_id in 1:n_net){
    for(para_id in 1:n_para){
      for(rep_id in 1:3){
        for(seed_id in 1:1){
          # load NI result 
          if(T){
            # load exon result
            load(file = paste(dir,'Network_inference/genie3_output_exon_',data_num,'_net',net_id,"_para",para_id,"_",
                              rep_id,"_",seed_id,'.RData',sep = ""))
            exon.mat = weightMatrix
            
            # load inex result
            load(file = paste(dir,'Network_inference/genie3_output_inex_',data_num,'_net',net_id,"_para",para_id,"_",
                              rep_id,"_",seed_id,'.RData',sep = ""))
            mat = weightMatrix
            # extract intron targets
            targets = colnames(mat)[stringr::str_detect(colnames(mat),"intron")]
            mat = mat[,targets]
            # change intron names
            colnames(mat) = unlist(strsplit(colnames(mat),"-intron"))
            # remove self regulation
            for(i in 1:nrow(mat)){
              if(is.element(rownames(mat)[i],colnames(mat))){
                mat[i,rownames(mat)[i]]=0
              }
            }
            inex.mat = mat
          }
          
          # check NaN
          id = unique(c(which(is.nan(exon.mat[1,])),which(is.nan(inex.mat[1,]))))
          if(length(id)!=0){
            exon.mat = exon.mat[,-id]
            inex.mat = inex.mat[,-id]
          }
          
          
          # ground-truth network
          if(T){
            #dyngen network
            load(file = paste(dir,'dataset',data_num,'_net',net_id,"_para",para_id,'.RData',sep = ''))
            net = model[["feature_network"]]
            reference = as.matrix(net[,1:2])
            gt = list_to_mat(reference)
          }
          
          
          # change size -- select submatrix according to TFs
          if(T){
            dim(exon.mat)
            dim(inex.mat)
            dim(gt)
            feature_info = model[["feature_info"]]
            TFs = as.character(as.matrix(feature_info[which(feature_info$is_tf==TRUE),1]))
            TFs = intersect(rownames(exon.mat),TFs)
            TFs = intersect(rownames(gt),TFs)
            targets = intersect(colnames(gt),colnames(exon.mat))
            gt = gt[TFs,targets]
            dim(gt)
            sum(as.vector(gt))
            exon.mat = exon.mat[TFs,targets]
            inex.mat = inex.mat[TFs,targets]
          }
          
          # Evaluation
          print(paste(data_num,net_id,para_id,rep_id,seed_id,"Start!!!",sep = " "))
          mat = exon.mat
          if(T){
            A1 = mat
            #construct Agr matrix, regulation exists:1, non-exists:0
            Agr1 <- gt
            
            #calculate precisino, recall and FPR
            interval= 100
            precision = recall = FPR =  rep(0,times = interval+1)
            
            for(i in 1:(interval+1)){			#equal interval for threshold
              thr = (i-1 )/interval*max(A1)
              result = Evaluation(A1,Agr1,thr)
              precision[i] = result[1]
              recall[i] = result[2]
              FPR[i] = result[3]
              # print(i/interval)
            }
            w <- cbind((0:interval)/interval*max(A1),precision,recall)		#compare precision and recall in different threshold
            recall.exon = recall
            precision.exon = precision
            
            AUPR = 0
            for(i in 1:(interval)){AUPR = AUPR+ (recall[i] - recall[i+1])*(precision[i]+precision[i+1])/2}
            AUROC = 0
            for(i in 1:interval){AUROC = AUROC+ (FPR[i] - FPR[i+1])*(recall[i]+recall[i+1])/2}
            random_precision = precision[1]
            
            # calculate EPR
            density_gt = sum(gt)/length(gt)
            thr = quantile(as.vector(A1),1-density_gt)
            result = Evaluation(A1,Agr1,thr)
            EP = result[1]   #precision
            EPR = EP / precision[1]
            # points(result[1],result[2],col = 'red',pch =16,cex = 2)
            
            
            # df:  data_num,net_id,num_genes,rep_id, seed_id,input, EPR, AUROC, AUPR, EP, random precision,TP
            df.tmp = c(data_num,net_id,num_genes,para_id,rep_id,seed_id,"exon", EPR, AUROC, AUPR, EP, random_precision)
            df = rbind(df,df.tmp)
            
            print(paste(data_num,net_id,para_id,rep_id,seed_id,"exon","Done!!!",sep = " "))
          }
          
          # Evaluation for inex.mat
          mat = inex.mat
          if(T){
            A1 = mat
            #construct Agr matrix, regulation exists:1, non-exists:0
            Agr1 <- gt
            
            #calculate precisino, recall and FPR
            interval= 100
            precision = recall = FPR =  rep(0,times = interval+1)
            # thr =  (0:interval )/interval*max(A1)
            # result = sapply(thr,Evaluation,A1 = A1, Agr1 = Agr1)
            
            for(i in 1:(interval+1)){			#equal interval for threshold
              thr = (i-1 )/interval*max(A1)
              result = Evaluation(A1,Agr1,thr)
              precision[i] = result[1]
              recall[i] = result[2]
              FPR[i] = result[3]
              # print(i/interval)
            }
            w <- cbind((0:interval)/interval*max(A1),precision,recall)		#compare precision and recall in different threshold
            recall.inex = recall
            precision.inex = precision
            
            # plot(recall.inex,precision.inex,type = 'b',main = paste(data_num,num_cells,rep_id,sep = " "),col = "red",
            #      ylim = c(0,max(precision.inex,precision.exon)))
            # abline(h = precision[1], col = 'blue')
            # points(recall.exon,precision.exon,type = "b",col = "black")
            # abline(h = precision.exon[1], col = 'blue')
            
            AUPR = 0
            for(i in 1:(interval)){AUPR = AUPR+ (recall[i] - recall[i+1])*(precision[i]+precision[i+1])/2}
            AUROC = 0
            for(i in 1:interval){AUROC = AUROC+ (FPR[i] - FPR[i+1])*(recall[i]+recall[i+1])/2}
            random_precision = precision[1]
            
            # calculate EPR
            density_gt = sum(gt)/length(gt)
            thr = quantile(as.vector(A1),1-density_gt)
            result = Evaluation(A1,Agr1,thr)
            EP = result[1]   #precision
            EPR = EP / precision[1]
            x.inex = result[1]
            y.inex = result[2]
            # points(x.inex,y.inex,col = 'red',pch =16,cex = 2)
            # points(x.exon,y.exon,col = 'black',pch =16,cex = 2)
            
            # df:  data_num,net_id,num_genes,rep_id, seed_id,input, EPR, AUROC, AUPR, EP, random precision,TP
            df.tmp = c(data_num,net_id,num_genes,para_id,rep_id,seed_id,"inex", EPR, AUROC, AUPR, EP, random_precision)
            df = rbind(df,df.tmp)
            
            print(paste(data_num,net_id,para_id,rep_id,seed_id,"inex","Done!!!",sep = " "))
            
          }
        } 
      }
    }
  }
}


colnames(df) = c("data_num","net_id","num_genes","para_id","rep_id", "seed_id","input", "EPR", "AUROC", "AUPR", "EP",
                 "random precision")
rownames(df) = NULL
df = data.frame(data_num = as.numeric(df[,1]),net_id = as.numeric(df[,2]),num_genes = as.numeric(df[,3]),
                para_id = as.numeric(df[,4]),rep_id = as.numeric(df[,5]),seed_id = as.numeric(df[,6]),
                input = as.character(df[,7]),EPR = as.numeric(df[,8]),
                AUROC = as.numeric(df[,9]),AUPR = as.numeric(df[,10]), EP = as.numeric(df[,11]),
                random_precision = as.numeric(df[,12]))

save(df,file = paste(dir,"evaluation/variation_of_config_para_all_backbones.RData",sep = ""))




