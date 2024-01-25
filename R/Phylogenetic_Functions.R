build_alignment <- function(datalist, fasta_output, 
                            ARB_DB = "~/PhD/Data_Storage/AlignmentDBs/SILVA138.1/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb", 
                            SINA_path = "/Users/felixmilke/opt/anaconda3/pkgs/sina-1.7.2-hadaa689_0/bin",
                            BMGE_path = "/Users/felixmilke/PhD/Software/BMGE/src/BMGE.jar",
                            BMGE = F) {
  
  require(reticulate)
  
  if (BMGE) {
    
    fasta_output_tmp <- str_replace_all(fasta_output, pattern = ".fasta", replacement = "_tmp.fasta")
    
    write_fasta(datalist, "output/tmp_fasta.fasta")
    
    writeLines(paste0("#!/bin/bash -i\n",
                      "eval \"$(conda shell.bash hook)\"\n",
                      "conda activate sina\n",
                      "\n",
                      "sina -i output/tmp_fasta.fasta",
                      " -r ", ARB_DB,
                      " -o ", fasta_output_tmp, "\n",
                      "\n",
                      "conda deactivate",
                      "\n\n",
                      "java -jar ", BMGE_path,
                      " -i ", fasta_output_tmp, " -t NT",
                      " -o ", fasta_output), "output/tmp_alignment_script.sh")
    
     system("bash output/tmp_alignment_script.sh")
    
  } else {
    
    write_fasta(datalist, fasta_output =  "output/tmp_fasta.fasta")
    
    writeLines(paste0("#!/bin/bash -i\n",
                      "eval \"$(conda shell.bash hook)\"\n",
                      "conda activate sina\n",
                      "\n",
                      "sina -i output/tmp_fasta.fasta",
                      " -r ", ARB_DB,
                      " -o ", fasta_output, "\n",
                      "\n",
                      "conda deactivate"), "output/tmp_alignment_script.sh")
    
    system("bash output/tmp_alignment_script.sh")
    
  }
  
  import_fasta(datalist, fasta_output, Type = Alignment) %>%
    return()
  
}

build_tree <- function(datalist, alignment_input = NULL, tree_output) {
  
  if (is.null(alignment_input)) {
    write_fasta(datalist, "output/tmp_fasta.fasta", Alignment)
    alignment_input = "output/tmp_fasta.fasta"
  }
  
  system(paste0("~/PhD/Software/FastTree/FastTree -gtr -nt ", alignment_input,
                " > ", tree_output))
  
  datalist$tree <- ape::read.tree(tree_output)
  
  return(datalist)
  
}

import_fasta <- function(datalist, fasta_location, Type = Sequence) {
  
  require(tidyverse)
  
  Type <- enquo(Type)
  
  tmp_seqs <- readLines(fasta_location)
  
  fasta_table <- tibble(Seq_ID = tmp_seqs[seq(1, length(tmp_seqs)-1, 2)],
                        !!quo_name(Type) := tmp_seqs[seq(2, length(tmp_seqs), 2)]) %>%
    mutate(Seq_ID = str_replace_all(Seq_ID, pattern = ">", replacement = "")) %>%
    dplyr::slice(match(datalist$Count_Data$OTU_ID, Seq_ID))
    
  if (nrow(fasta_table) != nrow(datalist$Count_Data)) {
    warning(paste0(Type, " file not complete: Expected ", nrow(datalist$Count_Data), 
                   " sequences, but ", nrow(fasta_table), " sequences matched to count-table."))
  }
  
  datalist$Sequence <- fasta_table
  
  return(datalist)

}

import_tree <- function(datalist, tree_location, root = T) {
  
  require(ape)
  
  tree_tmp <- ape::read.tree(tree_location) %>%
    keep.tip(tip = datalist$Count_Data$OTU_ID)
  
  if (root) {
    # tablify parts of tree that we need and
    # Take the longest terminal branch as outgroup
    new.outgroup <- tibble(Edge_1 = tree_tmp$edge[,1],
                     Edge_2 = tree_tmp$edge[,2],
                     Edge_Length = tree_tmp$edge.length) %>%
      slice(1:Ntip(tree_tmp)) %>%
      mutate(ID = tree_tmp$tip.label) %>%
      slice(which.max(Edge_Length)) %>%
      .$ID
    
    tree_tmp <- ape::root(tree_tmp, outgroup = new.outgroup, resolve.root = T)
    
    cat("Rooted tree with ASV ", new.outgroup, " as new root.")
    
  }
  
  return(tree_tmp)
  
}

write_fasta <- function(datalist, fasta_output, Type = Sequence) {
  
  require(tidyverse)
  
  Type <- enquo(Type)
  
  datalist$Sequence %>%
    select(Seq_ID, !!Type) %>%
    mutate(Seq_ID = paste0(">", Seq_ID)) %>%
    with(., do.call(rbind, lapply(seq(nrow(.)), function(i) t(.[i, ])))) %>%
    write.table(., file = fasta_output, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  return()
  
}

cophenetic_distance <- function(datalist, tree_location, use.cores=4) {
  
    require(furrr)
    
    plan(multisession, workers = use.cores)
  
    #Get tree in "postorder" order
    tree <- import_tree(datalist, tree_location, root = T) %>%
      ape::reorder.phylo(order="postorder")
    tip.order <- order(tree$tip.label)
    
    #Get "tip ages", or horizontal position of edges in re-ordered phylo object
    node.ages = ape::node.depth.edgelength(tree)
    my.combs <- gRbase::combn_prim(tip.order, 2) %>%
      t()
    
    #Get node paths to root for all tips
    my.node_paths <- furrr::future_map(seq(1, ape::Ntip(tree)), function (tip) {ape::nodepath(tree, from=tip, to = 1)})
    
    #Loop through all pairs of tips (in parallel) and calculate pairwise distance
    #=> d.cophenetic = (d[i, root] + d[j, root]) - 2 * d[mrca.ij, root]
    split.sums <- furrr::future_map2(my.combs[,1], my.combs[,2], function(g, h) {node.ages[g]+node.ages[h]})
    split.mrca <- furrr::future_map2_int(my.combs[,1], my.combs[,2], function(a, b) {
      max(intersect(my.node_paths[[a]], my.node_paths[[b]]))})
    
    dist.vec <- unlist(split.sums) - 2*node.ages[split.mrca]
    
    #Coerce into matrix
    mat.out <- matrix(nrow = ape::Ntip(tree), ncol = ape::Ntip(tree))
    diag(mat.out) <- 0
    rownames(mat.out) <- colnames(mat.out) <- tree$tip.label[tip.order]
    mat.out[upper.tri(mat.out)] <- dist.vec
    mat.out[lower.tri(mat.out)] <- t(mat.out)[lower.tri(mat.out)]
    
    return(mat.out)
}

UniFrac <- function(datalist, tree_location, distance="unweighted", blocksize=1000, use.cores = 4) {
  
  require(furrr)
  
  plan(multisession, workers = use.cores)
  
  tree <- import_tree(datalist, tree_location, root = T)
  
  Count_Data <- datalist$Count_Data %>%
    select_if(is.numeric) %>%
    as.matrix() %>%
    magrittr::set_rownames(datalist$Count_Data$OTU_ID) %>%
    .[tree$tip.label,]
  
  c.samples <- colnames(Count_Data)
  c.sample.sizes <- colSums(Count_Data)
  n.sample <- length(c.samples)
  
  #Get N x 2 matrix of pairwise combinations of samples
  combis <- utils::combn(colnames(Count_Data), 2, simplify=FALSE)
  
  # Create a list of descendants, starting from the first internal node (root)
  # Add the terminal edge descendants (tips). By definition, can only have one descendant
  desc_list <- c(as.list(1:length(tree$tip.label)), ape::prop.part(tree, check.labels = FALSE))
  
  # Convert `descList` to `edge_array` that matches the order of things in `tree$edge`
  # For each entry in the tree$edge table, sum the descendants for each sample
  # `tree$edge[i, 2]` is the node ID.
  edge.array <- furrr::future_map(1:nrow(tree$edge), function(i) {
    colSums(Count_Data[desc_list[[tree$edge[i, 2]]], ,drop=F], na.rm=T)
  }) %>%
    bind_rows()
  
  edge.occ <- (edge.array > 0) - 0
  edge.occ.list <- apply((edge.array > 0), 2, which)
  
  z = ape::reorder.phylo(tree, order="postorder")
  #Get "tip ages", or horizontal position of edges in re-ordered phylo object
  #Keep only tips
  tip.ages = ape::node.depth.edgelength(tree)[1:length(tree$tip.label)] %>%
    magrittr::set_names(z$tip.label) %>%
    .[rownames(Count_Data)]
  
  #Unweighted UniFrac
  size.chunk <- blocksize
  if (distance == "unweighted") {
    cs.list <- furrr::future_map(split(1:length(combis), ceiling(seq_along(1:length(combis))/size.chunk)), function(i) {
      
      map(combis[i], function(sp) {
        #Get occurrence vectors for current samples
        a <- edge.occ.list[[sp[1]]]
        b <- edge.occ.list[[sp[2]]]
        curr.shared_edges <- intersect(a, b)
        curr.union_edges <- union(a, b)
        #Return: 1 - (shared branch length) / (total branch length)
        1 - sum(tree$edge.length[curr.shared_edges], na.rm=TRUE) / sum(tree$edge.length[curr.union_edges], na.rm=TRUE)
      
        }) %>%
        unlist()
    })
  } else if (distance == "weighted") {
    cs.list <- furrr::future_map(split(1:length(combis), ceiling(seq_along(1:length(combis))/size.chunk)), function(i) {
      
      map(combis[i], function(sp) {
        #Get current samples and total sample sizes
        a <- sp[1]
        b <- sp[2]
        a.T <- c.sample.sizes[a]
        b.T <- c.sample.sizes[b]
        #Calculate branchweights
        curr.bw <- abs(edge.array[, a]/a.T - edge.array[, b]/b.T)
        #Return: (unshared branch length) / (total branch length)
        sum(tree$edge.length * curr.bw, na.rm = TRUE) / sum(tip.ages * (Count_Data[,a]/a.T + Count_Data[,b]/b.T), na.rm = TRUE)
      }) %>%
        unlist()
    })
  }
  
  #Pass output
  mat.UF <- matrix(NA_real_, n.sample, n.sample)
  diag(mat.UF) <- 0
  rownames(mat.UF) <- colnames(mat.UF) <- c.samples
  
  # Matrix-assign lower-triangle of UniFracMat. Then coerce to dist and return.
  mat.idx <- do.call(rbind, combis)[, 2:1]
  #Coerce results into matrices
  mat.UF[mat.idx] <- unlist(cs.list)
  mat.UF[mat.idx[,2:1]] <- unlist(cs.list)
  
  return(as.dist(mat.UF))
}
