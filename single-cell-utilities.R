
# quick view of large matrixes
corner = function (x, y = c(1:5)){
  x[y, y]
}

#imports peak x cell matrix as a dgCMatrix (sparse matrix)
import10X_peaks = function(tenx_peak_matrix_dir){
  mex_dir_path = tenx_peak_matrix_dir
  
  mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '/')
  feature_path <- paste(mex_dir_path, "peaks.bed", sep = '/')
  barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '/')
  
  features <- readr::read_tsv(feature_path, col_names = F) %>% tidyr::unite(feature)
  barcodes <- readr::read_tsv(barcode_path, col_names = F) %>% tidyr::unite(barcode)
  
  mtx <- Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$feature) %>%
    magrittr::set_colnames(barcodes$barcode)
  
  # return dgCmatrix
  peaks_mtx <- as(mtx, "dgCMatrix")
  peaks_mtx
}

# imports tf x cell matrix as a dgCMatrix (sparse matrix)
import10X_tf = function(tenx_tf_matrix_dir){
  tf_dir_path = tenx_tf_matrix_dir
  
  tf_path <- paste(tf_dir_path, "matrix.mtx", sep = '/')
  feature_path_tf <- paste(tf_dir_path, "motifs.tsv", sep = '/')
  barcode_path_tf <- paste(tf_dir_path, "barcodes.tsv", sep = '/')
  
  features_tf <- readr::read_tsv(feature_path_tf, col_names = c('feature', 'common_name'))
  barcodes_tf <- readr::read_tsv(barcode_path_tf, col_names = F) %>% tidyr::unite(barcode)
  
  tf_mtx <- Matrix::readMM(tf_path) %>%
    magrittr::set_rownames(features_tf$common_name) %>%
    magrittr::set_colnames(barcodes_tf$barcode)
  
  # return dgCmatrix
  tf_mtx <- as(tf_mtx, "dgCMatrix")
  tf_mtx
}

# SCTransform works on "RNA" assay.
# so make temp seurat object, add RNA slot, perform SCTransform.
# Returns Assay object [package "Seurat"] which can manually be added to existing Seurat object.
sctransform_assay = function(seurat_obj, input_mat, key_name){
  seurat_obj[["RNA"]] = CreateAssayObject(counts = input_mat)
  seurat_obj <- SCTransform(seurat_obj, return.only.var.genes = FALSE)
  #sct=GetAssay(seurat_obj, assay = "SCT")
  new_assay=seurat_obj[["SCT"]]
  Key(object = new_assay) <- paste0(key_name, "_")
  return(new_assay)
}

# takes as input a large filtered peak_by_barcode matrix and divides into n bins of evenly sized submatrices
# allows for performing matrix operations on matrices too large for memory, without having to used DelayedArrays
# output is the same matrix class (i.e. dense, CSC, etc) as input matrix
# output variable names are 'input_matirx' character with bin appended as suffix
split_matrix = function (input_matrix, bins){
  l = nrow(input_matrix)
  bin_size = floor( l / bins )
  
  for ( i in 0:(bins-1) ){
    print(i)
    nam <- paste(deparse(substitute(input_matrix)), i, sep = "_")
    print(nam)
    if ( i == 0 ) {
      start = 1
      end = bin_size - 1
      print( c(start, end) )
      mat = input_matrix[ c(start:end), ]
      assign(nam, mat, envir = .GlobalEnv)
    }
    if ( i > 0 & i < bins-1 ) {
      start = i * bin_size
      end = ( start + bin_size ) - 1
      print( c(start, end) )
      mat = input_matrix[ c(start:end), ]
      assign(nam, mat, envir = .GlobalEnv)
    }
    if ( i == bins-1) {
      start = i * bin_size
      end = max(l)
      print( c(start, end) )
      mat = input_matrix[ c(start:end), ]
      assign(nam, mat, envir = .GlobalEnv)
    }
  }
}

#input sparse matrices with duplicate rownames, add columns together/reduce to single unique rowname
aggregates = function(x){
  agg <- Matrix.utils::aggregate.Matrix(x,rownames(x),fun='sum')
  sp_mat <- as(agg, "dgCMatrix")
  sp_mat
}

#input granges with + or - strand, no * in strand col. Returns grange with coordinates reduced.
#if * in strand, breaks (debug later).
#useful for merging reference coordinates (i.e. transcripts or exons), but prob. not for experimental reads
#"merged" mcol returns either union of all overlapping coord if any are found
#otherwise, if no overlaps found, "merged" mcols returns row irange (self coord).
reduce_overlaping_granges = function(gr){
  merged = GenomicRanges::reduce(gr)
  merged
  keep.overlaps = GenomicRanges::findOverlaps(gr, merged)
  keep.overlaps
  subjectOverlaps = merged[S4Vectors::subjectHits(x = keep.overlaps)]
  subjectOverlaps
  GenomicRanges::mcols(subjectOverlaps) <- GenomicRanges::mcols(gr)
  GenomicRanges::mcols(subjectOverlaps)
  GenomicRanges::mcols(subjectOverlaps)$merged = paste(subjectOverlaps$gene_id,
                                                       GenomicRanges::seqnames(subjectOverlaps), 
                                                       GenomicRanges::start(subjectOverlaps), 
                                                       GenomicRanges::end(subjectOverlaps), sep = "_")
  # range/distance of merged coords should always be equal or greater than original interval
  subjectOverlaps
}

#used by Seurat::CreateGeneActivityScore() and myCreateGeneActivityMatrix()
Extend <- function(x, upstream = 0, downstream = 0) {
  if (any(GenomicRanges::strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- GenomicRanges::strand(x = x) == "+" | GenomicRanges::strand(x = x) == "*"
  new_start <- GenomicRanges::start(x = x) - ifelse(test = on_plus, yes = upstream, no = downstream)
  new_end <- GenomicRanges::end(x = x) + ifelse(test = on_plus, yes = downstream, no = upstream)
  IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
  x <- GenomicRanges::trim(x = x)
  return(x)
}

# opposite of Extend, for myCreateGeneActivityMatrix()
Narrow <- function(x, from.start = 0, from.end = 0) {
  if (any(GenomicRanges::strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- GenomicRanges::strand(x = x) == "+" | GenomicRanges::strand(x = x) == "*"
  new_start <- GenomicRanges::start(x = x) + ifelse(test = on_plus, yes = from.start, no = from.end)
  new_end <- GenomicRanges::end(x = x) - ifelse(test = on_plus, yes = from.end, no = from.start)
  IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
  x <- GenomicRanges::trim(x = x)
  return(x)
}

# originally loaded from Seurat::CreateGeneActivityScore
#later, calling PlanThreads() was remowed from myCreateGeneActivityScore
# PlanThreads <- function() {
#   nthreads <- eval(expr = formals(fun = future::plan())$workers)
#   # changed %||% to %/% (integer division)
#   return(nthreads %/% 1)
# }

# used by Seurat::CreateGeneActivityScore and myCreateGeneActivityMatrix()
PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find ",
      paste(pkgs[!package.installed], collapse = ', '),
      "; please install"
    )
  }
  invisible(x = package.installed)
}

# takes a sparse matrix and gtf file as input, returns a sparse matrix containing counts of features x cells
# to count all reads within a gene, choose feature_type = "gene"
# to count all transcript level reads or TSS start sites, choose feature_type = "transcript"
# to count all exonic reads, choose feature_type = "exon"
# Can also select to include promoter region 
#     i.e. Promoter sums would be: feature_type = "transcript", upstream = 1000, downstream = 100)
# Can also combine 2 calls to myCreateGeneActivityMatrix() followed by aggregates() (line 63 this file)
#     i.e. Promoters plus exons would be:
#     promoters = feature_type = "transcript", include.body = FALSE, upstream = 1000, downstream = 0)
#     exons = feature_type = "exon", include.body = TRUE, upstream = 0, downstream = 0)
#     exonActivityMatrix = aggregates( rbind(promoters, exons) )
myCreateGeneActivityMatrix <- function(
  peak.matrix,
  annotation.file,
  #note chormosome levels are for mouse!
  seq.levels = paste0("chr", c(1:19, "X", "Y")),
  #added feature type selection
  feature_type = "gene",
  include.body = TRUE,
  upstream = 2000,
  downstream = 0,
  verbose = TRUE,
  gene_type = "most",
  select_exons = "all",
  use_nearest_dist = TRUE
) {
  if (!PackageCheck('GenomicRanges', error = FALSE)) {
    stop("Please install GenomicRanges from Bioconductor.")
  }
  if (!PackageCheck('rtracklayer', error = FALSE)) {
    stop("Please install rtracklayer from Bioconductor.")
  }
  
  # convert peak matrix to GRanges object
  peak.df <- rownames(x = peak.matrix)
  # added "_"
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":|_", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", 'start', 'end')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  
  # get annotation file, select genes
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')
  GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
  
  if (gene_type == "all") {
    print("using 'all' 'gene_type' option.")
    gtf.genes = gtf
  } else if (gene_type == "most") {
    print("using 'most' 'gene_type' option.")
    # Use all gene/transcript types except non-experimentally confirmed and pseudogenes:
    pseudogene.types = unique(gtf$gene_type)[grep(".*pseudogene.*", unique(gtf$gene_type) )]
    remove.genes = c("TEC", pseudogene.types)
    remove.index = grepl(paste(remove.genes, collapse="|"), gtf$gene_type)
    gtf.genes <- gtf[ !remove.index ]
  } else if (gene_type == "some") {
    print("using 'some' 'gene_type' option.")
    # Alternatively, filter genes as 10X genomics does (i.e. protein coding, lincRNA, antisense)
    keep.genes = c("protein_coding", "lincRNA", "antisense")
    keep.index = grepl(paste(keep.genes, collapse="|"), gtf$gene_type)
    gtf.genes <- gtf[ keep.index ]
  } else if (gene_type == "protein") {
    # Protein coding only
    print("using 'protein' 'gene_type' option.")
    keep.genes = c("protein_coding")
    keep.index = grepl(paste(keep.genes, collapse="|"), gtf$gene_type)
    gtf.genes <- gtf[ keep.index ]
  } else {
    print("Please choose a valid 'gene_type' option.")
  }
  
  if (feature_type == "transcript") {
    #currently only use for promtoer sums/ TSS quantification, does not currently support multiple isoforms 
    gtf.genes <- gtf.genes[ gtf.genes$type == 'transcript' ]
  } else if (feature_type == "exon") {
    if (select_exons == "first") {
      print("using first exon only")
      gtf.genes <- gtf.genes[ gtf.genes$type == 'exon' ]
      gtf.genes <- gtf.genes[ gtf.genes$exon_number == '1' ]
      # gtf.genes = reduce_overlaping_granges(gtf.genes)
      # #remove dup exons (1 unique exon per coordinate) to avoid bias for genes with multiple isoforms
      # gtf.genes = gtf.genes[ !duplicated(gtf.genes$merged) ]
    } else if (select_exons == "all") {
      print("using all exons")
      gtf.genes <- gtf.genes[ gtf.genes$type == 'exon' ]
      # gtf.genes = reduce_overlaping_granges(gtf.genes)
      # #remove dup exons (1 unique exon per coordinate) to avoid bias for genes with multiple isoforms
      # gtf.genes = gtf.genes[ !duplicated(gtf.genes$merged) ]
    } else {
      print("please choose either 'first' or 'all' for select_exons")
    }
  } else {
    gtf.genes <- gtf.genes[ gtf.genes$type == 'gene' ]
  }
  
  # Extend definition up/downstream
  if ((include.body) & (feature_type=="exon") & (select_exons=="all")) {
    print("using all exons and promoter")
    # interexon dist can be less than "upstream" distance, so only extend first exon coord and merge with all other exons.
    exon_first <- gtf.genes[ gtf.genes$exon_number == '1' ]
    exon_first.body_prom <- Extend(x = exon_first, upstream = upstream, downstream = downstream)
    exons_remaining <- gtf.genes[ gtf.genes$exon_number != '1' ]
    gtf.body_prom = c(exon_first.body_prom, exons_remaining)
    
  } else if ((include.body) & (feature_type=="exon") & (select_exons=="first")) {
    print("using first exon with or without promoter")
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
    
  } else if ((include.body) & !(feature_type=="exon")) {
    print("using gene body/transcript with or without promoter")
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  } else {
    print("using only promoter.")
    # retain only first bp 
    gtf.genes = GenomicRanges::resize(gtf.genes, 1)
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
  }

  require(parallel)
  detectBatchCPUs <- function() { 
    ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
    if (is.na(ncores)) { 
        ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
    } 
    if (is.na(ncores)) { 
        return(4) # for helix
    } 
    return(ncores)
  }
  numCores = detectBatchCPUs()

  # reduce coordinates and remove duplicates.
  # Iterate per gene, otherwise adjacent genes on same strand are reduced into single coordinates
  if (feature_type != "gene"){
    # reduce coordinates and remove duplicates.
    # Iterate per gene, otherwise adjacent genes on same strand are reduced into single coordinates
    gene_targets = unique(gtf.body_prom$gene_id)
    print(length(gene_targets))
    print(head(gene_targets))

    gtf.list = parallel::mclapply(seq(gene_targets), function(i){
      current_gene = gtf.body_prom[ gtf.body_prom$gene_id %in% gene_targets[i] ];
      current_overlaps = reduce_overlaping_granges( current_gene )
      dedup_gene = current_overlaps[ !duplicated(current_overlaps$merged) ]
      dedup_gene
    },
    mc.cores = numCores);
    gtf.gr = do.call("c", gtf.list)
  } else {
    gtf.gr = gtf.body_prom
  }
  
  #for all exons, do distanceToNearest, otherwise, counts will be dependent on # of exons per gene
  if (use_nearest_dist){
    print("using distanceToNearest method")
    gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.gr)
    keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  } else {
    # for feature_type == "gene" or "transcript", deleted old gene.distances/keep.overlaps 
    # (does not count peaks overlapping 2 genes/promoters)
    #added new keep.overlaps (retains peaks that map to 2 or more genes (i.e. Dlx1/Dlx2))
    print("using findOverlaps method")
    keep.overlaps = GenomicRanges::findOverlaps(peaks.gr, gtf.gr)
  }

  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.gr[S4Vectors::subjectHits(x = keep.overlaps)]

  # Some GTF rows will not have gene_name attribute
  # Replace it by gene_id attribute
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]

  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  # changed ":" and "-" to "_"
  peak.ids$peak <- paste0(peak.ids$seqnames, "_", peak.ids$start, "_", peak.ids$end)
  #new way
  #not nec in current myCreateGeneActivity() because peak.matrix names were changed while
  #loading
  #peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
  annotations <- peak.ids[, c('peak', 'gene.name')]
  colnames(x = annotations) <- c('feature', 'new_feature')
  
  # collapse into expression matrix
  peak.matrix <- as(object = peak.matrix, Class = 'matrix')
  all.features <- unique(x = annotations$new_feature)
  
  # changed/removed  PlanThreads check
  mysapply <- ifelse(test = verbose, yes = pbapply::pbsapply, no = sapply)
  newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x){
    features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature
    submat <- peak.matrix[features.use, ]
    if (length(x = features.use) > 1) {
      return(Matrix::colSums(x = submat))
    } else {
      return(submat)
    }
  })
  newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  colnames(x = newmat) <- colnames(x = peak.matrix)
  return(as(object = newmat, Class = 'dgCMatrix'))
}

overwrite_gmat = function(snap, seurat, assay){
  print( c("dim of seurat matrix", dim(seurat[[assay]]@data)) )
  print( c("dim of snap matrix", dim(snap@gmat)) )
  
  print( "seurat head:" )
  print( seurat[[assay]]@data[1:5,1:5] )
  print( "snap head:" )
  print( snap@gmat[1:5,1:5] )
  
  # how many genes are repeated in activity matrices of Seurat and SnapATAC?
  print("num duplicates seurat:")
  print( length(rownames(seurat[[assay]]@data)[duplicated(rownames(seurat[[assay]]@data))]) )
  print("num duplicates snap:")
  print( length(colnames(x.sp.snap@gmat)[duplicated(rownames(x.sp.snap@gmat))]) )
  
  new_gmat = seurat[[assay]]@data
  colnames(new_gmat) = seurat@meta.data$barcode
  snap@gmat = t(new_gmat)
  return(snap)
}
