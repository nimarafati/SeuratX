#' Read zUMIs Expression Data
#'
#' Reads and loads zUMIs expression data from RDS files, based on metadata provided in a tab-separated file.
#'
#' @param metadata_path The directory path where the metadata file is located.
#' @param sample_name_column The column name in the metadata file that contains sample identifiers.
#' @param data_path_column The column name in the metadata file that contains the path to sample data.
#' @return A list containing the expression data for each sample.
#' @export
#' @examples
#' expression_list <- read_zUMIs_expression("metadata.tsv", "path/to/metadata", "sample_id", "data_path")

read_zUMIs_expression <- function(metadata_path, sample_name_column, data_path_column) {

  # Read metadata
  metadata <- tryCatch(
    read.csv(metadata_path, sep = '\t'),
    error = function(e) stop("Error reading metadata file: ", e$message)
  )
  # Check passed parameters

  ## Check if metadata file exists
  if (!file.exists(metadata_path)) {
    stop("Error: Metadata file not found at ", metadata_path)
  }

  ## Check if sample_name_column exists
  if (!sample_name_column %in% colnames(metadata)) {
    stop("Error: Column '", sample_name_column, "' not found in metadata.")
  }

  ## Check if data_path_column exists
  if (!data_path_column %in% colnames(metadata)) {
    stop("Error: Column '", data_path_column, "' not found in metadata.")
  }



  # Extract sample names
  samples <- metadata[[sample_name_column]]

  # Initialize an empty list to store expression data
  expression_list <- list()

  # Read expression data for each sample
  lapply(samples, function(x) { cat('\r',x)
    tmp_path <- file.path(metadata[[data_path_column]][metadata[[sample_name_column]] == x], paste0("zUMIs_output/expression/", x, ".dgecounts.rds"))

    if (!file.exists(tmp_path)) {
      warning("Warning: File not found: ", tmp_path)
      return(NULL)
    }

    tmp <- tryCatch(
      readRDS(tmp_path),
      error = function(e) {
        warning("Error reading RDS file for sample ", x, ": ", e$message)
        return(NULL)
      }
    )

    expression_list[[x]] <<- tmp
  })

  return(expression_list)
}


#' Export zUMIs Expression Matrix
#'
#' Extracts quantification matrices from a list of loaded zUMIs expression data and generates a summary dataframe.
#'
#' @param expression_list A list containing zUMIs expression data for different samples.
#' @param quant_type The type of quantification data to extract. Must be one of \code{"umicount"} or \code{"readcount"}.
#' @param feature_type A character vector specifying which feature types to extract (e.g., \code{"exon"}, \code{"inex"}, \code{"intron"}).
#' @return A list containing:
#'   \item{df_nwell}{A summary dataframe with the number of wells per feature.}
#'   \item{quant_list}{A nested list containing quantification matrices per sample and feature.}
#' @export
#' @examples
#' quant_list <- export_expression_matrix(expression_list, "umicount", c("exon", "intron", "inex"))
export_expression_matrix <- function(expression_list, quant_type, feature_type) {
  # Check passed parameters

  # expression_list
  if (length(expression_list) == 0) {
    stop("Error: expression_list is empty. Provide a list of loaded zUMIs RDS files from read_zUMIs_expression.")
  }

  # quant_type
  quant_type_vec <- c("umicount", "readcount")
  if (!quant_type %in% quant_type_vec) {
    stop("Error: Specified quantification type '", quant_type, "' does not exist. It should be one of: ", paste(quant_type_vec, collapse=", "))
  }

  # feature_type
  feature_type_vec <- c("exon", "inex", "intron")
  if (!all(feature_type %in% feature_type_vec)) {
    stop("Error: Specified feature type '", feature_type, "' is incorrect. It should be one of: ", paste(feature_type_vec, collapse=", "))
  }

  # Initializing lists
  quant_list <- list()

  # Loop over samples
  for (smpl in names(expression_list)) {
    cat("Parsing sample", smpl, "\n")
    tmp_zumis <- expression_list[[smpl]]

    for (feat in feature_type) {
      quant_list[[smpl]][[feat]] <- tmp_zumis[[quant_type]][[feat]][['all']]
    }
  }

  # Return both outputs as a list
  return(quant_list = quant_list)
}


#' Read and Process a GTF File
#'
#' This function reads a GTF file, extracts gene annotations, and saves the processed data to a CSV file.
#' If a preprocessed gene info file exists, it is loaded instead. `data_genes` is derived dynamically
#' from `quant_list` based on the provided `feature_type`.
#'
#' @param gtf_file Path to the GTF file.
#' @param gene_info_file Path to the gene info file (if available, it will be used instead of parsing the GTF file).
#' @param quant_list A list containing gene expression matrices for different samples.
#' @param feature_type Character vector specifying which feature types to extract from `quant_list`.
#' @return A data frame containing processed gene annotations.
#' @export
#' @examples
#' gtf_data <- read_gtf_file("path/to/gtf.gtf", "path/to/gene_info.csv", quant_list, feature_type = c("exon", "intron"))
read_gtf_file <- function(gtf_file, gene_info_file = 'zUMIs_annot.csv', quant_list, feature_type) {
  # Validate input parameters
  if (!is.list(quant_list) || length(quant_list) == 0) {
    stop("Error: quant_list must be a non-empty list.")
  }

  if (!is.character(feature_type) || length(feature_type) == 0) {
    stop("Error: feature_type must be a non-empty character vector.")
  }

  # Derive unique gene names from quant_list
  data_genes <- unique(unlist(lapply(names(quant_list), function(smpl) {
    unique(unlist(lapply(feature_type, function(ft) {
      rownames(quant_list[[smpl]][[ft]])
    })))
  })))

  # Check if gene info file exists
  if (file.exists(gene_info_file)) {
    cat("Reading gene info file:", gene_info_file, "\n")
    gtf <- tryCatch(
      read.csv(gene_info_file),
      error = function(e) stop("Error reading gene info file: ", e$message)
    )
  } else {
    cat("Gene info file not found. Processing GTF file:", gtf_file, "\n")

    # Check if GTF file exists
    if (!file.exists(gtf_file)) {
      stop("Error: GTF file not found at ", gtf_file)
    }

    # Read the GTF file
    gtf <- tryCatch(
      read.table(gtf_file, header = FALSE, sep = "\t", quote = ""),
      error = function(e) stop("Error reading GTF file: ", e$message)
    )

    # Filter only "gene" entries or ERCC genes
    gtf <- gtf[gtf[, 3] == "gene" | grepl("ERCC", gtf[, 1]), ]
    colnames(gtf) <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

    # Helper function to extract attributes
    extract_attributes <- function(gtf_attributes, att_of_interest) {
      att <- strsplit(gtf_attributes, "; ")
      att <- gsub("\"", "", unlist(att))
      att <- gsub(";", "", att)
      match_idx <- grep(att_of_interest, att)

      if (length(match_idx) > 0) {
        return(unlist(strsplit(att[match_idx], " "))[2])
      } else {
        return(NA)
      }
    }

    # Extract gene attributes
    gene_biotype <- unname(sapply(gtf$attributes, extract_attributes, "gene_biotype"))
    gtf$gene_biotype <- unname(sapply(gtf$attributes, extract_attributes, "gene_biotype"))
    gtf$gene_name <- unname(sapply(gtf$attributes, extract_attributes, "gene_name"))
    gtf$gene_id <- unname(sapply(gtf$attributes, extract_attributes, "gene_id"))


    # Fix ERCC gene names
    gtf[grepl("ERCC", gtf$chr), "gene_name"] <- gtf[grepl("ERCC", gtf$chr), "chr"]

    # Ensure unique gene names
    gtf$uniq.name <- NA
    matched_genes <- match(data_genes, gtf$gene_id)
    unique_gene_names <- make.unique(gtf$gene_name[matched_genes])
    gtf$uniq.name[matched_genes] <- unique_gene_names
    gtf$gene_biotype <- gene_biotype
    # Save processed GTF to file
    write.csv(gtf, file = gene_info_file)
  }

  return(gtf)
}


#' Read zUMIs UMI Counts
#'
#' Reads UMI count files for each sample specified in the metadata and returns them as a list.
#'
#' @param metadata_path Path to the metadata file (TSV format).
#' @param sample_name_column The column name in the metadata file that contains sample identifiers.
#' @param data_path_column The column name in the metadata file that contains the path to sample data.
#' @return A named list where each element is a data frame of UMI counts per sample.
#' @export
#' @examples
#' umi_counts_list <- read_zUMIs_UMIcounts("metadata.tsv", "sample_id", "data_path")
read_zUMIs_UMIcounts <- function(metadata_path, sample_name_column, data_path_column) {

  # Read metadata
  metadata <- tryCatch(
    read.csv(metadata_path, sep = '\t'),
    error = function(e) stop("Error reading metadata file: ", e$message)
  )

  # Validate metadata columns
  if (!sample_name_column %in% colnames(metadata)) {
    stop("Error: Column '", sample_name_column, "' not found in metadata.")
  }

  if (!data_path_column %in% colnames(metadata)) {
    stop("Error: Column '", data_path_column, "' not found in metadata.")
  }

  # Extract sample names
  samples <- metadata[[sample_name_column]]

  # Initialize an empty list to store UMI counts
  umi_counts_list <- list()

  # Read UMI counts for each sample
  lapply(samples, function(x) { cat('\r', x)
    umi_path <- file.path(metadata[[data_path_column]][metadata[[sample_name_column]] == x],
                          paste0("zUMIs_output/stats/", x, ".UMIcounts.txt"))

    if (!file.exists(umi_path)) {
      warning("Warning: File not found: ", umi_path)
      return(NULL)
    }

    umi_counts <- tryCatch(
      read.table(umi_path, header = TRUE, sep = "\t"),
      error = function(e) {
        warning("Error reading UMI counts file for sample ", x, ": ", e$message)
        return(NULL)
      }
    )

    umi_counts_list[[x]] <<- umi_counts
  })

  return(umi_counts_list)
}

#' Read zUMIs Gene Counts
#'
#' Reads gene count files for each sample specified in the metadata and returns them as a list.
#'
#' @param metadata_path Path to the metadata file (TSV format).
#' @param sample_name_column The column name in the metadata file that contains sample identifiers.
#' @param data_path_column The column name in the metadata file that contains the path to sample data.
#' @return A named list where each element is a data frame of gene counts per sample.
#' @export
#' @examples
#' gene_counts_list <- read_zUMIs_genecounts("metadata.tsv", "sample_id", "data_path")
read_zUMIs_genecounts <- function(metadata_path, sample_name_column, data_path_column) {

  # Read metadata
  metadata <- tryCatch(
    read.csv(metadata_path, sep = '\t'),
    error = function(e) stop("Error reading metadata file: ", e$message)
  )

  # Validate metadata columns
  if (!sample_name_column %in% colnames(metadata)) {
    stop("Error: Column '", sample_name_column, "' not found in metadata.")
  }

  if (!data_path_column %in% colnames(metadata)) {
    stop("Error: Column '", data_path_column, "' not found in metadata.")
  }

  # Extract sample names
  samples <- metadata[[sample_name_column]]

  # Initialize an empty list to store gene counts
  gene_counts_list <- list()

  # Read gene counts for each sample
  lapply(samples, function(x) { cat('\r', x)
    gene_path <- file.path(metadata[[data_path_column]][metadata[[sample_name_column]] == x],
                           paste0("zUMIs_output/stats/", x, ".genecounts.txt"))

    if (!file.exists(gene_path)) {
      warning("Warning: File not found: ", gene_path)
      return(NULL)
    }

    gene_counts <- tryCatch(
      read.table(gene_path, header = TRUE, sep = "\t"),
      error = function(e) {
        warning("Error reading gene counts file for sample ", x, ": ", e$message)
        return(NULL)
      }
    )

    gene_counts_list[[x]] <<- gene_counts
  })

  return(gene_counts_list)
}

#' Read zUMIs Well Barcodes
#'
#' Reads well barcode files for each sample specified in the metadata and returns them as a list.
#'
#' @param metadata_path Path to the metadata file (TSV format).
#' @param sample_name_column The column name in the metadata file that contains sample identifiers.
#' @param data_path_column The column name in the metadata file that contains the path to sample data.
#' @return A named list where each element is a vector of well barcodes per sample.
#' @export
#' @examples
#' well_barcodes_list <- read_zUMIs_wellbarcodes("metadata.tsv", "sample_id", "data_path")
read_zUMIs_wellbarcodes <- function(metadata_path, sample_name_column, data_path_column) {

  # Read metadata
  metadata <- tryCatch(
    read.csv(metadata_path, sep = '\t'),
    error = function(e) stop("Error reading metadata file: ", e$message)
  )

  # Validate metadata columns
  if (!sample_name_column %in% colnames(metadata)) {
    stop("Error: Column '", sample_name_column, "' not found in metadata.")
  }

  if (!data_path_column %in% colnames(metadata)) {
    stop("Error: Column '", data_path_column, "' not found in metadata.")
  }

  # Extract sample names
  samples <- metadata[[sample_name_column]]

  # Initialize an empty list to store well barcodes
  well_barcodes_list <- list()
  #uuu77b/6b&
  # Readu well barcodes for each sample
  lapply(samples, function(x) { cat('\r', x)
    barcode_path <- file.path(metadata[[data_path_column]][metadata[[sample_name_column]] == x],
                              paste0("zUMIs_output/stats/", x, ".well_barcodes.txt"))
    #bv7b
    if (!file.exists(barcode_path)) {
      warning("Warning: File not found: ", barcode_path)
      return(NULL)
    }

    well_barcodes <- tryCatch(
      well_barcodes <- read.csv(barcode_path, sep = '\t', skip = 1),
      error = function(e) {
        warning("Error reading well barcode file for sample ", x, ": ", e$message)
        return(NULL)
      }
    )
    colnames(well_barcodes) <- c('SampleID', 'WellID')


    well_barcodes_list[[x]] <<- well_barcodes
  })

  return(well_barcodes_list)
}

#' Create and merge Seurat objects from zUMIs data
#'
#' This function creates individual Seurat objects from a list containing quantification matrices, UMI counts, and gene counts
#' for each sample. It then merges these objects into a single Seurat object. Note that this function requires a global data frame
#' \code{gtf} with columns \code{gene_id} and \code{uniq.name} to correctly map gene identifiers.
#'
#' @param data_list A list containing three elements:
#'   \item{quant_list}{A named list of quantification matrices (one per sample).}
#'   \item{umicount_list}{A named list of UMI count data frames (one per sample).}
#'   \item{genecount_list}{A named list of gene count data frames (one per sample).}
#' @param feature_type A character string specifying the feature type. Valid options are "exon", "intron", or "inex".
#'
#' @return A merged Seurat object combining the data from all samples.
#'
#' @details The function extracts the column corresponding to the specified \code{feature_type} from each sample's quantification data.
#' It then remaps the row names using the global \code{gtf} data frame. For each sample, the function filters the UMI and gene counts
#' based on the selected feature type and combines these into sample-specific metadata, which is added to the Seurat object.
#'
#' @examples
#' \dontrun{
#'   data_list <- list(
#'     quant_list = quantifications,
#'     umicount_list = umi_counts,
#'     genecount_list = gene_counts
#'   )
#'   merged_obj <- create_merged_seurat(data_list, "exon")
#' }
#'
#' @importFrom Seurat CreateSeuratObject merge
#' @export
create_merged_seurat <- function(data_list, feature_type, gtf_file, gene_info_file) {
  quant_list <- data_list$quant_list
  umicount_list <- data_list$umicount_list
  genecount_list <- data_list$genecount_list

  gtf <- read_gtf_file(gtf_file = gtf_file, gene_info_file = gene_info_file, quant_list = quant_list, feature_type = feature_type)

  # Determine the appropriate feature type label for filtering
  if (feature_type == 'exon') {
    feat_label <- 'Exon'
  } else if (feature_type == 'intron') {
    feat_label <- 'Intron'
  } else if (feature_type == 'inex') {
    feat_label <- 'Intron+Exon'
  } else {
    stop("Invalid feature_type. Must be one of: 'exon', 'intron', 'inex'.")
  }

  seurat_list <- list()

  for (smpl in names(quant_list)) {
    # Subset quantification data based on feature_type column
    tmp_quant <- data.frame(quant_list[[smpl]][[feature_type]])


    # Map gene IDs to unique names using the global gtf data frame
    m <- match(rownames(tmp_quant), gtf$gene_id)
    if (any(is.na(m))) {
      warning("Some gene IDs were not found in gtf for sample ", smpl)
    }
    rownames(tmp_quant) <- gtf$uniq.name[m]

    # Retrieve sample-specific gene and UMI count data
    genecounts_smpl <- genecount_list[[smpl]]
    umicounts_smpl <- umicount_list[[smpl]]

    # Filter the counts based on the specified feature type
    genecounts_smpl <- genecounts_smpl[(genecounts_smpl$type == feat_label),]
    names(genecounts_smpl)[names(genecounts_smpl) == 'Count'] <- 'gene_count'
    umicounts_smpl <- umicounts_smpl[(umicounts_smpl$type == feat_label),]
    names(umicounts_smpl)[names(umicounts_smpl) == 'Count'] <- 'umi_count'

    # Combine the UMI and gene count data into sample metadata
    smpl_metadata <- cbind(umicounts_smpl, genecounts_smpl)
    smpl_metadata <- subset(smpl_metadata, select = -c(type, SampleID))
    rownames(smpl_metadata) <- smpl_metadata$SampleID

    # Create a Seurat object using the quantification data and metadata
    tmp_obj <- CreateSeuratObject(counts = tmp_quant, project = smpl, meta.data = smpl_metadata)
    seurat_list[[smpl]] <- tmp_obj
  }

  # Merge the individual Seurat objects
  if (length(seurat_list) == 0) {
    stop("No Seurat objects were created.")
  } else if (length(seurat_list) == 1) {
    ser_obj <- seurat_list[[1]]
  } else {
    ser_obj <- merge(x = seurat_list[[1]], y = seurat_list[-1])
  }

  rm(tmp_quant, genecounts_smpl, umicounts_smpl, smpl_metadata, tmp_obj)
  gc()
  genes.table <- gtf[match(rownames(ser_obj), gtf$uniq.name),]
  result_list <- list(adata=ser_obj, genes.table=genes.table)
  return(result_list)
}


#' Read zUMIs Data
#'
#' This function reads expression, UMI count, gene count, and well barcode data from zUMIs output files
#' based on metadata information and exports the expression matrix.
#'
#' @param metadata_file_name Character. Name of the metadata file.
#' @param metadata_dir Character. Directory where the metadata file is located.
#' @param sample_name_column Character. Column name in the metadata file containing sample names.
#' @param data_path_column Character. Column name in the metadata file containing data paths.
#'
#' @return A list containing:
#' \item{expression_list}{List. Loaded expression data.}
#' \item{umicount_list}{List. Loaded UMI count data.}
#' \item{genecount_list}{List. Loaded gene count data.}
#' \item{wellbarcode_list}{List. Loaded well barcode data.}
#' \item{quant_list}{List. Exported expression matrix.}
#'
#' @examples
#' \dontrun{
#' data_list <- read_zUMIs_data(
#'   metadata_file_name = "metadata.csv",
#'   metadata_dir = "path/to/metadata",
#'   sample_name_column = "SampleID",
#'   data_path_column = "DataPath"
#' )
#' }
#'
#' @export
read_zUMIs_data <- function(metadata_file_name,
                            metadata_dir,
                            sample_name_column,
                            data_path_column){
  print('Loading expression files...')
  expression_list <- read_zUMIs_expression(metadata_path = file.path(metadata_dir,  metadata_file_name),
                                           sample_name_column = sample_name_column,
                                           data_path_column = data_path_column)
  print('Exporting expression matrix...')
  quant_list <- export_expression_matrix(expression_list, "readcount", feature_type = c("exon", "intron", "inex"))


  print('Loading umicount files...')
  umicount_list <- read_zUMIs_UMIcounts(metadata_path = file.path(metadata_dir,  metadata_file_name),
                                        sample_name_column = sample_name_column,
                                        data_path_column = data_path_column)
  print('Loading genecount files...')
  genecount_list <- read_zUMIs_genecounts(metadata_path = file.path(metadata_dir,  metadata_file_name),
                                          sample_name_column = sample_name_column,
                                          data_path_column = data_path_column)
  print('Loading wellbarcodes files...')
  wellbarcode_list <- read_zUMIs_wellbarcodes(metadata_path = file.path(metadata_dir,  metadata_file_name),
                                              sample_name_column = sample_name_column,
                                              data_path_column = data_path_column)
  data_list <- list(expression_list = expression_list,
                    umicount_list = umicount_list,
                    genecount_list = genecount_list,
                    wellbarcode_list = wellbarcode_list,
                    quant_list = quant_list)
  return(data_list)
}

#' Calculate Statistics for Single-Cell RNA-seq Data
#'
#' This function calculates statistics for single-cell RNA sequencing data,
#' including percentages of ERCC spike-ins, mitochondrial reads, and rRNA reads.
#' It also optionally filters out these gene categories.
#'
#' @param adata A Seurat object containing single-cell RNA-seq data.
#' @param filter_ercc Logical; if TRUE, ERCC spike-ins are removed. Default is TRUE.
#' @param filter_rRNA Logical; if TRUE, rRNA genes are removed. Default is TRUE.
#' @param filter_mito Logical; if TRUE, mitochondrial genes are removed. Default is TRUE.
#' @param mito_chr_name Character; the chromosome name corresponding to mitochondrial genes. Default is "MT".
#' @param genes_table A data frame containing gene metadata, with at least
#'   columns "gene_name" and "chr" (for mitochondrial genes) and "gene_biotype" (for rRNA genes).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{adata}{A modified Seurat object with updated metadata.}
#'   \item{genes.table}{A filtered gene annotation table corresponding to the retained genes.}
#' }
#'
#' @examples
#' # Create a mock dataset for testing
#' library(Seurat)
#' adata <- CreateSeuratObject(counts = matrix(rnorm(100), nrow = 10))
#' genes_table <- data.frame(Gene = paste0("Gene", 1:10))
#'
#' # Example usage:
#' result <- calculate_stats(adata, filter_ercc = TRUE, filter_rRNA = TRUE, filter_mito = TRUE, mito_chr_name = "MT", genes_table)
#' @export
calculate_stats <- function(adata, filter_ercc = TRUE, filter_rRNA = TRUE, filter_mito = TRUE, mito_chr_name = "MT", genes_table){
  ## To do:
  ## - Add percent intron

  ## Calculate stats
  nC <- colSums(adata@assays$RNA@counts)
  M <- data.frame(percentERCC = rep(0,ncol(adata)),row.names = colnames(adata))

  ## ERCC
  ercc <- grepl("^ERCC", rownames(adata), ignore.case = T)
  nERCC <- colSums(adata@assays$RNA@counts[ercc,])
  if(sum(nERCC) > 0){
    percentERCC <- nERCC/nC*100
  }


  ## Mito
  mito <- genes.table$gene_name[genes.table$chr == mito_chr_name]
  if(sum(!is.na(mito)) >0){
    mito <- mito[!is.na(mito)]
    nMT <- colSums(adata@assays$RNA@counts[(rownames(adata@assays$RNA@counts) %in% mito), ])
    percentMT <- nMT/(nC-nERCC)*100
  }else{
    nMT <- rep(0,ncol(adata))
    percentMT <- 0
  }



  ## rRNA
  #rRNA <- genes.table$gene_name[genes.table$gene_biotype == 'rRNA']
  rRNA <- genes.table$gene_name[grep('rRNA', genes.table$gene_biotype)]
  if(sum(!is.na(rRNA)) >0){
    rRNA <- rRNA[!is.na(rRNA)]
    rRNA_genes <- genes.table$gene_name[genes.table$gene_biotype == "rRNA" & !(genes.table$gene_name %in% mito)]
    rRNA_genes <- intersect(rRNA_genes, rownames(adata@assays$RNA@counts))  # Ensure genes exist in matrix
    nR <- colSums(adata@assays$RNA@counts[rRNA_genes, ])

    percentrRNA <- nR/(nC-nMT-nERCC)*100
  }else{
    nR <- rep(0,ncol(adata))
    percentrRNA <- 0
  }


  ## Protein coding
  protein_coding <- genes.table$gene_name[genes.table$gene_biotype == 'protein_coding']

  protein_coding_genes <- protein_coding[!protein_coding %in% mito]
  protein_coding_genes <- intersect(protein_coding_genes, rownames(adata@assays$RNA@counts))
  nPC <- colSums(adata@assays$RNA@counts[protein_coding_genes, ])

  percentPC <- nPC/(nC-nMT-nERCC)*100

  ## Ribosomal protein.
  rpc <- grepl("^Rp[ls]", rownames(adata), ignore.case = T)
  nRPC <- colSums(adata@assays$RNA@counts[rpc, ])
  percentRiboPC <- nRPC/(nC-nMT-nERCC)*100

  ##

  if(filter_ercc == TRUE && sum(nERCC) >0){
    adata <- adata[!ercc,]
  }

  if(filter_mito == TRUE && sum(!is.na(mito)) >0){
    adata <- adata[!rownames(adata) %in% mito,]
  }

  if(filter_rRNA == TRUE && sum(!is.na(rRNA)) >0){
    adata <- adata[!rownames(adata) %in% rRNA,]
  }
  ## stats matrix
  M$percentERCC <- percentERCC
  M$percentrRNA <- percentrRNA
  M$percentMT <- percentMT
  M$nCount_ERCC <- nERCC
  M$read_counts <- nC
  M$Plate <- adata$orig.ident
  M$percentPC <- percentPC
  M$percentRiboPC<- percentRiboPC

  adata <- AddMetaData(adata, M)


  genes.table <- genes.table[genes.table$gene_name %in% rownames(adata),]
  output <- list(adata=adata, genes.table=genes.table)
  return(output)
  gc()
}

#' Compute and Visualize Relative Gene Expression Per Cell
#'
#' This function computes the relative expression of each gene per cell using sparse matrix operations
#' and generates a boxplot of the most expressed genes.
#'
#' @param seurat_obj A Seurat object containing single-cell RNA-seq data.
#' @param output_path Character; the directory path where the output PNG image will be saved.
#' @param top_gene Integer; the number of top expressed genes to visualize. Default is 20.
#'
#' @return None (the function generates and saves a boxplot). The function also prints the output file path
#'   and the gene with the highest variation among the selected top genes.
#'
#' @examples
#' # Example usage:
#' compute_relative_expression(seurat_obj, output_path = "../../results/01-QC", top_gene = 20)
#'
#' @export
compute_relative_expression <- function(seurat_obj, output_path, top_gene = 20) {

  # Extract count matrix
  C <- seurat_obj@assays$RNA@counts

  # Compute relative expression per cell
  C <- Matrix::t(Matrix::t(C) / Matrix::colSums(C)) * 100

  # Identify the top 20 most expressed genes
  most_expressed <- order(apply(C, 1, median), decreasing = TRUE)[top_gene:1]

  # Save the boxplot as a PNG file
  file_path <- file.path(output_path, 'Top_Expressed_genes.png')
  png(file_path, width = 2000, height = 1000, res = 200)
  par(mar = c(4, 10, 2, 1))

  # Generate and print the boxplot
  p <- boxplot(t(as.matrix(C[most_expressed, ])),
               cex = 0.1,
               las = 1,
               xlab = "% total count per cell",
               col = (scales::hue_pal())(top_gene)[top_gene:1],
               horizontal = TRUE)
  dev.off()

  # Garbage collection
  gc()

  top_var_gene_name <- rownames(C)[most_expressed[length(most_expressed)]]
  print(paste("Boxplot was saved in:", file_path))
  print(paste(top_var_gene_name, "shows the highest variation among the top", top_gene, "genes."))
}

#' Generate and Save QC Violin Plots with Optional Cutoff Lines
#'
#' This function generates violin plots for specified features in a Seurat object,
#' optionally including cutoff values as horizontal lines, and saves them as PNG images.
#' The plots are arranged in grids of 2x2 per image for efficient visualization.
#'
#' @param seurat_obj A Seurat object containing the single-cell RNA-seq data.
#' @param features A character vector of feature names to be plotted.
#' @param output_path A string specifying the directory where the output images
#'        should be saved. If the directory does not exist, it will be created.
#' @param show_cutoff A logical value indicating whether to display cutoff values
#'        as horizontal error bars (default: FALSE).
#' @param x_var A string specifying the variable to be used for the x-axis (e.g., sample identity).
#' @param fill_var A string specifying the variable to be used for the fill color.
#'        If empty (`''`), it defaults to `x_var`.
#' @param cutoff_val A numeric vector specifying the cutoff values corresponding
#'        to each feature. If empty, no cutoffs are applied.
#'
#' @return A named list of ggplot objects corresponding to the generated violin plots.
#' @export
#'
#' @examples
#' \dontrun{
#' qc_plot(seurat_obj, features = c("nFeature_RNA", "percent.mt"),
#'         output_path = "qc_plots", show_cutoff = TRUE,
#'         x_var = "orig.ident", fill_var = "sample_group",
#'         cutoff_val = c(1000, 5))
#' }
qc_plot <- function(seurat_obj, features, output_path, show_cutoff = FALSE, x_var = orig.ident, fill_var = '', cutoff_val = ''){
  library(gridExtra)
  plt_list <- list()
  final_features_vec <- c()
  meta <- seurat_obj@meta.data
  if(fill_var == ''){
    fill_var <- x_var
  }
  # Create feature_cutoff_df if both provided
  if(any(cutoff_val != '') && show_cutoff == TRUE && is.numeric(cutoff_val)){
    cutoff_df <- data.frame(features = features, cutoff_val = cutoff_val)
    rownames(cutoff_df) <- cutoff_df$features
  }else{
    cutoff_df <- data.frame(features = features)
    rownames(cutoff_df) <- cutoff_df$features
  }

  for(feat in cutoff_df$features){
    if(sum(!is.na(seurat_obj@meta.data[,feat]))>0){
      final_features_vec <- c(final_features_vec,feat)
      if(show_cutoff == TRUE && is.numeric(cutoff_df[feat,'cutoff_val'])){
        cutoff_val_feat <- cutoff_df[feat, "cutoff_val"]

        meta$cutoff_val <- cutoff_val_feat # Adding cutoff_val_feat to meta df because neither passing the values directly nor saved in new variable helped. It was using the last value for geom_errorbar.
        plt_list[[feat]] <- ggplot(meta, aes(x = .data[[x_var]], y = .data[[feat]], fill = .data[[fill_var]])) +
          geom_violin(scale = "width", trim=TRUE, adjust = 1) +
          geom_errorbar(width = 0.8, aes(ymax = cutoff_val, ymin = cutoff_val)) +
          geom_jitter(height = 0, size = .1) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle(feat) +
          theme_classic() +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          NoLegend()


        #plt_list[[feat]] <- VlnPlot(seurat_obj, features = feat, pt.size = 0.1) + NoLegend()
      }else{
        plt_list[[feat]] <- ggplot(meta, aes(x = .data[[x_var]], y = .data[[feat]], fill = .data[[fill_var]])) +
          geom_violin(scale = "width", trim=TRUE, adjust = 1) +
          geom_jitter(height = 0, size = .1) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          ggtitle(feat) +
          theme_classic() +
          theme(plot.title = element_text(hjust = 0.5)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          NoLegend()

        #plt_list[[feat]] <- VlnPlot(seurat_obj, features = feat, pt.size = 0.1) + NoLegend()
      }
    }
  }

  n_features <- length(final_features_vec)
  if (n_features > 0) {
    # Split into groups of 4
    split_indices <- split(seq_len(n_features), ceiling(seq_along(seq_len(n_features)) / 4))

    # Create output directory if it doesn't exist
    if (!dir.exists(output_path)) {
      dir.create(output_path, recursive = TRUE)
    }

    # Iterate through each group of 4 plots and save them
    for (i in seq_along(split_indices)) {
      feature_subset <- final_features_vec[split_indices[[i]]]  # Select features for this group
      plot_subset <- plt_list[feature_subset]  # Select corresponding plots

      # Arrange plots into a grid
      grid_plot <- marrangeGrob(plot_subset, nrow = 2, ncol = 2, top = paste("QC Plots Group", i))

      # Save as PNG
      output_file <- file.path(output_path, paste0("qc_plot_group_", i, ".png"))
      ggsave(output_file, grid_plot, width = 10, height = 8, dpi = 300)
    }
  }

  return(plt_list)
}

#' Plot QC on Plate
#'
#' This function generates heatmap plots from metadata for multiple features and saves them as PNG files.
#' The filenames are determined based on the feature names provided.
#'
#' @param adata An object containing metadata with `WellID` and `Plate` information.
#' @param features A character vector representing the feature names, used for coloring the heatmap and naming the files.
#' @param output_dir A character string specifying the directory where the files will be saved (default: "plots").
#' @param width Numeric value for the width of the saved plot in inches (default: 8).
#' @param height Numeric value for the height of the saved plot in inches (default: 6).
#' @param dpi Numeric value specifying the resolution in dots per inch (default: 300).
#' @param facet_var A character string specifying the column name to be used for faceting.
#'
#' @return The function does not return a value; it saves the plots as PNG files.
#'
#' @examples
#' # plot_qc_on_plate(adata, c("percentMT", "nFeature_RNA"), facet_var = "Plate")
#'
#' @export
plot_qc_on_plate <- function(adata, features, output_dir = "plots", width = 8, height = 6, dpi = 300, facet_var = "Plate") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  m <- adata@meta.data
  m$column <- as.numeric(gsub("\\D", "", adata@meta.data$WellID))
  m$row <- gsub("\\d", "", adata@meta.data$WellID)

  for (feat in features) {
    if (!is.numeric(m[[feat]])) {
      warning("Skipping feature '", feat, "' as it is not numeric.")
      next
    }

    p <- ggplot(m, aes(x = column, y = row, fill = .data[[feat]])) +
      geom_tile() +
      scale_fill_gradientn(colors = c("grey", "yellow", "red", "black")) +
      facet_wrap(as.formula(paste("~", facet_var))) +
      theme_classic()

    file_name <- paste0(output_dir, "/plate_", feat, ".png")

    ggsave(filename = file_name, plot = p, width = width, height = height, dpi = dpi)

    message("Plot saved as: ", file_name)
  }
}

#' Filter Cells Based on Multiple Cut-off Criteria
#'
#' This function filters cells from a metadata data frame by applying thresholds on several metrics.
#' If a given cut-off is NA, the function calculates a dynamic threshold based on the median and MAD
#' (multiplied by \code{n.sd}) of the metric. The filtering criteria (greater or smaller than the cut-off)
#' is determined by the associated type.
#'
#' @param meta A data frame containing cell metadata, with row names representing cell IDs.
#' @param percentMT_cut Numeric value for the mitochondrial percentage threshold (default is 4).
#' @param percentMT_type Character indicating filtering direction for mitochondrial percentage ("greater" to filter cells with values above the cut-off; default is "greater").
#' @param percentERCC_cut Numeric value for the ERCC percentage threshold (default is 3).
#' @param percentERCC_type Character indicating filtering direction for ERCC percentage ("greater"; default is "greater").
#' @param percentPC_cut Numeric value for the percentPC threshold (default is 90).
#' @param percentPC_type Character indicating filtering direction for percentPC ("smaller" filters cells with values below the threshold; default is "smaller").
#' @param percentrRNA_cut Numeric value for the rRNA percentage threshold (default is 2).
#' @param percentrRNA_type Character indicating filtering direction for rRNA percentage ("greater"; default is "greater").
#' @param nFeature_RNA_cut Numeric value for the number of features threshold (default is 5000).
#' @param nFeature_RNA_type Character indicating filtering direction for number of features ("smaller"; default is "smaller").
#' @param gene_counts_cut Numeric value for the UMI counts threshold (default is 5000).
#' @param gene_counts_type Character indicating filtering direction for UMI counts ("smaller"; default is "smaller").
#' @param n.sd Numeric value representing the multiplier for the median absolute deviation (MAD) used to calculate the dynamic threshold when the cut-off is NA (default is 2).
#' @param plate Character specifying a plate identifier. Only cells belonging to this plate (i.e. where \code{meta$orig.ident == plate}) will be considered.
#'
#' @return A list with two elements:
#' \item{filtered}{A vector of cell IDs (row names from \code{meta}) that meet the filtering criteria.}
#' \item{settings}{A list of final cut-off values used for each metric.}
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   # Assume 'meta' is a data frame with the necessary columns.
#'   result <- cut.cells(meta, plate = "Plate1")
#'   print(result$filtered)
#'   print(result$settings)
#' }
#'
#' @export
cut_cells <- function(meta,
                      percentMT_cut = 10, percentMT_type = "greater",
                      percentERCC_cut = 25, percentERCC_type = "greater",
                      percentPC_cut = 75, percentPC_type = "smaller",
                      percentrRNA_cut = 2, percentrRNA_type = "greater",
                      nFeature_RNA_cut = 1000, nFeature_RNA_type = "smaller",
                      gene_count_cut = 1000, gene_count_type = "smaller",
                      n.sd = 2, plate = NULL) {

  # Restrict meta to a specific plate if provided
  meta <- meta[meta$orig.ident == plate, ]

  filt.cells <- list()
  settings.out <- list()

  # Helper function to process each metric
  process_metric <- function(metric, cut, type) {
    if (is.na(cut)) {
      m <- median(meta[, metric])
      s <- mad(meta[, metric])
      cut <- m + n.sd * s
      if (type == "smaller") {
        cut <- m - n.sd * s
      }
    }
    settings.out[[metric]] <<- cut  # Store computed cut-off
    if (type == "smaller") {
      filt.cells[[metric]] <<- rownames(meta)[meta[, metric] < cut]
    } else {
      filt.cells[[metric]] <<- rownames(meta)[meta[, metric] > cut]
    }
  }

  # Process each metric individually
  process_metric("percentMT", percentMT_cut, percentMT_type)
  process_metric("percentERCC", percentERCC_cut, percentERCC_type)
  process_metric("percentPC", percentPC_cut, percentPC_type)
  process_metric("percentrRNA", percentrRNA_cut, percentrRNA_type)
  process_metric("nFeature_RNA", nFeature_RNA_cut, nFeature_RNA_type)
  process_metric("gene_count", gene_count_cut, gene_count_type)

  all.filt <- unique(unlist(filt.cells))
  cat(sprintf("For plate %s filtering %d cells\n\n", plate, length(all.filt)))

  return(list(filtered = all.filt, settings = settings.out))
}
