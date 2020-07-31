#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# When training classifiers, avoid memory overflow by weighted down-sampling of cells.
# The most prevalent cell types are filtered out first, so less-represented cells are not under risk of being removed. 

option_list = list(
    make_option(
        c("-e", "--expression-data"),
        action = "store",
        default = NA,
        type = 'character',
        help = '10xGenomics-type directory holding expression matrix, genes, and barcodes'
    ),
    make_option(
        c("-m", "--metadata"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Metadata file mapping cells to cell types'
    ),
    make_option(
        c("-i", "--cell-id-field"),
        action = "store",
        default = "id",
        type = 'character',
        help = 'Name of cell id column in metada file'
    ),
    make_option(
        c("-c", "--cell-type-field"),
        action = "store",
        default = "inferred.cell.type",
        type = 'character',
        help = 'Name of cell type column in metada file'
    ),
    make_option(
        c("-l", "--array-size-limit"),
        action = "store",
        default = 2000000000,
        type = 'numeric',
        help = 'Maximum length of R array'
    ),
    make_option(
        c("-o", "--output-dir"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output directory for downsampled expression data'
    ),
    make_option(
        c("-n", "--metadata-upd"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Updated metadata file output path'
    )
)

opt = wsc_parse_args(option_list, mandatory = c('expression_data', 'metadata', 'output_dir', 'metadata_upd'))

if (is.na(opt$expression_data) || ! dir.exists(opt$expression_data)){
  stop("Provide a valid directory for --expression-data")
}else if (is.na(opt$metadata) || ! file.exists(opt$metadata)){
  stop("Provide a valid cell metadata file")
}

# Don't parse whole matrix up-front, we might not need to do anything...

print("Checking inputs for downsampling...")
expr_data = opt$expression_data
genes = read.csv(paste(expr_data, "genes.tsv", sep="/"), sep="\t", stringsAsFactors = FALSE, header = FALSE)
barcodes = read.csv(paste(expr_data, "barcodes.tsv", sep="/"), sep="\t", stringsAsFactors = FALSE, header = FALSE)
cell_num_limit = floor(opt$array_size_limit / nrow(genes))
current_cell_num = nrow(barcodes)

# if no down-samling is required return a special status code to let the user know
if(current_cell_num <= cell_num_limit){
  write("No downsampling required", stderr())
  quit(status = 2)
}
print("... done input checks, proceeding to downsampling")

# Okay, we do have to do something, so parse matrix properly

print("Parsing full matrix and metadata...")
suppressPackageStartupMessages(require(DropletUtils))
sce <- read10xCounts(opt$expression_data)

metadata <- read.csv(opt$metadata, sep = "\t", stringsAsFactors = FALSE)

for (field in c(opt$cell_id_field, opt$cell_type_field)){
  if (! field %in% colnames(metadata)){
    write(paste0("Supplied ID field: ", opt$cell_id_field, " not in metadata frame"), stderr())
    quit(status = 1)
  }
}
print("Done full parse")

# Put the metadata into the object so we can subset both together
colData(sce) <- merge(colData(sce), metadata, by.x='Barcode', by.y=opt$cell_id_field, all.x=TRUE, sort=FALSE)

# First candidates for removal are those without a label at all
print("Checking unlablled")
sce <- sce[, sce[[opt$cell_type_field]] != '']

# If we still have too many after removing unlabelled...

if (ncol(sce) > cell_num_limit ){ 

    print("... unlabelled removed (where applicable), we still need to downsample")

    # Only down-sample the most frequent cell types. Identify the ones to downsample
    # by progressively resetting the proportion of each type to that of the next
    # least abundant until total cell number falls below the limit.

    cell_type_freqs <- sampling_freqs <- sort(table(sce[[opt$cell_type_field]]), decreasing = TRUE)

    props <- cell_type_freqs/ sum(cell_type_freqs)
    classes_to_downsample <- c()

    for (i in 1:length(cell_type_freqs)){
      
      # Set the proportion for this cell type (and any preceding ones) to that
      # of the subsequent cell type

      props[1:i] <- props[i+1]
      classes_to_downsample <- c(classes_to_downsample, names(cell_type_freqs[i]))
      
      if ( sum(props * sum(cell_type_freqs)) < cell_num_limit ){
        break 
      }
    }

    no_cells_to_remove <- ncol(sce) - cell_num_limit
    cells_in_downsampled_groups <- sum(cell_type_freqs[classes_to_downsample])
    sampling_rate <- (cells_in_downsampled_groups - no_cells_to_remove)/ cells_in_downsampled_groups

    sampling_freqs[classes_to_downsample] <- floor(sampling_freqs[classes_to_downsample] * sampling_rate)

    # Now derive a cells list

    # These are the cells for gropus we don't need to sample
    unsampled <- sce$Barcode[sce[[opt$cell_type_field]] %in% names(cell_type_freqs)[! names(cell_type_freqs) %in% classes_to_downsample]]

    sampled <- unlist(lapply(classes_to_downsample, function(cd){
      sample(sce$Barcode[sce[[opt$cell_type_field]] == cd ], sampling_freqs[[cd]])
    }))

    selected_barcodes <- c(unsampled, sampled)
    sce <- sce[,sce$Barcode %in% selected_barcodes]
}else{
    print("... unlabelled removed (where applicable), no further downsampling required")
}

# write data
print("Writing outputs")
write10xCounts(opt$output_dir, sce, barcodes = sce$Barcode)
write.table(colData(sce)[,c(-1, -2)], opt$metadata_upd, sep="\t")
print("Outputs written successfully")
