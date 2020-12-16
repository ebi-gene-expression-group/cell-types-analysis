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
        help = '10xGenomics-type directory holding expression matrix, genes, 
                and barcodes'
    ),
    make_option(
        c("-m", "--metadata"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Metadata file mapping cells to cell types'
    ),
    make_option(
        c("-x", "--exclusions"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Path to the yaml file with excluded terms for 
        initial matrix filtering"
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
        c("-r", "--cell-count-threshold"),
        action = "store",
        default = 5,
        type = 'numeric',
        help = 'Threshold number of cells to keep a cell type in the matrix'
    ),
    make_option(
        c("-n", "--metadata-upd"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Updated metadata file output path'
    )
)

opt = wsc_parse_args(option_list, mandatory = c('expression_data', 'metadata',
                                                'output_dir', 'metadata_upd'))

if (is.na(opt$expression_data) || ! dir.exists(opt$expression_data)){
  stop("Provide a valid directory for --expression-data")
}else if (is.na(opt$metadata) || ! file.exists(opt$metadata)){
  stop("Provide a valid cell metadata file")
}

# Don't parse whole matrix up-front, we might not need to do anything...
print("Checking inputs for downsampling...")
expr_data = opt$expression_data
genes = read.csv(paste(expr_data, "genes.tsv", sep="/"), sep="\t", 
                 stringsAsFactors = FALSE, header = FALSE, check.names=FALSE)
barcodes = read.csv(paste(expr_data, "barcodes.tsv", sep="/"), sep="\t", 
                 stringsAsFactors = FALSE, header = FALSE, check.names=FALSE)
cell_num_limit = floor(opt$array_size_limit / nrow(genes))
current_cell_num = nrow(barcodes)

print(paste("Matrix limit of", opt$array_size_limit, 'for', nrow(genes),
            'genes implies cell number limit of', cell_num_limit))

# parse metadata file & check the fields
metadata <- read.csv(opt$metadata, sep = "\t", 
                     stringsAsFactors = FALSE, check.names=FALSE)

for (field in c(opt$cell_id_field, opt$cell_type_field)){
  if (! field %in% colnames(metadata)){
    write(paste0("Supplied ID field: ", opt$cell_id_field,
          " not in metadata frame"), stderr())
    quit(status = 1)
  }
}

metadata[[opt$cell_id_field]] <- gsub(' ', '_', metadata[[opt$cell_id_field]])

# check whether cells with low counts are present
cell_count_threshold <- opt$cell_count_threshold
cell_type_counts <- sort(table(metadata[[opt$cell_type_field]]), decreasing = FALSE)
minor_cell_types_present <- cell_type_counts[1] < cell_count_threshold

# check whether we need to down-sample at all
downsampling_required <- current_cell_num > cell_num_limit

# if no down-samling is required return a special status code to let the user know
if(!(downsampling_required | minor_cell_types_present)){
  write("No downsampling required; no minor cell types found", stderr())
  quit(status = 2)
}
print("... done input checks, proceeding to downsampling")

# Okay, we do have to do something, so parse matrix properly

print("Parsing full matrix...")
suppressPackageStartupMessages(require(DropletUtils))
suppressPackageStartupMessages(require(yaml))
sce <- read10xCounts(opt$expression_data)
print("Done full parse")

# Put the metadata into the object so we can subset both together

if (! any(sce$Barcode %in% metadata[[opt$cell_id_field]])){
  write(paste("Cannot match any cells to metadata using", 
               opt$cell_id_field), stderr())
  quit(status = 1)
}
colData(sce) <- merge(colData(sce), 
                      metadata[!duplicated(metadata[[opt$cell_id_field]]),],
                      by.x='Barcode', by.y=opt$cell_id_field, 
                      all.x=TRUE, sort=FALSE)

# If there are any NAs in the cell type field, set to empty string
colData(sce)[[opt$cell_type_field]][ is.na(colData(sce)[[opt$cell_type_field]]) ] <- ''

# Source function definitions
script_dir = dirname(strsplit(commandArgs()[grep('--file=', commandArgs())], '=')[[1]][2])
source(file.path(script_dir, 'cell_types_utils.R'))

if(downsampling_required){
    # First candidates for removal are those without a label at all
    print("Checking unlabelled")
    unlabelled <- ''
    if(! is.na(opt$exclusions)){
        e = yaml.load_file(opt$exclusions)
        unlabelled = c(unlabelled, tolower(e$unlabelled))
    }
    sce <- sce[, ! tolower(colData(sce)[[opt$cell_type_field]]) %in% unlabelled]

    # If we still have too many after removing unlabelled...

    if (ncol(sce) > cell_num_limit ){ 

        print("... unlabelled removed (where applicable), we still need to downsample")

        # Only down-sample the most frequent cell types. Identify the ones to downsample
        # by progressively resetting the proportion of each type to that of the next
        # least abundant until total cell number falls below the limit.

        cell_type_freqs <- sampling_freqs <- sort(table(colData(sce)[[opt$cell_type_field]]),
                                                  decreasing = TRUE)
        print("Starting cell type frequencies:")
        print(cell_type_freqs)

        props <- checkprops <- cell_type_freqs/ sum(cell_type_freqs)
        classes_to_downsample <- c()
        over_abundance <- c()

        for (i in 1:length(cell_type_freqs)){
          
          # Set the proportion for this cell type (and any preceding ones) to that
          # of the subsequent cell type

          checkprops[1:i] <- checkprops[i+1]
          over_abundance <- (props - checkprops)[1:i]

          classes_to_downsample <- c(classes_to_downsample, names(cell_type_freqs[i]))
          
          if ( sum(checkprops * sum(cell_type_freqs)) < cell_num_limit ){
            break 
          }
        }

        # For the cell types we downsample, do so in proportion to their relative
        # over-abundance

        no_cells_to_remove <- ncol(sce) - cell_num_limit
        cells_in_downsampled_groups <- sum(cell_type_freqs[classes_to_downsample])
        
        remove_props <- over_abundance / sum(over_abundance)    
        remove_freqs <- floor(no_cells_to_remove * remove_props)

        # Now derive a cells list

        # These are the cells for gropus we don't need to sample
        unsampled <- sce$Barcode[colData(sce)[[opt$cell_type_field]] %in% 
                     names(cell_type_freqs)[! names(cell_type_freqs) %in%
                     classes_to_downsample]]

        # Remove cells in proportion to their over-abundance 

        sampled <- unlist(lapply(classes_to_downsample, function(cd){
          cell_type_cells <- sce$Barcode[colData(sce)[[opt$cell_type_field]] == cd ]
          cell_type_cells[! cell_type_cells %in% sample(cell_type_cells, remove_freqs[[cd]])]
        }))

        selected_barcodes <- c(unsampled, sampled)
        sce <- sce[,sce$Barcode %in% selected_barcodes]
    }
    else{
        print("... unlabelled removed (where applicable), no further downsampling required")
    }
}
if(minor_cell_types_present){
    cell_type_freqs <- sort(table(colData(sce)[[opt$cell_type_field]]), decreasing = FALSE)
    cell_types_to_drop = cell_type_freqs[cell_type_freqs < cell_count_threshold]
    # need to make sure those rare cells are still present after down-sampling
    if(length(cell_types_to_drop) > 0){
        print("Removing cell types with low frequency...")
        barcodes_to_drop <- colData(sce)[[opt$cell_type_field]] %in% cell_types_to_drop
        sce <- sce[, !barcodes_to_drop]
    }
}

cell_type_freqs <- sort(table(colData(sce)[[opt$cell_type_field]]), decreasing = TRUE)
print("Ending cell type frequencies:")
print(cell_type_freqs)
    
print(paste('Final object has', ncol(sce), 'cells'))

# write data
print("Writing outputs")
write10xCounts(opt$output_dir, assays(sce)[[1]], barcodes = sce$Barcode, 
                                                 gene.id=rownames(sce))
# rename cell id field and remove redundant column
colnames(colData(sce))[1] <- opt$cell_id_field
colData(sce) <- colData(sce)[, -2] 
write.table(colData(sce), opt$metadata_upd, sep="\t")
print("Outputs written successfully")
