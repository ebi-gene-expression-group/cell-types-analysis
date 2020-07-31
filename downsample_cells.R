#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# When training classifiers, avoid memory overflow by weighted down-sampling of cells.
# The most prevalent cell types are filtered out first, so less-represented cells are not under risk of being removed. 

option_list = list(
    make_option(
        c("-d", "--expression-data"),
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
        c("-f", "--cell-id-field"),
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
        c("-s", "--sampling-rate"),
        action = "store",
        default = 0.25,
        type = 'numeric',
        help = 'Percantage of cells to be removed in a single iteration'
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
suppressPackageStartupMessages(require(Matrix))

expr_data = opt$expression_data
genes = read.csv(paste(expr_data, "genes.tsv", sep="/"), sep="\t", stringsAsFactors = FALSE, header = FALSE)
barcodes = read.csv(paste(expr_data, "barcodes.tsv", sep="/"), sep="\t", stringsAsFactors = FALSE, header = FALSE)
cell_num_limit = floor(opt$array_size_limit / nrow(genes))
current_cell_num = nrow(barcodes)

# if no down-samling is required, simply return the output
if(current_cell_num <= cell_num_limit){
    system(paste("mv", opt$expression_data, opt$output_dir, sep=" "))
    system(paste("mv", opt$metadata, opt$metadata_upd, sep=" "))
    cat(paste("Matrix not large enough, down-sampling is not required.\nData are moved to ",
               opt$output_dir, " and ", opt$metadata_upd, "\n", sep=""))
    quit(status = 0)
}

# parse the remaining data
cell_labels = opt$cell_type_field
cell_id_col = opt$cell_id_field
matrix = Matrix::readMM(paste(expr_data, "matrix.mtx", sep="/"))
metadata = read.csv(opt$metadata, sep = "\t", stringsAsFactors = FALSE)

# remove technical duplicate rows
metadata = metadata[which(!duplicated(metadata[, cell_id_col])), ]
# get indices of all cell types; starting with most abundant
num_per_cell_type = sort(table(metadata[, cell_labels]), decreasing = TRUE)
grouped_indices = lapply(seq_along(num_per_cell_type),
                   function(idx) which(metadata[, cell_labels] == names(num_per_cell_type[idx])))

# build a vector of cell indices to be removed 
total_cells_to_remove = c()
current_cell_type_idx = 1

while(current_cell_num >= cell_num_limit){
    sample_size = floor(length(grouped_indices[[current_cell_type_idx]]) * opt$sampling_rate)
    # sample by index, then select cells by that index
    i = sample(1:length(grouped_indices[[current_cell_type_idx]]), sample_size)
    current_cells_to_remove = grouped_indices[[current_cell_type_idx]][i]
    total_cells_to_remove = append(total_cells_to_remove, current_cells_to_remove)

    # remove sampled cells 
    grouped_indices[[current_cell_type_idx]] = grouped_indices[[current_cell_type_idx]][-i]
    # continue sampling from the same cell type or switch to next one    
    if(current_cell_type_idx < length(num_per_cell_type)){
        if(length(grouped_indices[[current_cell_type_idx]]) <= length(grouped_indices[[current_cell_type_idx + 1]])){
            current_cell_type_idx = current_cell_type_idx + 1
        }
    } else{
        current_cell_type_idx = 1
    }
    current_cell_num = current_cell_num - sample_size
}

# subset matrix, barcodes and metadata by the generated index
matrix = matrix[, -total_cells_to_remove]
barcodes = data.frame(barcodes[-total_cells_to_remove, ])
metadata = metadata[-total_cells_to_remove, ]

# write data
dir.create(opt$output_dir)
Matrix::writeMM(matrix, paste(opt$output_dir, "matrix.mtx", sep="/"))
write.table(genes, paste(opt$output_dir, "genes.tsv", sep="/"), sep="\t", col.names = FALSE)
write.table(barcodes, paste(opt$output_dir, "barcodes.tsv", sep="/"), sep="\t", row.names=FALSE, col.names = FALSE)
write.table(metadata, opt$metadata_upd, sep="\t", row.names = FALSE)
