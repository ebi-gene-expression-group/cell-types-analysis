#!/usr/bin/env Rscript 

# This script accepts a collection of (condensed) SDRF files as input and returns a 
# dictionary containing cell labels as keys and ontology terms as values 
# for all cells in corresponding SDRF files. In case a cell does not have 
# ontology mapping, label will map to NA value. 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(tools))

option_list = list(
    make_option(
        c("-i", "--input-dir"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the directory with condensed SDRF files'
    ),
    make_option(
        c("-k", "--condensed-sdrf"),
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = 'Boolean: is the provided SDRF file in a condensed form? Default: TRUE'
    ),
    make_option(
        c("-b", "--barcode-col-name"),
        action = "store",
        default = "id",
        type = 'character',
        help = 'Name of the barcode column in SDRF files (must be identical across all files)'
    ),
    make_option(
        c("-l", "--cell-label-col-name"),
        action = "store",
        default = "inferred cell type",
        type = 'character',
        help = 'Name of the cell label column in SDRF files (must be identical across all files)'
    ),
    make_option(
        c("-c", "--cell-ontology-col-name"),
        action = "store",
        default = "cell.type.ontology",
        type = 'character',
        help = 'Name of the cell ontology terms column in SDRF files (must be identical across all files)'
    ),
    make_option(
        c("-o", "--output-dict-path"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path for serialised object containing the dictionary'
    ),
    make_option(
        c("-t", "--output-text-path"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path for txt version of label - term mapping'
    )
)

# parse arguments 
opt = wsc_parse_args(option_list, mandatory = c("input_dir", "output_dict_path", "output_text_path"))
# source function definitions 
script_dir = dirname(strsplit(commandArgs()[grep('--file=', commandArgs())], '=')[[1]][2])
source(file.path(script_dir, 'cell_types_utils.R'))
cell_label = opt$cell_label_col_name

# import the rest of dependencies 
suppressPackageStartupMessages(require(hash))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))
# parse input SDRF files
condensed = opt$condensed_sdrf
file_names = list.files(opt$input_dir, full.names=TRUE)
if(condensed){
    sdrf_tables = lapply(file_names, function(file) data.frame(fread(file, header=FALSE, stringsAsFactors = FALSE, fill = TRUE)))
} else{
    sdrf_tables = lapply(file_names, function(file) read.csv(file, sep="\t", stringsAsFactors = FALSE))
}

# extract cell labels and CL terms from SDRF files (condensed or non-condensed)
# then map each label to corresponding CL term in a hash table 
.map_cell_labels = function(df, hash_table, condensed){
    if(condensed){
        # select rows which have cell type 
        df = df[df[, 5] == cell_label, ]
        cell_labs = tolower(df[, 6])
        cl_terms = df[, 7]
    } else {
        cell_labs = tolower(df[, cell_label])
        cl_terms = df[, opt$cell_ontology_col_name]
    }
    
    for(idx in seq_along(cell_labs)){
        # make sure CL terms not duplicated
        if(grepl(",", cl_terms[idx])){
            cl_terms[idx] = unlist(strsplit(cl_terms[idx], ","))[1]
        }
        # add key-value pair to dict 
        if(!(is.na(cell_labs[idx]) | is.na(cl_terms[idx]) | cell_labs[idx]=="" | cl_terms[idx]=="")){
            hash_table[cell_labs[idx]] = cl_terms[idx]
        }
    }
}

# populate dictionary with key-value pairs 
cell_type_id_mapping = hash()
for(df in sdrf_tables){
    .map_cell_labels(df, cell_type_id_mapping, condensed)
}
out_path = opt$output_dict_path
saveRDS(cell_type_id_mapping, out_path)
# save a human-readable version of the dictionary 
label_cl_table = data.frame(cell_labels=keys(cell_type_id_mapping), CL_terms=values(cell_type_id_mapping))
write.table(label_cl_table, file=opt$output_text_path, sep="\t", row.names=FALSE)
