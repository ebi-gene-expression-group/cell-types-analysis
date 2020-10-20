#!/usr/bin/env Rscript

# Remove non-aplhanumeric characters 
# set to lower-case, if scpecified

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(data.table))

option_list = list(
    make_option(
        c("-i", "--input-file"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to input metadata file in .tsv format'
    ),
    make_option(
        c("-l", "--label-field"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Name of label field in metadata file'
    ),
    make_option(
        c("-k", "--condensed"),
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = 'Is the provided metadata file in condensed format? Default: False'
    ),
    make_option(
        c("-t", "--attribute-type-col-num"),
        action = "store",
        default = 5,
        type = 'numeric',
        help = 'Number of the attribute type field in condensed metadata file. Default: 5'
    ),
    make_option(
        c("-v", "--variable-col-num"),
        action = "store",
        default = 6,
        type = 'numeric',
        help = 'Number of the label field in condensed metadata file. Default: 6'
    ),
    make_option(
        c("-a", "--avoid-lowercase"),
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = 'Should setting the labels to lowercase be skipped? Default: False'
    ),
    make_option(
        c("-o", "--output-path"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output for updated file'
    )
)

# parse arguments 
opt = wsc_parse_args(option_list, mandatory = c("input_file", "label_field", "output_path"))
# source function definitions 
script_dir = dirname(strsplit(commandArgs()[grep('--file=', commandArgs())], '=')[[1]][2])
source(file.path(script_dir, 'cell_types_utils.R'))

if(!file.exists(opt$input_file)){
    stop("Input file is not provided")
}

reg_exp = "[^-A-Za-z0-9>+ ]"
cond = opt$condensed
lab_field = opt$label_field

if(cond){
    data = data.frame(fread(opt$input_file, header=FALSE, stringsAsFactors = FALSE, fill = TRUE))
    # extract rows with labels 
    idx = grep(lab_field, data[, opt$attribute_type_col_num])
    if(length(idx) < 1){
        stop("No labels found in supplied file")
    }
    labels = data[idx, opt$variable_col_num]
} else{
    data = read.csv(opt$input_file, sep = "\t", stringsAsFactors = FALSE)
    if(!lab_field %in% colnames(data)){
    stop("Provided label field not found in metada file")
    }
    labels = data[, lab_field]
}

# filter labels
labels = filter_labels(labels, reg_exp)
if(! opt$avoid_lowercase){
    labels = tolower(labels)
}

# update metadata file
if(cond){
    data[idx, 6] = labels
    cols = FALSE
    
} else{
    data$lab_field = labels
    cols = TRUE
}
write.table(data, opt$output_path, sep = "\t", row.names=FALSE, col.names = cols)
