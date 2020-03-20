#!/usr/bin/env Rscript 

# Combine output tables from classifiers trained on multiple data sets
# using the same cell type prediction tool. Select top n labels for further
# downstream analysis. 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

option_list = list(
    make_option(
        c("-i", "--input-dir"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the directory with standardised output .tsv files from multiple
                classifiers. It is expected that input files follow the format: A_B_final-labs.tsv,
                where A is dataset or origin and B is classifier used to obtain predictions.'
    ),
    make_option(
        c("-n", "--top-labels-num"),
        action = "store",
        default = 3,
        type = 'numeric',
        help = 'Number of top labels to keep'
    ),
    make_option(
        c("-e", "--exclusions"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Path to the yaml file with excluded terms. Must contain fields 'unlabelled' and 'trivial_terms'"
    ),
    make_option(
        c("-s", "--scores"),
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = 'Boolean: Are prediction scores available for the given method? Default: FALSE'
    ),
     make_option(
        c("-o", "--output-table"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the output table in text format. The following naming is expected: tool-X_combined.tsv'
    )
)

# parse arguments 
opt = wsc_parse_args(option_list, mandatory = c("input_dir", "output_table"))
# source function definitions 
script_dir = dirname(strsplit(commandArgs()[grep('--file=', commandArgs())], '=')[[1]][2])
source(file.path(script_dir, 'cell_types_utils.R'))
# import the rest of dependencies 
suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(yaml))

# parse input tables
file_names = list.files(opt$input_dir, full.names=TRUE)
predicted_labs_tables = lapply(file_names, function(file) read.csv(file, sep="\t", stringsAsFactors=FALSE))
# extract datasets of origin
file_names = sapply(file_names, function(p) file_path_sans_ext(basename(p)))
datasets = sapply(file_names, function(name) unlist(strsplit(name, "_"))[1])

# check cell ids are identical across tables
cell_ids = get_unq_cell_ids(predicted_labs_tables)
# extract labels as lists of vectors and transform them into 2d arrays 
labels = lapply(predicted_labs_tables, function(tbl) tolower(tbl[, "predicted_label"]))

# read in exclusions config file if provided
if(!is.na(opt$exclusions)){
    e = yaml.load_file(opt$exclusions)
    unlabelled = tolower(e$unlabelled)
}

# coerce unlabelled cells to NAs
.filter_labs = function(lab_vec){
    lab_vec[lab_vec %in% unlabelled] = NA
    return(lab_vec)
}

labels = lapply(labels, .filter_labs)
labels = do.call(cbind, labels)
top_n = opt$top_labels_num
if(opt$scores){ 
    # select top n predictions based on prediction scores
    scores = lapply(predicted_labs_tables, function(tbl) tbl[, "score"])
    scores = do.call(cbind, scores)

    .get_top_scores = function(score_list, top_n){
        sort_lst = sort(score_list, index.return=TRUE, na.last=TRUE)
        # restrict to top_n items
        top_items = lapply(sort_lst, function(v) v[1:top_n])
        return(top_items)
    }

    # get array of top n scores and corresponding column indices for each cell
    top_items = apply(scores, 1, function(score_list) .get_top_scores(score_list, top_n))

    # extract list of labels corresponding to top score indices
    # top_items contains both sorted elements and corresponding indices
    # select columns of the label array and corresponding datasets by this index 
    sorted_lab_idx = lapply(seq_along(top_items), function(idx) as.numeric(unlist(top_items[[idx]][2])))
    top_labels = lapply(seq_along(nrow(labels)), function(idx) labels[ idx, sorted_lab_idx[[idx]] ])
    top_labels = data.frame(do.call(rbind, top_labels))
    names = paste("label", c(1:top_n), sep="_")
    colnames(top_labels) = names

    # determine which datasets produced top predicitons for each cell
    top_datasets = lapply(sorted_lab_idx, function(idx_lst) datasets[idx_lst])
    top_datasets = data.frame(do.call(rbind, top_datasets))
    names = paste("dataset", c(1:top_n), sep="_")
    colnames(top_datasets) = names

    # extract list of top scores per cell
    top_scores = lapply(top_items, function(l) l[[1]])
    top_scores = data.frame(do.call(rbind, top_scores))
    names = paste("score", c(1:top_n), sep="_")
    colnames(top_scores) = names

    # combine output columns
    output_table = cbind(cell_id = cell_ids, top_labels, top_datasets, top_scores)
    write.table(output_table, opt$output_table, sep="\t", row.names=FALSE)
} else{
    # scores are not available, therefore simply return all labels 
    labels = data.frame(labels)
    colnames(labels) = paste("label", c(1:ncol(labels)), sep="_")
    datasets = lapply(1:nrow(labels), function(idx) return(datasets))
    datasets = data.frame(do.call(rbind, datasets))
    colnames(datasets) = paste("dataset", c(1:ncol(labels)), sep="_")
    output_table = cbind(cell_id=cell_ids, labels, datasets)
    write.table(output_table, opt$output_table, sep="\t", row.names=FALSE)
}
