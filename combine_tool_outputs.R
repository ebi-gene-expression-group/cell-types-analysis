#!/usr/bin/env Rscript 

# combine output tables from multiple tools and select top n labels for further analysis
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

option_list = list(
    make_option(
        c("-i", "--input-dir"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the directory with standardised output .tsv files from multiple
                methods'
    ),
    make_option(
        c("-n", "--top-labels-num"),
        action = "store",
        default = 3,
        type = 'numeric',
        help = 'Number of top labels to keep'
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
        help = 'Path to the output table in text format'
    )
)

# parse arguments 
opt = wsc_parse_args(option_list, mandatory = c("input_dir", "output_table"))
# source function definitions 
#p = system("which cell_types_utils.R", intern = TRUE)
#source(p)

# parse input tables
file_names = list.files(opt$input_dir, full.names=TRUE)
predicted_labs_tables = lapply(file_names, function(file) read.csv(file, sep="\t", stringsAsFactors=FALSE))
top_n = opt$top_labels_num

# check cell ids are identical across tables
cell_ids = lapply(predicted_labs_tables, function(tbl) tbl[ , "cell_id"])
if(! length(unique(cell_ids)) == 1){
    stop("Inconsistent cell ids provided")
}
cell_ids = unlist(unique(cell_ids))
# extract labels as lists of vectors and transform them into 2d arrays 
labels = lapply(predicted_labs_tables, function(tbl) tbl[, "predicted_label"])
labels = do.call(cbind, labels)

if(opt$scores){
    # select top n predictions based on scores

    scores = lapply(predicted_labs_tables, function(tbl) tbl[, "score"])
    scores = do.call(cbind, scores)

    .get_top_scores = function(score_list, top_n){
        sort_lst = sort(score_list, index.return=TRUE, na.last=TRUE)
        # select only top n items
        top_items = lapply(sort_lst, function(v) v[1:top_n])
        return(top_items)
    }

    # get 3d array of top n scores and corresponding indices
    top_items = apply(scores, 1, function(score_list) .get_top_scores(score_list, top_n))

    # extract list of labels corresponding to top score indices
    top_labels = lapply(seq_along(top_items), function(idx) labels[ idx, as.numeric(unlist(top_items[[idx]][2])) ])
    top_labels = data.frame(do.call(rbind, top_labels))
    names = paste("label", c(1:top_n), sep="_")
    colnames(top_labels) = names

    # extract list of top scores per cell
    top_scores = lapply(top_items, function(l) l[[1]])
    top_scores = data.frame(do.call(rbind, top_scores))
    names = paste("score", c(1:top_n), sep="_")
    colnames(top_scores) = names

    # combine output columns
    output_table = cbind(cell_id = cell_ids, top_labels, top_scores)
    write.table(output_table, opt$output_table, sep="\t", row.names=FALSE)
} else{
    # in case a tool does not allow to retrieve prediction scores, 
    # filter out unlabelled predictions
    unlabelled = c("unassigned", "Unassigned", "unknown",
                'Unknown','rand','Node','ambiguous', NA)

    labels = apply(labels, 1, function(row) row[which(!row %in% unlabelled)])
    if(length(labels) == 0){ 
        stop("There are no labelled cells in the data set")
    }
    max_lab = max(unlist(sapply(labels, length)))

    .fill_vectors = function(vec, max_len){
        d = max_len - length(vec)
        print(d)
        return(append(vec, rep(NA, d)))
    }
    # fill empty spots with NAs 
    labels = lapply(labels, function(lab_vec) .fill_vectors(lab_vec, max_lab))
    labels = data.frame(do.call(rbind, labels))
    names = paste("label", c(1:max_lab), sep="_")
    colnames(labels) = names
    output_table = cbind(cell_id=cell_ids, labels)
    write.table(output_table, opt$output_table, sep="\t", row.names=FALSE)
}
