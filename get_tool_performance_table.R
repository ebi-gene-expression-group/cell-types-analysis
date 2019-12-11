#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

#### Create a table for evaluation metrics of multiple methods ####
#### Inputs: 
####    1) Text file with reference cell types
####    2) Directory path containing a list of output files from multiple methods.
####       A standard format is assumed: first column - cell_id; second column: pred_label
####


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
        c("-r", "--ref-file"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the file with reference cell type assignments'
    ),
    make_option(
        c("-s", "--cell-ontology-col"),
        action = "store",
        default = 'CL_term',
        type = 'character',
        help = 'Name of CL id column in ref dataset'
    ),
    make_option(
        c("-l", "--label-column-ref"),
        action = "store",
        default = 'cell_type',
        type = 'character',
        help = 'Name of the label column name in reference file'
    ),
    make_option(
        c("-p", "--label-column-pred"),
        action = "store",
        default = 'pred_label',
        type = 'character',
        help = 'Name of the label column name in prediction file'
    ),
    make_option(
        c("-o", "--output-path"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the output table in .tsv format'
    )
)

# source function definitions 
source("Utils.R")
# parse arguments 
opt = wsc_parse_args(option_list, mandatory = c("input_dir", "ref_file", "output_path"))

# file names must start with the tool name 
file_names = list.files(opt$input_dir, full.names=TRUE)
predicted_labs_tables = lapply(file_names, function(file) read.csv(file, sep="\t"))
pred_labs_col = opt$label_column_pred

# read reference file 
reference_labs_df = read.csv(opt$ref_file, sep=" ") #TODO: decide on separator
ref_labs_col = opt$label_column_ref
# NB: keep this relevant to the tools' output
unlabelled = c("unassigned", "Unassigned", "unknown")
# proportion of unknowns in reference cell types
prop_unlab_reference = get_unlab_rate(reference_labs_df[, ref_labs_col], unlabelled)

output_table = list()
# iterate through tools' outputs and calculate relevant statistics per tool
for(idx in 1:length(predicted_labs_tables)){
    predicted_labs_df = predicted_labs_tables[[idx]]
    tool = unlist(strsplit(basename(file_names[idx]), "_"))[1]
    print(tool)

    # check reference cell IDs match predicted cell IDs
    if(all(as.character(predicted_labs_df[, "cell_id"]) == as.character(predicted_labs_df[, "cell_id"]))) {
        # extract label vectors
        print("extract label vectors")
        predicted_labs = as.character(predicted_labs_df[, pred_labs_col])
        reference_labs = as.character(reference_labs_df[, ref_labs_col])
    } else {
        stop(paste("Error: cell id mismatch for tool: ", tool))
    }

    ############################################
                # obtain metrics
    ############################################
    # proportion of unlabelled cells in predicted labels
    prop_unlab_predicted = get_unlab_rate(predicted_labs, unlabelled)
    
    # ratio of unlabelled cells proportions
    unlab_ratio = get_unlab_rate_ratio(prop_unlab_reference, prop_unlab_predicted)

    # propotion of exact matches 
    exact_match_prop = get_exact_matches(reference_labs, predicted_labs)

    # match based on shared terms
    trivial_terms = c("cell", "of", "tissue") # add common words here
    mean_shared_terms = get_shared_terms_prop(reference_labs, predicted_labs, trivial_terms)

    # calculate F1 scores and accuracy
    metrics = get_preformance_metrics(reference_labs, predicted_labs, unlabelled)
    median_F1 = metrics$MedF1
    accuracy = metrics$Acc

    # cell ontology similarity
    CL_col = opt$cell_ontology_col 
    siml = get_CL_similarity(reference_labs_df, ref_labs_col, CL_col, predicted_labs)

    # combine metrics
    row = cbind(Tool = tool,
                Unlab_ref = prop_unlab_reference,
                Unlab_pred = prop_unlab_predicted,
                Unlab_ratio = unlab_ratio,
                Exact_match_prop = exact_match_prop,
                Mean_partial_match = mean_shared_terms,
                Med_F1 = median_F1,
                Accuracy = accuracy,
                CL_similarity = siml)
    output_table[[idx]] = row
}

output_table = data.frame(do.call(rbind, output_table))
write.table(output_table, file = opt$output_path, sep ="\t")

