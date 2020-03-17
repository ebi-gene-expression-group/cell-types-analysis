#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(hash)) # must be version 2.2.6.1
suppressPackageStartupMessages(require(foreach))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(doParallel))

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
        help = 'Path to the file with reference, "true" cell type assignments'
    ),
     make_option(
        c("-c", "--num-cores"),
        action = "store",
        default = NA,
        type = 'double',
        help = 'Number of cores to run the process on. Default: all available cores'
    ),
    make_option(
        c("-f", "--ontology-graph"),
        action = "store",
        default = "data/cl-basic.obo",
        type = 'character',
        help = 'Path to the ontology graph in .obo or .xml format'
    ),
    make_option(
        c("-m", "--lab-cl-mapping"),
        action = "store",
        default = NA, 
        type = 'character',
        help = 'Path to serialised object containing cell label - CL terms mapping'
    ),
    make_option(
        c("-b", "--barcode-col-ref"),
        action = "store",
        default = 'cell_id',
        type = 'character',
        help = 'Name of the cell id field in reference file'
    ),
    make_option(
        c("-a", "--barcode-col-pred"),
        action = "store",
        default = 'cell_id',
        type = 'character',
        help = 'Name of the cell id field in prediction file'
    ),
    make_option(
        c("-l", "--label-column-ref"),
        action = "store",
        default = 'cell_type',
        type = 'character',
        help = 'Name of the label column in reference file'
    ),
    make_option(
        c("-p", "--label-column-pred"),
        action = "store",
        default = 'pred_label',
        type = 'character',
        help = 'Name of the label column in prediction file'
    ),
    make_option(
        c("-s", "--semantic-sim-metric"),
        action = "store",
        default = 'edge_resnik',
        type = 'character',
        help = 'Semantic similarity scoring method. Must be supported by Onassis
                package. See listSimilarities()$pairwiseMeasures for a list
                of accepted options'
    ),
    make_option(
        c("-o", "--output-path"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the output table in .tsv format'
    )
)

# parse arguments 
opt = wsc_parse_args(option_list, mandatory = c("input_dir", "ref_file", "lab_cl_mapping", "output_path"))
# source function definitions 
p = system("which cell_types_utils.R", intern = TRUE)
source(p)

# file names must start with the tool name 
file_names = list.files(opt$input_dir, full.names=TRUE)
predicted_labs_tables = lapply(file_names, function(file) read.csv(file, sep="\t", stringsAsFactors=FALSE))
pred_labs_col = opt$label_column_pred
barcode_ref = opt$barcode_col_ref
barcode_pred = opt$barcode_col_pred

# read reference file 
reference_labs_df = read.csv(opt$ref_file, sep="\t", stringsAsFactors=FALSE)
# make sure there are no duplicate rows in reference SDRF file 
reference_labs_df = reference_labs_df[which(!duplicated(reference_labs_df[, barcode_ref])), ]
ref_labs_col = opt$label_column_ref

# parameters for semantic similarity 
ontology = opt$ontology_graph
lab_cl_mapping = readRDS(opt$lab_cl_mapping)
sim_metric = opt$semantic_sim_metric 

# find proportion of unknowns in reference cell types
prop_unlab_reference = get_unlab_rate(reference_labs_df[, ref_labs_col], unlabelled)

# iterate through tools' outputs and calculate relevant statistics per tool
.get_metrics = function(idx){
    predicted_labs_df = predicted_labs_tables[[idx]]
    tool = unlist(strsplit(basename(file_names[idx]), "_"))[1]
    print(paste("Evaluating tool:", tool, sep = " "))

    # check reference cell IDs match predicted cell IDs
    # if so, extract reference and predicted labels as vectors
    if(!all(as.character(reference_labs_df[, barcode_ref]) == as.character(predicted_labs_df[, barcode_pred]))) {
        stop(paste("Error: cell id mismatch for tool: ", tool))        
    }
    # extract label vectors
    predicted_labs = as.character(predicted_labs_df[, pred_labs_col])
    reference_labs = as.character(reference_labs_df[, ref_labs_col])
    # run evaluation function for a tool
    row = obtain_metrics_list(tool=tool, 
                              reference_labs=reference_labs,
                              predicted_labs=predicted_labs, 
                              prop_unlab_reference=prop_unlab_reference,
                              lab_cl_mapping=lab_cl_mapping,
                              unlabelled=unlabelled, 
                              trivial_terms=trivial_terms,
                              ontology=ontology,
                              sim_metric=sim_metric)
    return(do.call(cbind, row))
}

# set number of cores for computation
if(is.na(opt$num_cores)){
    n_cores = detectCores()
} else {
    n_cores = opt$num_cores
}
n_tools = length(predicted_labs_tables)
registerDoParallel(n_cores)
metrics_lst = foreach (iter=1:n_tools) %dopar% {
    .get_metrics(iter)
}
output_table = data.frame(do.call(rbind, metrics_lst))
output_table = output_table[order(output_table$Combined_score), ]
write.table(output_table, file = opt$output_path, sep ="\t", row.names=FALSE)
