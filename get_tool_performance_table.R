#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(hash)) # NB: must be version 2.2.6.1

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
        c("-f", "--ontology-graph"),
        action = "store",
        default = "data/cl-basic.obo",
        type = 'character',
        help = 'Path to the ontology graph in .obo or .xml format'
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
        c("-m", "--semantic-sim-metric"),
        action = "store",
        default = 'edge_resnik',
        type = 'character',
        help = 'Semantic similarity scoring method. 
                Must be supported by Onassis package.
                See listSimilarities()$pairwiseMeasures for a list of accepted options'
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
opt = wsc_parse_args(option_list, mandatory = c("input_dir", "ref_file", "output_path"))
# source function definitions 
source("Utils.R")

# file names must start with the tool name 
file_names = list.files(opt$input_dir, full.names=TRUE)
predicted_labs_tables = lapply(file_names, function(file) read.csv(file, sep="\t"))
pred_labs_col = opt$label_column_pred

# read reference file 
reference_labs_df = read.csv(opt$ref_file, sep="\t")
ref_labs_col = opt$label_column_ref

# NB: keep these relevant to the tools' output
unlabelled = c("unassigned", "Unassigned", "unknown",
                'Unknown','rand','Node','ambiguous', NA)
trivial_terms = c("cell", "of", "tissue", "entity", "type") # add common words here

# extract ontology terms for cell types 
ontology = opt$ontology_graph
CL_col = opt$cell_ontology_col
ref_CL_terms = as.character(reference_labs_df[, CL_col])
sim_metric = opt$semantic_sim_metric 

# find proportion of unknowns in reference cell types
prop_unlab_reference = get_unlab_rate(reference_labs_df[, ref_labs_col], unlabelled)

# iterate through tools' outputs and calculate relevant statistics per tool
output_table = list()
for(idx in 1:length(predicted_labs_tables)){
    predicted_labs_df = predicted_labs_tables[[idx]]
    tool = unlist(strsplit(basename(file_names[idx]), "_"))[1]
    print(paste("Evaluating tool:", tool, sep = " "))

    # check reference cell IDs match predicted cell IDs
    # if so, extract reference and predicted labels as vectors
    if(!all(as.character(predicted_labs_df[, "cell_id"]) == as.character(predicted_labs_df[, "cell_id"]))) {
        stop(paste("Error: cell id mismatch for tool: ", tool))        
    }
    # extract label vectors
    predicted_labs = as.character(predicted_labs_df[, pred_labs_col])
    reference_labs = as.character(reference_labs_df[, ref_labs_col])
    # run evaluation functions
    row = obtain_metrics_list(tool, reference_labs,
                              predicted_labs, prop_unlab_reference,
                              unlabelled, trivial_terms, ref_CL_terms, ontology, sim_metric)
    row = do.call(cbind, row)
    output_table[[idx]] = row
}

output_table = data.frame(do.call(rbind, output_table))
output_table = output_table[order(output_table$Combined_score), ]
write.table(output_table, file = opt$output_path, sep ="\t")
