#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

#### Create a table for evaluation metrics of multiple methods ####

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
        c("-n", "--parallel"),
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = 'Boolean: Should computation be run in parallel? Default: FALSE'
    ),
     make_option(
        c("-c", "--num-cores"),
        action = "store",
        default = NA,
        type = 'double',
        help = 'Number of cores to run the process on. Default: all available cores. --parallel must be set to "true" for this to take effect'
    ),
     make_option(
        c("-e", "--exclusions"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Path to the yaml file with excluded terms. Must contain fields 'unlabelled' and 'trivial_terms'"
    ),
     make_option(
        c("-d", "--tmpdir"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Cache directory path'
    ),
    make_option(
        c("-f", "--ontology-graph"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the ontology graph in .obo or .xml format. Import link can also be provided.'
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
        default = 'lin',
        type = 'character',
        help = 'Semantic similarity scoring method. Must be supported by Onassis
                package. See listSimilarities()$pairwiseMeasures for a list
                of accepted options. 
                NB: if included in combined score calculation, make sure to select a metric with values in the [0;1] range.'
    ),
    make_option(
        c("-k", "--include-sem-siml"),
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = 'Should semantic similarity be included into combined score calculation? Default: FALSE.
                If setting to TRUE, note that this confines the options on semantic similarity metric
                to those with range in the [0;1] interval only.'
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
script_dir = dirname(strsplit(commandArgs()[grep('--file=', commandArgs())], '=')[[1]][2])
source(file.path(script_dir, 'cell_types_utils.R'))

# import the rest of dependencies 
suppressPackageStartupMessages(require(hash)) # must be version 2.2.6.1
suppressPackageStartupMessages(require(foreach))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(doParallel))
suppressPackageStartupMessages(require(yaml))


# file names must start with the tool name 
file_names = list.files(opt$input_dir, full.names=TRUE)
predicted_labs_tables = lapply(file_names, function(file) read.csv(file, sep="\t", stringsAsFactors=FALSE, comment.char = "#"))
# extract corresponding tools 
tools = as.character(sapply(file_names, function(file) extract_metadata(file)[['tool']]))
pred_labs_col = opt$label_column_pred
barcode_ref = opt$barcode_col_ref
barcode_pred = opt$barcode_col_pred
include_sem_siml = opt$include_sem_siml

# read reference file 
reference_labs_df = read.csv(opt$ref_file, sep="\t", stringsAsFactors=FALSE)
# make sure there are no duplicate rows in reference SDRF file 
reference_labs_df = reference_labs_df[which(!duplicated(reference_labs_df[, barcode_ref])), ]
ref_labs_col = opt$label_column_ref

# parameters for semantic similarity
ontology = import_ontology_graph(opt$tmpdir, opt$ontology_graph)
lab_cl_mapping = readRDS(opt$lab_cl_mapping)
sim_metric = opt$semantic_sim_metric 

# read in exclusions file, if provided
if(! is.na(opt$exclusions)){
    e = yaml.load_file(opt$exclusions)
    unlabelled = tolower(e$unlabelled)
    trivial_terms = tolower(e$trivial_terms)
}

# find proportion of unknowns in reference cell types
prop_unlab_reference = get_unlab_rate(reference_labs_df[, ref_labs_col], unlabelled)

# iterate through tools' outputs and calculate relevant statistics per tool
.get_metrics = function(idx){
    predicted_labs_df = predicted_labs_tables[[idx]]
    #tool = unlist(strsplit(basename(file_names[idx]), "_"))[1]
    tool = tools[idx]
    print(paste("Evaluating tool:", tool, sep = " "))

    # check reference cell IDs match predicted cell IDs
    # if so, extract reference and predicted labels as vectors
    if(!all(as.character(reference_labs_df[, barcode_ref]) == as.character(predicted_labs_df[, barcode_pred]))) {
        stop(paste("Error: cell id mismatch for tool: ", tool))        
    }
    # extract label vectors
    predicted_labs = tolower(as.character(predicted_labs_df[, pred_labs_col]))
    reference_labs = tolower(as.character(reference_labs_df[, ref_labs_col]))
    # run evaluation function for a tool
    row = obtain_metrics_list(tool=tool, 
                              reference_labs=reference_labs,
                              predicted_labs=predicted_labs, 
                              prop_unlab_reference=prop_unlab_reference,
                              lab_cl_mapping=lab_cl_mapping,
                              unlabelled=unlabelled, 
                              trivial_terms=trivial_terms,
                              ontology=ontology,
                              sim_metric=sim_metric,
                              include_sem_siml=include_sem_siml)
    return(do.call(cbind, row))
}

n_tools = length(predicted_labs_tables)
# run parallel computation, if specified 
if(opt$parallel){
    # set number of cores for computation
    if(is.na(opt$num_cores)){
        n_cores = detectCores()
    } else {
        n_cores = opt$num_cores
    }
    registerDoParallel(n_cores)
    metrics_lst = foreach (iter=1:n_tools) %dopar% {
        .get_metrics(iter)
    }
} else{
    # run computations sequentially 
    metrics_lst = lapply(1:n_tools, function(idx) .get_metrics(idx))
}

output_table = data.frame(do.call(rbind, metrics_lst))
output_table = output_table[order(output_table$Combined_score), ]
write.table(output_table, file = opt$output_path, sep ="\t", row.names=FALSE)
