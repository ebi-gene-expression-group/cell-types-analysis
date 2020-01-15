#!/usr/bin/env Rscript 
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(hash))
suppressPackageStartupMessages(require(Onassis))
suppressPackageStartupMessages(require(matrixStats))

### Concatenate labels for specific cell obtained from different tools.
### Calculate statistics: agreement between tools, CL similarity among predictions,  
### and provide aggregate scores across predicted labels 

option_list = list(
    make_option(
        c("-i", "--input-dir"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the directory with standardised .tsv files from multiple
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
        c("-t", "--tool-table"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the tool evaluation table in text format'
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
        help = 'Name of CL id column in reference dataset'
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

# parse arguments 
opt = wsc_parse_args(option_list, mandatory = c("input_dir", "ref_file", "output_path"))
# source function definitions 
p = system("which cell_types_utils.R", intern = TRUE)
source(p)

# file names must start with the tool name 
file_names = list.files(opt$input_dir, full.names=TRUE)
predicted_labs_tables = lapply(file_names, function(file) read.csv(file, sep="\t"))
pred_labs_col = opt$label_column_pred
barcode_ref = opt$barcode_col_ref
barcode_pred = opt$barcode_col_pred

# read reference data 
reference_labs_df = read.csv(opt$ref_file, sep="\t")
ref_labs_col = opt$label_column_ref
ref_CL_terms = as.character(reference_labs_df[, opt$cell_ontology_col]) 
reference_labs = as.character(reference_labs_df[, ref_labs_col])
output_table = reference_labs_df[, c(barcode_ref, ref_labs_col)]
# NB: keep these relevant to the tools' output
unlabelled = c("unassigned", "Unassigned", "unknown",
                'Unknown','rand','Node','ambiguous', NA)

# iterate through tools' outputs and combine predictions per cell
tools = c() 
for(idx in 1:length(predicted_labs_tables)){
    predicted_labs_df = predicted_labs_tables[[idx]]
    tool = unlist(strsplit(basename(file_names[idx]), "_"))[1]
    tools[idx] = tool

    # check reference cell IDs match predicted cell IDs
    # if so, extract reference and predicted labels as vectors
    if(!all(as.character(reference_labs_df[, barcode_ref]) == as.character(predicted_labs_df[, barcode_pred]))) {
         stop(paste("Error: cell id mismatch for tool: ", tool))
    }    
    # extract label vectors
    predicted_labs = as.character(predicted_labs_df[, pred_labs_col])
    output_table = cbind(output_table, predicted_labs)
}

# add tool names to columns 
names(output_table)[3:length(names(output_table))] = tools
# find agreement between tools for a specific cell 
agreement_rate = apply(output_table, 1, get_agreement_rate)

# find semantic similarity among labels
siml_object = new('Similarity')
ontology = opt$ontology_graph
ontology(siml_object) = ontology
metric = opt$semantic_sim_metric
# configure similarity measurement metric
idx = which(listSimilarities()$pairwiseMeasures == metric)
pairwiseConfig(siml_object) = listSimilarities()$pairwiseMeasures[idx]
cell_type_id_mapping = hash()
i = which(!(reference_labs %in% unlabelled | ref_CL_terms=='')) # filter unmapped cells 
.set(cell_type_id_mapping, keys=reference_labs[i], values=ref_CL_terms[i])
group_siml = apply(output_table, 1, function(row) get_cell_CL_siml(siml_object, row, cell_type_id_mapping))
output_table = cbind(output_table, agreement_rate = agreement_rate, group_siml = group_siml)

# obtain scores according to combined tool performance
tool_table = read.delim(opt$tool_table)
tool_scores = tool_table$Combined_score
names(tool_scores) = tool_table$Tool
tool_labels = output_table[, tools]
tool_labels = apply(tool_labels, 1, function(row) get_weighted_score(row, tool_scores))
output_table[, tools] = t(tool_labels)

# put cells with weakest evidence on top of the list 
output_table = output_table[order(agreement_rate, group_siml) ,]
write.table(output_table, file=opt$output_path, sep="\t")
