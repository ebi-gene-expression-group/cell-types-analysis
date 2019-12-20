#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(hash))
suppressPackageStartupMessages(require(foreach))
suppressPackageStartupMessages(require(doParallel))

### Generate a set of emprirical distributions for metrics defined in Utils.R 
### The script takes reference dataset as an input and shuffles it a specified number of times 
### for each iteration of shuffling, metrics are calculated; CDFs are produced after that  

option_list = list(
    make_option(
        c("-i", "--input-ref-file"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to file with reference cell types'
    ),
    make_option(
        c("-l", "--label-column-ref"),
        action = "store",
        default = 'cell_type',
        type = 'character',
        help = 'Name of the label column in reference file'
    ),
    make_option(
        c("-s", "--cell-ontology-col"),
        action = "store",
        default = 'CL_term',
        type = 'character',
        help = 'Name of CL id column in reference dataset'
    ),
    make_option(
        c("-n", "--num-iterations"),
        action = "store",
        default = 5,
        type = 'numeric',
        help = 'Number of sampling iterations to construct empirical distribution'
    ),
    make_option(
        c("-c", "--num-cores"),
        action = "store",
        default = 2,
        type = 'numeric',
        help = 'Number of cores to run the process on. Default: 2'
    ),
    make_option(
        c("-g", "--ontology-graph"),
        action = "store",
        default = "data/cl-basic.obo",
        type = 'character',
        help = 'Path to the ontology graph in .obo or .xml format'
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
        help = 'Path to the output object in .rds format'
    )
)

opt = wsc_parse_args(option_list, mandatory = c("input_ref_file", "output_path", "ontology_graph"))
source("Utils.R")
reference_labs_df = read.csv(opt$input_ref_file, sep="\t")
reference_labs = as.character(reference_labs_df[, opt$label_column_ref])
num_iter = opt$num_iterations
ref_CL_terms = as.character(reference_labs_df[, opt$cell_ontology_col])
ontology = opt$ontology_graph
sim_metric = opt$semantic_sim_metric

unlabelled = c("unassigned", "Unassigned", "unknown",
                'Unknown','rand','Node','ambiguous', NA)
trivial_terms = c("cell", "of", "tissue", "entity", "type") # add common words here

# generate empirical distribution
.run_simulations = function(siml_num){
    print(paste("Running simulation ", siml_num, "out of ", num_iter))
    predicted_labs = sample(reference_labs)
    exact_match_prop = get_exact_matches(reference_labs, predicted_labs)
    mean_shared_terms = get_shared_terms_prop(reference_labs, predicted_labs, trivial_terms)
    metrics = get_f1(reference_labs, predicted_labs, unlabelled)
    median_F1 = metrics$MedF1
    accuracy = metrics$Acc
    siml = get_CL_similarity(reference_labs, ref_CL_terms, predicted_labs, ontology, sim_metric, unlabelled)

    metric_list = list(Exact_match_prop = exact_match_prop,
                       Mean_partial_match = mean_shared_terms,
                       Med_F1 = median_F1,
                       Accuracy = accuracy,
                       CL_similarity = siml)
    return(metric_list)
}

# run simulations 
n_cores = opt$num_cores
registerDoParallel(n_cores)
emp_samples = foreach (iter=1:num_iter) %dopar% {
  .run_simulations(iter)
}

emp_samples = do.call(rbind, emp_samples)
print(emp_samples)
# generate cumulative empirical distributions
emp_dist = apply(emp_samples, 2, function(col) ecdf(as.numeric(col)))
saveRDS(emp_dist, file = opt$output_path)
