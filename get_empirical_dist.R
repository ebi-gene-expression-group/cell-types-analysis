#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

### Generate a set of emprirical distributions for metrics defined for tool performance table
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
        c("-e", "--exclusions"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Path to the yaml file with excluded terms. Must contain fields 'unlabelled' and 'trivial_terms'"
    ),
    make_option(
        c("-l", "--label-column-ref"),
        action = "store",
        default = 'cell_type',
        type = 'character',
        help = 'Name of the label column in reference file'
    ),
    make_option(
        c("-m", "--lab-cl-mapping"),
        action = "store",
        default = NA, 
        type = 'character',
        help = 'Path to serialised object containing cell label to CL terms mapping'
    ),
    make_option(
        c("-p", "--parallel"),
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = 'Boolean: Should computation be run in parallel? Default: FALSE'
    ),
    make_option(
        c("-n", "--num-iterations"),
        action = "store",
        default = 5,
        type = 'integer',
        help = 'Number of sampling iterations to construct empirical distribution'
    ),
    make_option(
        c("-c", "--num-cores"),
        action = "store",
        default = NA,
        type = 'integer',
        help = 'Number of cores to run the process on. Default: all available cores. --parallel must be set to "true" for this to take effect'
    ),
    make_option(
        c("-g", "--ontology-graph"),
        action = "store",
        default = "data/cl-basic.obo",
        type = 'character',
        help = 'Path to the ontology graph in .obo or .xml format'
    ),
    make_option(
        c("-s", "--semantic-sim-metric"),
        action = "store",
        default = 'edge_resnik',
        type = 'character',
        help = 'Semantic similarity scoring method. Must be supported by Onassis package. See listSimilarities()$pairwiseMeasures for a list of accepted options'
    ),
    make_option(
        c("-o", "--output-path"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the output CDF list object in .rds format'
    )
)

opt = wsc_parse_args(option_list, mandatory = c("input_ref_file", "output_path", "lab_cl_mapping", "ontology_graph"))
script_dir = dirname(strsplit(commandArgs()[grep('--file=', commandArgs())], '=')[[1]][2])
source(file.path(script_dir, 'cell_types_utils.R'))
# import the rest of dependencies 
suppressPackageStartupMessages(require(hash))
suppressPackageStartupMessages(require(foreach))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(doParallel))
suppressPackageStartupMessages(require(yaml))

# read in exclusions file, if provided
if(! is.na(opt$exclusions)){
    e = yaml.load_file(opt$exclusions)
    unlabelled = tolower(e$unlabelled)
    trivial_terms = tolower(e$trivial_terms)
}

reference_labs_df = read.csv(opt$input_ref_file, sep="\t", stringsAsFactors=FALSE)
reference_labs = reference_labs_df[, opt$label_column_ref]
num_iter = opt$num_iterations
ontology = opt$ontology_graph
lab_cl_mapping = readRDS(opt$lab_cl_mapping)
sim_metric = opt$semantic_sim_metric

# generate empirical distribution by running simulations
.run_simulations = function(siml_num){
    print(paste("Running simulation ", siml_num, "out of ", num_iter))
    predicted_labs = sample(reference_labs)
    exact_match_prop = get_exact_matches(reference_labs, predicted_labs)
    mean_shared_terms = get_shared_terms_prop(reference_labs, predicted_labs, trivial_terms)
    metrics = get_f1(reference_labs, predicted_labs, unlabelled)
    median_F1 = metrics$MedF1
    accuracy = metrics$Acc
    sim_vec = c()
    for(idx in seq_along(reference_labs)){
        lab_1 = reference_labs[idx]
        lab_2 = predicted_labs[idx]
        siml = get_CL_similarity(lab_1, lab_2, 
                                         lab_cl_mapping=lab_cl_mapping, 
                                         ontology=ontology, 
                                         sim_metric=sim_metric, 
                                         unlabelled=unlabelled)
        sim_vec[idx] = log10(siml+1)
    }
    siml = mean(sim_vec, na.rm=TRUE)
    metric_list = list(Exact_match_prop = exact_match_prop,
                       Mean_partial_match = mean_shared_terms,
                       Med_F1 = median_F1,
                       Accuracy = accuracy,
                       CL_similarity = siml)

    return(metric_list)
}

# run simulations in parallel, if specified
if(opt$parallel){
    # run simulations 
    if(is.na(opt$num_cores)){
        n_cores = detectCores()
    } else {
        n_cores = opt$num_cores
    }
    registerDoParallel(n_cores)
    emp_samples = foreach (iter=1:num_iter) %dopar% {
      .run_simulations(iter)
    }
} else {
    # run simulations sequentially
    emp_samples = lapply(1:num_iter, function(idx) .run_simulations(idx))
}

emp_samples = do.call(rbind, emp_samples)
# generate cumulative empirical distributions
emp_dist = apply(emp_samples, 2, function(col) ecdf(as.numeric(col)))
saveRDS(emp_dist, file = opt$output_path)
