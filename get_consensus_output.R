#!/usr/bin/env Rscript 
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

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
        c("-t", "--tool-table"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the tool evaluation table in text format'
    ),
    make_option(
        c("-p", "--parallel"),
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = 'Boolean: Should computation be run in parallel? Default: FALSE'
    ),
     make_option(
        c("-c", "--num-cores"),
        action = "store",
        default = NA,
        type = 'integer',
        help = 'Number of cores to run the process on. Default: all available cores. --parallel must be set to "true" for this to take effect'
    ),
    make_option(
        c("-d", "--cl-dictionary"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the mapping between labels and CL terms in .rds format'
    ),
     make_option(
        c("-e", "--exclusions"),
        action = "store",
        default = NA,
        type = 'character',
        help = "Path to the yaml file with excluded terms. Must contain fields 'unlabelled' and 'trivial_terms'"
    ),
     make_option(
        c("-g", "--tmpdir"),
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
        c("-m", "--semantic-sim-metric"),
        action = "store",
        default = 'lin',
        type = 'character',
        help = 'Semantic similarity scoring method. 
                Must be supported by Onassis package.
                See listSimilarities()$pairwiseMeasures for a list of accepted options.
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
        c("-s", "--sort-by-agg-score"),
        action = "store_true",
        default = TRUE,
        type = 'logical',
        help = 'Should cells be sorted by their aggregated scores? Default: TRUE'
    ),
    make_option(
        c("-l", "--true-labels"),
        action = "store",
        type = 'character',
        default = NA,
        help = '(OPTIONAL) Path to the true labels tsv file in case tool performance is evaluated. 
                Expected columns: cell_id, true_label, ontology_term'
    ),
    make_option(
        c("-o", "--summary-table-output-path"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the output table with top labels and per-cell metrics in .tsv format'
    ),
    make_option(
        c("-r", "--raw-table-output-path"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the output table with all labels in .tsv format'
    )

)

# parse arguments 
opt = wsc_parse_args(option_list, mandatory = c("input_dir", "cl_dictionary", 
                                                "summary_table_output_path",
                                                "raw_table_output_path"))
# source function definitions 
lab_cl_mapping = readRDS(opt$cl_dictionary)
include_siml = opt$include_sem_siml
sim_metric = opt$semantic_sim_metric
script_dir = dirname(strsplit(commandArgs()[grep('--file=', commandArgs())], '=')[[1]][2])
source(file.path(script_dir, 'cell_types_utils.R'))
# import the rest of dependencies
suppressPackageStartupMessages(require(hash))
suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(foreach))
suppressPackageStartupMessages(require(parallel))
suppressPackageStartupMessages(require(doParallel))
suppressPackageStartupMessages(require(yaml))

# retrieve tool scores if specified
if(!is.na(opt$tool_table)){
    tool_table = read.csv(opt$tool_table, sep="\t", stringsAsFactors=FALSE)
    tools = tool_table[, "Tool"]
    tool_scores = tool_table[, "Combined_score"]
    names(tool_scores) = tools
} else{
    tool_scores = NULL
}
ontology = import_ontology_graph(opt$tmpdir, opt$ontology_graph)

# check if true labels provided
true_labs_provided = !is.na(opt$true_labels)

# read in predicted labels; check barcode consistency 
file_names = list.files(opt$input_dir, full.names=TRUE)
predicted_labs_tables = lapply(file_names, function(f) read.csv(f, sep="\t", 
                                                       stringsAsFactors=FALSE,
                                                       comment.char = "#"))

# extract tools that produced given predictions
source_tools = as.character(sapply(file_names, function(f) extract_metadata(f)[['tool']]))
cell_ids = get_unq_cell_ids(predicted_labs_tables)

# read in exclusions file, if provided
if(! is.na(opt$exclusions)){
    e = yaml.load_file(opt$exclusions)
    unlabelled = tolower(e$unlabelled)
}

# process label columns to know which tool generated which label
# and what datasets a label maps to
lab_dataset_mapping = hash()
comb_ds = list()
for(idx in seq_along(predicted_labs_tables)){
    tbl = predicted_labs_tables[[idx]]
    tool_name = source_tools[idx]
    # select columns with labels
    labels = tbl[ , which(startsWith(colnames(tbl), "label"))]
    # shift to lowercase
    labels = sapply(labels, tolower)
    colnames(labels) = paste(tool_name, c(1:ncol(labels)), sep="_") 
    predicted_labs_tables[[idx]] = labels

    # extract datasets of origin 
    datasets = tbl[, which(startsWith(colnames(tbl), "dataset"))]
    colnames(datasets) = paste(tool_name, "dataset", c(1:ncol(datasets)), sep="_") 
    comb_ds[[idx]] = datasets
    # map labels to corresponding datasets 
    l = as.vector(t(labels))
    d = as.vector(t(datasets))
    lab_dataset_mapping = build_label_dataset_mapping(l, d, lab_dataset_mapping) 
}
# combine predicted labels and corresponding datasets into single data frames
labels = do.call(cbind, predicted_labs_tables)
comb_ds = do.call(cbind, comb_ds)

###################################################################
# Calculate metrics for combined labels and select top candidates
###################################################################
# top labels based on frequency (accounted for tool scores, if provided)
top_labs = apply(labels, 1, function(row) get_top_labels(row, tool_scores=tool_scores))
top_labs = t(top_labs)

# semantic similarity across predicted labels 
.get_sem_sim = function(iter){
    label_vec = as.character(top_labs[iter, c(1:3) ])
    sem_sim = matrix(nrow=length(label_vec), ncol=length(label_vec))
    for(i in 1:length(label_vec)){
        label_i = label_vec[i]
        for(j in i:length(label_vec)){
            label_j = label_vec[j]
            # skip if NA values present
            if(!is.na(label_i) & !is.na(label_j)){
                sem_sim[i,j] = get_CL_similarity(label_i, label_j, 
                                            lab_cl_mapping=lab_cl_mapping,
                                            ontology=ontology,
                                            sim_metric=sim_metric,
                                            unlabelled=unlabelled)
            } else{
                sem_sim[i,j] = NA
            }
        }
    }
    return(round(mean(sem_sim, na.rm=TRUE),3))
}

# semantic similarity to true labels, if provided
.get_sem_sim_true_labs = function(iter){
    pred_lab = as.character(top_labs[[iter]][1]) # select first consensus label
    true_lab = true_labels[iter, "true_label"]
	
    if(!is.na(pred_lab) & !is.na(true_lab)){
        sem_siml = get_CL_similarity(pred_lab, true_lab, 
                                        lab_cl_mapping=lab_cl_mapping,
                                        ontology=ontology,
                                        sim_metric=sim_metric,
                                        unlabelled=unlabelled)
    } else{
        sem_siml = NA
    }
    return(round(sem_siml, 3))
}


# agreement among predicted labels
agreement_rate = apply(labels, 1, get_agreement_rate)

# get proportion of unlabedlled cells (important to distinguish between cells
# with high agreement and poorply labelled cells)
unlab_rate = apply(labels, 1, function(lab_vec) get_unlab_rate(lab_vec, unlabelled))

n_cells = nrow(labels)

# include true labels & update dictionary
if(true_labs_provided){ 
    true_labels = read.csv(opt$true_labels, sep="\t", stringsAsFactors=FALSE)  

    if(! all(c("cell_id", "true_label", "ontology_term") %in% colnames(true_labels))){
        stop("Incorrect field names in true labels file provided.")
    }

    # make sure cell ids correspond
    if(! setequal(true_labels$cell_id, cell_ids)){ 
        stop("Provided true labels do not match those in predicted data")
    }

    # make sure order is the same
    true_labels = true_labels[ match(true_labels$cell_id, cell_ids), ]

    # update the label-CL terms dictionary 
    .upd_dictionary = function(key, val, hash){      
	if(has.key(key, hash)){
            return()
        }
        hash[tolower(key)] = val
    }
    for(idx in  1:n_cells){
	.upd_dictionary(true_labels$true_label[idx],
                        true_labels$ontology_term[idx],
                        lab_cl_mapping)
	}
}

sem_sim = NA
# if specified, calculate semantic similarity across labels
if(include_siml){
    # run computations in parallel, if specified 
    if(opt$parallel){
        if(is.na(opt$num_cores)){
            n_cores = detectCores()
        } else {
            n_cores = opt$num_cores
        }
        registerDoParallel(n_cores)
        # average semantic similarity across labels
        avg_siml = foreach (iter=1:n_cells) %dopar% {
            .get_sem_sim(iter)
        }
        sem_sim = do.call(rbind, avg_siml)
        
        if(true_labs_provided){
	    registerDoParallel(n_cores)
            siml_to_true_labs = foreach(iter=1:n_cells) %dopar% {
                .get_sem_sim_true_labs(iter)
            }
	}

    } else{
        # sequential execution
        avg_siml = lapply(1:n_cells, function(idx) .get_sem_sim(idx))
        sem_sim = do.call(cbind, avg_siml)
        if(true_labs_provided){
            siml_to_true_labs = sapply(1:n_cells, function(idx) .get_sem_sim_true_labs(idx))
        }
    }
}

# extract dataset(s) of origin for top labels
ds_tbl = apply(top_labs, 1, function(row) extract_datasets(row, lab_dataset_mapping))
ds_tbl = data.frame(t(ds_tbl))
colnames(ds_tbl) = paste("dataset", c(1:3), sep="_")

# calculate combined score for each cell
combined_score = sapply(1:n_cells, function(idx) mean(c(agreement_rate[idx], (1 - unlab_rate[idx]), sem_sim[idx]), na.rm=TRUE))
top_labs_tbl = data.frame(cbind(cell_id=cell_ids, 
                                top_labs,
                                agreement_rate=agreement_rate,
                                unlab_rate=unlab_rate,
                                comb_score=combined_score, 
                                ds_tbl))

if(include_siml){
    top_labs_tbl["mean_sem_siml"] = sem_sim
}

if(true_labs_provided){
    top_labs_tbl["siml_to_true_labs"] = unlist(siml_to_true_labs)
}

if(opt$sort_by_agg_score){
    top_labs_tbl = top_labs_tbl[ order(-top_labs_tbl$comb_score), ]
}
write.table(top_labs_tbl, file=opt$summary_table_output_path, sep="\t", row.names=FALSE)

comb_labels = data.frame(cbind(cell_id=cell_ids, labels, comb_ds))
write.table(comb_labels, file=opt$raw_table_output_path, sep="\t", row.names=FALSE)

#TODO: add p-values for semantic similarity? 
#TODO: finding most specific common ancestor of given terms
