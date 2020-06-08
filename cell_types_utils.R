#!/usr/bin/env Rscript

# Define global variables 
# list of unlabelled cell types
# NB: keep these relevant to the tools' output
unlabelled = c("unassigned", "unknown", 'rand','node','ambiguous', NA, "not available")

# add common words here to improve specificity of word-matching method 
trivial_terms = c("cell", "of", "tissue", "entity", "type") 

### find number of unassigned cells ###
get_unlab_rate = function(labels, unlabelled){
    labels = as.character(labels)
    prop_unknown_provided = length(labels[(labels %in% unlabelled)])/length(labels)
    prop_unknown_provided = round(prop_unknown_provided, 3)
    return(prop_unknown_provided)
}

#### find delta between unlabelled cells in reference dataset and those in predicted dataset ###
get_unlab_rate_delta = function(prop_unlab_reference, prop_unlab_predicted){
    unlab_delta = abs(prop_unlab_predicted - prop_unlab_reference)
    return(round(unlab_delta, 3))
}

#### find proportion of exact matches between reference cell types and predicted labels ###
get_exact_matches = function(reference_labels, predicted_labels){
    matches = length(reference_labels[reference_labels == predicted_labels])
    prop_exact_match = round(matches/length(reference_labels), 3)
    return(prop_exact_match)
}

### find mean proportion of shared terms between reference and predicted labels ###
get_shared_terms_prop = function(reference_labs, pred_labs, trivial_terms){
    # split predicted and provided annotations into words 
    split_ref = sapply(reference_labs, function(x) strsplit(as.character(x), " "))
    split_prediction = sapply(pred_labs, function(x) strsplit(as.character(x), " "))
    # find corresponding terms between word vectors 
    intersection = sapply(1:length(split_ref), 
                        function(x) intersect(split_ref[[x]], split_prediction[[x]]))
    # remove trivial terms from intersection 
    intersection = sapply(intersection, function(x) x[!x %in% trivial_terms])
    # find average proportion of matching words 
    intersect_ratio = mean(sapply(1:length(intersection), 
                        function(x) length(intersection[[x]])/length(split_prediction[[x]])))
    intersect_ratio = round(intersect_ratio, 3)
    return(intersect_ratio)
}

### find per-class/median F-1 scores and accuracy ###
# adapted from https://github.com/tabdelaal/scRNAseq_Benchmark/blob/master/evaluate.R
get_f1 = function(reference_labs, predicted_labs, unlabelled=unlabelled) {
    unique_ref = unique(reference_labs)
    unique_pred = unique(predicted_labs)
    conf = table(reference_labs, predicted_labs)
    pop_size = rowSums(conf)
    
    # confusion table for existing cell types, excluding unassigned 
    conf_F1 = table(reference_labs,predicted_labs, exclude = unlabelled)
    F1 = c()
    sum_acc = 0

    # iterate across classes and calculate per-class F1-score
    for (i in c(1:length(unique_ref))){
    # check if label is in confusion table
    find_label = colnames(conf_F1) == row.names(conf_F1)[i]
    if(sum(find_label, na.rm = TRUE)){
      prec = conf_F1[i, find_label] / colSums(conf_F1)[find_label]
      rec = conf_F1[i, find_label] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] = (2*prec*rec) / (prec + rec)
      }
      sum_acc = sum_acc + conf_F1[i, find_label]
    } else {
      F1[i] = 0
    }
  }
  # filter out empty classes
  pop_size = pop_size[pop_size > 0]
  names(F1) = names(pop_size)
  med_F1 = round(median(F1), 3)
  # calculate overall accuracy 
  acc = round(sum_acc/sum(conf_F1), 3)
  out = list(MedF1 = med_F1, Acc = acc)
  return(out)
}

# calculate similarity between a pair of labels 
get_CL_similarity = function(label_1, label_2, lab_cl_mapping, ontology, sim_metric, unlabelled) {
    suppressPackageStartupMessages(require(Onassis)) #NB: keep package import within function, otherwise parallelism is broken
    # initialise and configure Similarity object 
    siml_object = new('Similarity')
    ontology(siml_object) = ontology
    # configure similarity measurement metric
    if(! sim_metric %in% listSimilarities()$pairwiseMeasures){
        stop("Incorrect semantic similarity metric provided.")
    }
    idx = which(listSimilarities()$pairwiseMeasures == sim_metric)
    pairwiseConfig(siml_object) = listSimilarities()$pairwiseMeasures[idx]
    # map labels to CL terms; return 0 if unlabelled 
    if(label_1 %in% unlabelled | label_2 %in% unlabelled){
        return(0)
    }
    # check if labels in dict keys; return NA otherwise 
    if(! label_1 %in% keys(lab_cl_mapping) && label_2 %in% keys(lab_cl_mapping)){
        return(NA)
    }
    term_1 = lab_cl_mapping[[as.character(label_1)]]
    term_2 = lab_cl_mapping[[as.character(label_2)]]

    # semantic similarity 
    tryCatch({
        psim = pairsim(siml_object, term_1, term_2)
        return(psim)
        },
    error = function(cond){
        print(cond)
        print("returning 0")
        return(0)})
}

# get a combined score across all metrics 
get_tool_combined_score = function(unlab_delta,
                                   exact_match_prop,
                                   mean_shared_terms, 
                                   median_F1, 
                                   accuracy, 
                                   siml) { 
    
    score = exact_match_prop + mean_shared_terms + median_F1 + siml
    return(round(score, 3))

}

#########################################
# combine function calls  
#########################################
obtain_metrics_list = function(tool,
                               reference_labs,
                               predicted_labs,
                               prop_unlab_reference,
                               lab_cl_mapping,
                               unlabelled,
                               trivial_terms,
                               ontology,
                               sim_metric){

    # proportion of unlabelled cells in predicted labels
    prop_unlab_predicted = get_unlab_rate(predicted_labs, unlabelled)
    # ratio of unlabelled cells proportions
    unlab_delta = get_unlab_rate_delta(prop_unlab_reference, prop_unlab_predicted)
    # propotion of exact matches 
    exact_match_prop = get_exact_matches(reference_labs, predicted_labs)
    # match based on shared terms
    mean_shared_terms = get_shared_terms_prop(reference_labs, predicted_labs, trivial_terms)
    # calculate F1 scores and accuracy
    metrics = get_f1(reference_labs, predicted_labs, unlabelled)
    median_F1 = metrics$MedF1
    accuracy = metrics$Acc
    # cell ontology similarity
    sim_vec = c()
    for(idx in seq_along(reference_labs)){
        lab_1 = reference_labs[idx]
        lab_2 = predicted_labs[idx]
        sim_vec[idx] = get_CL_similarity(lab_1, lab_2, 
                                         lab_cl_mapping=lab_cl_mapping, 
                                         ontology=ontology, 
                                         sim_metric=sim_metric, 
                                         unlabelled=unlabelled)
    }
    siml = round(mean(sim_vec, na.rm=TRUE), 3)
    score = get_tool_combined_score(unlab_delta, exact_match_prop, mean_shared_terms, 
                                   median_F1, accuracy, siml)

    # combine metrics
    row = list(Tool = tool,
                Unlab_ref = prop_unlab_reference,
                Unlab_pred = prop_unlab_predicted,
                Unlab_delta = unlab_delta,
                Exact_match_prop = exact_match_prop,
                Mean_partial_match = mean_shared_terms,
                Med_F1 = median_F1,
                Accuracy = accuracy,
                CL_similarity = siml, 
                Combined_score = score)
    return(row)
}

################################################################################
# methods for per-cell statistics
################################################################################

# extract list of unique cell ids
get_unq_cell_ids = function(tables_list){
    cell_ids = lapply(tables_list, function(tbl) tbl[ , "cell_id"])
    if(!length(unique(cell_ids)) == 1){
        stop("Inconsistent cell IDs provided")
    }
    cell_ids = unlist(unique(cell_ids))
    return(cell_ids)
}

# get agreement across predicted labels 
get_agreement_rate = function(label_list){
    label_list = label_list[which(!is.na(label_list))]
    n_unq = length(unique(label_list))
    if(n_unq == 0){
        return(NA)
    }
    agreement = 1 / n_unq
    return(round(agreement, 2))
}

# find the most frequent labels across top candidates from each tool; if specified, get weighted scores 
get_top_labels = function(label_list, tool_scores=NULL){
    unq_labs = unique(label_list)
    # determine frequency of each unique label 
    .get_freq = function(lab, n_labs){
        freq = as.numeric(table(label_list)[lab] / n_labs)
        return(round(freq, 3))
    }
    freqs = sapply(unq_labs, function(lab) .get_freq(lab, n_labs=length(label_list)))

    # if tool scores provided, find weighted scores
    if(!is.null(tool_scores)){
        # extract tool name for each label and find corresponding score
        # NB: label vector must contrain corresponding tool in its name, e.g. tool-X_1
        source_tools = sapply(names(label_list), function(name) unlist(strsplit(name, "_"))[1])
        label_scores = sapply(source_tools, function(tool) tool_scores[tool])
        # find sum of scores corresponding to each label 
        score_sums = sapply(unq_labs, function(lab) sum(label_scores[which(label_list==lab)], na.rm=TRUE))
        # account for label frequencies and sort
        sorted = sort(score_sums * freqs, decreasing=TRUE, index.return=TRUE, na.last=TRUE)
        score_names = paste("weighted-score", c(1:3), sep="_")
        } else{
            # otherwise, simply return top labels by frequency
            sorted = sort(freqs, decreasing=TRUE, index.return=TRUE, na.last=TRUE)
            score_names = paste("frequency", c(1:3), sep="_")
        }

        sorted_scores = sorted$x
        sorted_scores = sapply(sorted_scores, function(s) round(s,2))
        sorted_labs = unq_labs[sorted$ix]
        # select top 3 labels; add NAs if there are fewer than 3
        while(length(sorted_scores) < 3){
            sorted_labs = append(sorted_labs, NA)
            sorted_scores = append(sorted_scores, NA)
        }
        sorted_labs = sorted_labs[1:3]
        names(sorted_labs) = paste("label", c(1:3), sep="_")
        sorted_scores = sorted_scores[1:3]
        names(sorted_scores) = score_names
        return(append(sorted_labs, sorted_scores))
}


# create a mapping between predicted labels and datasets of origin 
# it might be that one label comes from more than one dataset
# in this case, combine datasets names into vector.
# takes as input vectors of labels and corresponding datasets + dictionary to fill in
# returns updated dictionary   
build_label_dataset_mapping = function(labels, datasets, lab_ds_map){
    for(idx in seq_along(labels)){
        lab = labels[idx]
        ds = datasets[idx]
        if(is.na(lab)){
            next
        }
        if(lab %in% keys(lab_ds_map)){
            if(!ds %in% lab_ds_map[[lab]]){
                lab_ds_map[lab] = append(lab_ds_map[[lab]], ds)
            } else {
                next
            }
        } else {
            lab_ds_map[lab] = ds
        }
    }
    return(lab_ds_map)
}

# extract datasets of origin by building a string of accession codes
# take as input a vector of labels 
# return a vector of corresponding dataset(s)
extract_datasets = function(label_vec, lab_ds_map){
    label_vec = label_vec[1:3]
    s = c()
    for(idx in seq_along(label_vec)){
        l = label_vec[idx]
        if(l %in% keys(lab_ds_map)){
            ds_vec = lab_ds_map[[l]]
            } else{
            ds_vec = NA
            }
        if(length(ds_vec) > 1){
            s[idx] = paste(ds_vec, collapse="; ")
        } else{
            s[idx] = ds_vec
        }
    }
    return(s)
}

# extract metadata from file into a list 
extract_metadata = function(file){
    con = file(file, "r")
    vals = list()
    idx = 1
    while(TRUE){
        line = readLines(con, 1)
        if(!startsWith(line, "#")){
            return(vals)
        }
        data = unlist(strsplit(line, " "))
        # 1st item is '#', 2nd is field name 
        n = data[2]
        # extract values 
        vals[[idx]] = data[-c(1,2)]
        names(vals)[idx] = tolower(n)
        idx = idx + 1 
    }
}
