### find number of unassigned cells ###
get_unlab_rate = function(labels, unlabelled){
    labels = as.character(labels)
    prop_unknown_provided = length(labels[(labels %in% unlabelled)])/length(labels)
    prop_unknown_provided = round(prop_unknown_provided, 3)
    return(prop_unknown_provided)
}

#### find ratio of unlabelled cells in reference dataset to that in predicted dataset ###
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
get_f1 = function(reference_labs, predicted_labs, unlabelled) {
    unique_ref = unique(reference_labs)
    unique_pred = unique(predicted_labs)
    #unique_all = unique(c(unique(unique_ref, unique_pred)))
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

get_CL_similarity = function(reference_labs, ref_CL_terms, predicted_labs, ontology, sim_metric, unlabelled) {
    suppressPackageStartupMessages(require(Onassis)) #NB: keep package call within function, otherwise parallelism is broken
    # initialise and configure Similarity object 
    siml_object = new('Similarity')
    ontology(siml_object) = ontology
    # configure similarity measurement metric
    if(! sim_metric %in% listSimilarities()$pairwiseMeasures){
        stop("Incorrect semantic similarity metric provided.")
    }
    idx = which(listSimilarities()$pairwiseMeasures == sim_metric)
    pairwiseConfig(siml_object) = listSimilarities()$pairwiseMeasures[idx]
    #pairwiseConfig(siml_object) = listSimilarities()$pairwiseMeasures[5]
    # map cell types to CL terms 
    cell_type_id_mapping = hash()
    i = which(!(reference_labs %in% unlabelled | ref_CL_terms==''))
    .set(cell_type_id_mapping, keys=reference_labs[i], values=ref_CL_terms[i])
    # get CL terms for predicted labels
    pred_cl = sapply(predicted_labs, function(x) cell_type_id_mapping[[x]])

    # semantic similarity 
    .find_siml = function(idx){
        tryCatch({
            siml = pairsim(siml_object, as.character(pred_cl[idx]), ref_CL_terms[idx])
            return(siml)},
        error = function(cond){
            print(cond)
            print("returning NA")
            return(NA)})
    }

    similarity = sapply(1:length(ref_CL_terms), .find_siml)
    mean_siml = round(mean(similarity, na.rm = TRUE), 3)
    return(mean_siml)
}

# get a combined score across all metrics 
get_tool_combined_score = function(unlab_delta,
                                   exact_match_prop,
                                   mean_shared_terms, 
                                   median_F1, 
                                   accuracy, 
                                   siml) { 
    
    score = (1 - unlab_delta) + exact_match_prop + mean_shared_terms + median_F1 + siml
    return(round(score, 3))

}

#########################################
# combine function calls  
#########################################
obtain_metrics_list = function(tool,
                               reference_labs,
                               predicted_labs,
                               prop_unlab_reference,
                               unlabelled,
                               trivial_terms,
                               ref_CL_terms,
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
    siml = get_CL_similarity(reference_labs, ref_CL_terms, predicted_labs, ontology, sim_metric, unlabelled)
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
get_agreement_rate = function(label_list){
    # assume 1st column is cell id, 2nd col is label, and exclude them
    label_list = label_list[-c(1,2)]
    agreement = 1 / length(unique(label_list))
    return(round(agreement, 3))
}

# mean pairwise similarity among predicted labels 
get_cell_CL_siml = function(siml_object, label_list, labs_dict) {
    label_list = label_list[-c(1,2)]
    siml = c()
    for(i in 1:length(label_list)){
        term_i = labs_dict[[label_list[i]]]
        for(j in 1:length(label_list)){
            term_j = labs_dict[[label_list[j]]]
            tryCatch({
                siml = c(siml, pairsim(siml_object, as.character(term_i), as.character(term_j)))
            },
            error = function(cond){
                print(cond)
                siml = c(siml, NA)
            })
        }
    }
    mean_siml = mean(siml, na.rm = TRUE)
    return(round(mean_siml, 4))
}

# calculate aggregated score for predicted labels 
get_weighted_score = function(labels_list, tool_scores){
    if(!all(names(labels_list) == names(tool_scores))){
        stop("Name mismatch between predicted labels and tools")
    }
    .get_sums = function(lab){
        t = names(labels_list[labels_list == lab])
        s = sum(tool_scores[t])
        return(s)
    }

    sums = sapply(labels_list, function(label) .get_sums(label)) 
    labels_list = sapply(1:length(labels_list), function(idx) paste(labels_list[idx], sums[idx], sep=": "))
    return(labels_list)
}

