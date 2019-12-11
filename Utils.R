### find number of unassigned cells ###
get_unlab_rate = function(labels, unlabelled){
    labels = as.character(labels)
    prop_unknown_provided = length(labels[(labels %in% unlabelled)])/length(labels)
    prop_unknown_provided = round(prop_unknown_provided, 3)*100
    return(prop_unknown_provided)
}

#### find ratio of unlabelled cells in reference dataset to that in predicted dataset ###
get_unlab_rate_ratio = function(prop_unlab_reference, prop_unlab_predicted){
    if(prop_unlab_reference != 0){
        cat("Prop unlab != 0") 
        unlab_ratio = round((prop_unlab_predicted/prop_unlab_reference)*100,3)
    } else{
        cat("Prop unlab = 0")
        unlab_ratio = round(prop_unlab_predicted*100 / (prop_unlab_reference+1), 3)
    }
    return(unlab_ratio)
}

#### find proportion of exact matches between reference cell types and predicted labels ###
get_exact_matches = function(reference_labels, predicted_labels){
    matches = length(reference_labels[reference_labels == predicted_labels])
    prop_exact_match = round(matches/length(reference_labels), 3)*100
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
    intersect_ratio = round(intersect_ratio, 3)*100
    return(intersect_ratio)
}

### find per-class/median F-1 scores and accuracy ###
# adapted from https://github.com/tabdelaal/scRNAseq_Benchmark/blob/master/evaluate.R
get_preformance_metrics = function(reference_labs, predicted_labs, unlabelled) {
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
    findLabel = colnames(conf_F1) == row.names(conf_F1)[i]
    if(sum(findLabel)){
      prec = conf_F1[i, findLabel] / colSums(conf_F1)[findLabel]
      rec = conf_F1[i, findLabel] / rowSums(conf_F1)[i]
      if (prec == 0 || rec == 0){
        F1[i] = 0
      } else{
        F1[i] = (2*prec*rec) / (prec + rec)
      }
      sum_acc = sum_acc + conf_F1[i, findLabel]
    } else {
      F1[i] = 0
    }
  }
  # filter out empty classes
  pop_size = pop_size[pop_size > 0]
  names(F1) = names(pop_size)
  med_F1 = median(F1)
  # calculate overall accuracy 
  acc = sum_acc/sum(conf_F1)
  out = list(MedF1 = med_F1, Acc = acc)
  return(out)
}

get_CL_similarity = function(reference_labs_df, ref_labs_col, CL_col, predicted_labs) {
    suppressPackageStartupMessages(require(Onassis))
    suppressPackageStartupMessages(require(org.Hs.eg.db)) # TODO: need argument to control database 
    suppressPackageStartupMessages(require(hash))

    ontology_terms = as.character(reference_labs_df[, CL_col])
    cell_types = as.character(reference_labs_df[, ref_labs_col])
    
    # map cell types to CL terms 
    cell_type_id_mapping = hash()
    .set(cell_type_id_mapping, keys=cell_types, values=ontology_terms)
    # get CL terms for predicted labels
    obo = "data/cl-basic.obo"
    pred_cl = sapply(predicted_labs, function(x) cell_type_id_mapping[[x]])
    ref_cl = as.character(reference_labs_df[, CL_col])

    # semantic similarity 
    .find_siml = function(idx){
        tryCatch({
            siml = Similarity(obo, pred_cl[idx], ref_cl[idx])
            return(siml)},
        error = function(cond){
            return(NA)})
    }

    similarity = sapply(1:length(ref_cl), .find_siml)
    mean_siml = round(mean(similarity, na.rm = TRUE), 3)
    return(mean_siml)
}












