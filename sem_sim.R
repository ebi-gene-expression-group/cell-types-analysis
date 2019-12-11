#!/usr/bin/env Rscript

# TODO: add onassis to env specification 
suppressPackageStartupMessages(require(Onassis))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

# accept as input 2 sequences of CL ID terms for predicted and reference cell types
# return a vector of similarities of the same length 



# test semantic similarity library on data from TM 

# read in metadata file 
metadata = read.table("cell_types_project/data/TM/tm_cell_ontology.tsv", sep = "\t", header = TRUE)
ontology_terms = as.character(metadata$Sample.Characteristic.Ontology.Term.inferred.cell.type.)
cell_types = as.character(metadata$Sample.Characteristic.inferred.cell.type.)

cell_type_id_mapping = hash()
# populate the hash table 
for(idx in 1:length(cell_types)){
    cell_type_id_mapping[[cell_types[idx]]] = ontology_terms[idx]
}
saveRDS(cell_type_id_mapping, "cell_types_project/data/TM/tm_cell_ontology_mapping.rds")


test_data = readRDS("~/dev/reports/scmap_performance_analysis/data/scmap_analysis_data.rds")
predicitons = test_data[[2]][,2]
pred_cl = sapply(predicitons, function(x) cell_type_id_mapping[[x]])
reference = as.character(test_data[[3]])
ref_cl = sapply(reference, function(x) cell_type_id_mapping[[x]])
cl_table = data.frame(cbind(ref_cl, pred_cl))
saveRDS(cl_table, file = "~/cell_types_project/data/TM/cl_similarity_table.rds")


# semantic similarity 
table = readRDS("cl_similarity_table.rds")
find_siml = function(idx){
    tryCatch({
        Similarity(obo, table[idx,1], table[idx,2])},
    error = function(cond){
        return(NA)})
}
