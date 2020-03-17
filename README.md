# cell-types-analysis [![Anaconda-Server Badge](https://anaconda.org/ebi-gene-expression-group/atlas-fastq-provider/badges/installer/conda.svg)](https://anaconda.org/ebi-gene-expression-group/cell-types-analysis)
A suite of scripts for analysis of scRNA-seq cell type classification tool outputs. These scripts can be used both for evaluating the existing methods by running pipelines on labelled data and for analysing predicted labels for novel data sets.  

## Installation 
The package is installed via EBI Gene Expression Group conda channel:

```conda install -c ebi-gene-expression-group cell-types-analysis``` 

Use `run_post_install_tests.sh` script to make sure installation was successful. 

Run `<script_name>.R --help` to see help for corresponding script.  


## Usage

### Tool evaluation scenario 
For each tool, a list of statistics is generated. The following metrics are used:
* percentage of unlabelled cells
* change in the percentage of unlabelled cells in predicted labels vs reference labels
* proportion of exact matches between predicted and reference labels (exact matching)
* average proportion of shared words between reference and predicted labels (partial matching)
* accuracy & median F1-score
* average semantic similarity between predicted and reference labels
* combined score: combination of the above metrics expressed as a single number. This score is then used as a 'weight' for labels predicted by corresponding tool in production scenario. 

As a result, a table with a set of metrics per each tool is produced. 
Run `get_tool_performance_table.R` to generate this table. 

For a subset of statistics, empirical distribution function can be estimated and then used to determine the p-values of computed scores. To create empirical distribution and assign p-values to the scores, run `get_empirical_dist.R` and `get_tool_pvals.R`, respectively.    

A common output format across tools is assumed - a tab-separated table with 3 columns: cell id (or barcode), predicted labels, and corresponding prediction scores (e.g. p-values or distance metrics). In cases when scores cannot be retrieved, this column can be omitted. Outputs from multiple tools must be stored in a single directory, with file names prefixed by tool name, e.g. `toolX_output.tsv`. A reference table is required, with the following compulsory columns: cell ids, reference labels, and corresponding cell ontology terms. See the example snippet below: 

```
"cell_id" | "predicted_label" | "score"
--- | --- | --- 
"ERR2632411" | "memory B cell" | 0.8
"ERR2632412" | "memory B cell" | 0.8
"ERR2632413" | "memory B cell" | 0.8
```

### Production scenario (novel data)
In production scenario, we are interested in getting as accurate predictions for novel data as possible. To do so, we run a host of pre-trained classifiers against an incoming data set. Then, we aggregate predictions on a single tool basis and select top ones based on corresponding scores (`combine_tool_outputs.R`).

See output file snippet below: 
```
"cell_id" | "label_1" | "label_2" | "dataset_1" | "dataset_2" | "score_1" | "score_2"
--- | --- | --- | --- | --- | --- | ---
"ERR2632411" | "memory B cell" | "memory B cell" | "E-MTAB-6386" |   "OTHER-DATASET" | 0.8 | 0.9
"ERR2632412" | "memory B cell" | "memory B cell" | "E-MTAB-6386" |   "OTHER-DATASET" | 0.8 | 0.9
```

This is followed by a second filtering step where predictions across different tools are analysed for consistency and semantic similarity (`get_consensus_output.R`). Output file snippet: 

Example output tables can be found [here](https://github.com/ebi-gene-expression-group/cell-types-analysis/tree/master/example_output).

## Commands 
**build_cell_ontology_dict.R**: create a mapping between labels and ontology terms
```
build_cell_ontology_dict.R\
          --input-dir <path to the directory with condensed SDRF files>\
          --condensed-sdrf <Boolean: is the provided SDRF file in a condensed form? Default: TRUE>\
          --barcode-col-name <Name of the barcode column in SDRF files (must be identical across all files)>\
          --cell-label-col-name <Name of the cell label column in SDRF files (must be identical across all files)>\
          --cell-ontology-col-name <Name of the cell ontology terms column in SDRF files (must be identical across all files)>\
          --output-dict-path <Output path for serialised object containing the dictionary>
```

**get_tool_performance_table.R**: create a table of metrics for the analysed tools
```
get_tool_performance_table.R\
          --input-dir <path to directory with standardised tab-separated output files from analysed methods>\
          --ref-file <path to tab-delimited file with reference cell type labels>\
          --num-cores <number of cores to run the process on>\
          --ontology-graph <path to the Cell Ontology graph object>\
          --lab-cl-mapping <path to label - CL terms dictionary in .rds format>\
          --barcode-col-ref <name of the cell id field in reference file>\
          --barcode-col-pred <name of the cell id field in predictions file>\
          --label-column-ref <name of label column in reference file>\
          --label-column-pred <name of label column in prediction file>\
          --semantic-sim-metric <semantic similarity scoring method>\    
          --output-path <path to tab-delimited output table>
```
Note that semantic similarity scoring method must be supported by the Onassis package. To see a list of available methods, run `listSimilarities()$pairwiseMeasures` or refer to Onassis [documentation](https://bioconductor.org/packages/release/bioc/html/Onassis.html). 

**get_empirical_dist.R**: generate empirical cumulative distributions for a set of relevant statistics
```
get_empirical_dist.R\
          --input-ref-file <path to tab-delimited file with reference labels>\
          --label-column-ref <name of label column in reference file>\
          --lab-cl-mapping <path to label - CL terms dictionary in .rds format>\
          --num-iterations <number of sampling iterations to construct empirical distribution>\
          --num-cores <number of cores to run simulations>\
          --ontology-graph <path to the Cell Ontology graph object>\
          --semantic-sim-metric <semantic similarity scoring method>\  
          --output-path <path to serialised output object in .rds format>
```

**get_tool_pvals.R**: get empirical p-values from generated distributions
```
get_tool_pvals.R\
          --input-table <path to table of tool statistics produced by get_tool_performance_table.R>\
          --emp-dist-list <path to serialised list of empirical distributions in .rds format>\
          --output-table <path to tab-delimited output table with added p-values> 
```

**combine_tool_outputs.R**: aggregate predictions from multiple classifiers (for a single tool)
```
combine_tool_outputs.R\
          --input-dir <Path to the directory with standardised output .tsv files from multiple classifiers>\
          --top-labels-num <Number of top labels to keep>\
          --scores <Boolean: Are prediction scores available for the given method? Default: FALSE>\
          --output-table <Path to the output table in text format>
```
Note: files in the input directory are assumed to be of the following structire: `A_B_final-labs.tsv`, where A is dataset or origin and B is classifier used to obtain predictions.

**get_consensus_output.R**: Get most likely labels across all tools
```
get_consensus_output.R\
          --input-dir <Path to the directory with standardised .tsv files from multiple methods>
          --tool-table <Path to the tool evaluation table in text format>\
          --num-cores <Number of cores to run the process>\
          --cl-dictionary <Path to the mapping between labels and CL terms in .rds format>\
          --ontology-graph <Path to the ontology graph in .obo or .xml format>\
          --semantic-sim-metric <Semantic similarity scoring method>\
          --summary-table-output-path <Path to the output table with top labels and per-cell metrics in .tsv format>
          --raw-table-output-path <Path to the output table with all labels in .tsv format>
```
Note: the following naming is expected for input tables, e.g. `tool-X_combined.tsv`
