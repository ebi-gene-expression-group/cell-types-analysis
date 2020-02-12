# cell-types-analysis
A suite of scripts for analysis of cell type classification tools outputs. 

## Usage 
For each tool, a list of statistics is generated. The following metrics are used:
* percentage of unlabelled cells
* change in the percentage of unlabelled cells in predicted labels vs reference labels
* proportion of exact matches between predicted and reference labels (exact matching)
* average proportion of shared words between reference and predicted labels (partial matching)
* accuracy & median F1-score
* average semantic similarity between predicted and reference labels

For a subset of statistics, empirical distribution function can be estimated and then used to determine the p-value of obtained scores.  

A common output format across tools is assumed - a tab-separated table with 2 columns: `cell_id` (or barcode) and `predicted_label`. Outputs from multiple tools must be stored in a single directory, with file names prefixed by tool name, e.g. `toolX_output.tsv`. A reference table is required, with the following compulsory columns: cell ids, reference labels, and corresponding cell ontology terms.

Predicted labels for specific cells can be aggregated and analysed, providing evidence in favour of each label based on the overall performance of corresponding tools.

Example output tables can be found [here](https://github.com/ebi-gene-expression-group/cell-types-analysis/tree/master/example_output).

### Installation 
The package is installed via EBI Gene Expression Group conda channel:

```conda install -c ebi-gene-expression-group cell-types-analysis``` 

Use `run_post_install_tests.sh` script to make sure installation was successful. 

Run `<script_name>.R --help` to see documentation.  

### Commands 
**get_tool_performance_table.R**: create a table of metrics for the analysed tools
```
get_tool_performance_table.R\
          --input-dir <path to directory with standardised tab-separated output files from analysed methods>\
          --ref-file <path to tab-delimited file with reference cell type labels>\
          --label-column-ref <name of label column in reference file>\
          --label-column-pred <name of label column in prediction file>\
          --cell-ontology-col <name of cell ontology id column in reference dataset>\
          --semantic-sim-metric <semantic similarity scoring method>\    
          --output-path <path to tab-delimited output table>
```
Note that semantic similarity scoring method must be supported by the Onassis package. To see a list of available methods, run `listSimilarities()$pairwiseMeasures` or refer to Onassis [documentation](https://bioconductor.org/packages/release/bioc/html/Onassis.html). 

**get_empirical_dist.R**: generate empirical cumulative distributions for a set of relevant statistics
```
get_empirical_dist.R\
          --input-ref-file <path to tab-delimited file with reference labels>\
          --label-column-ref <name of label column in reference file>\
          --cell-ontology-col <name of cell ontology id column in reference dataset>\
          --num-iterations <number of sampling iterations to construct empirical distribution>\
          --num-cores <number of cores to run simulations>\
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
**get_cell_annotations_table.R**: calculate cell-level statistics
```
get_cell_annotations_table.R\
          --input-dir <path to directory with standardised tab-separated output files from analysed methods>\
          --ref-file <path to tab-delimited file with reference cell type labels>\
          --tool-table <path to table of tool statistics produced by get_tool_performance_table.R>\
          --cell-ontology-col <name of cell ontology id column in reference dataset>\
          --semantic-sim-metric <semantic similarity scoring method>\  
          --label-column-ref <name of label column in reference file>\
          --label-column-pred <name of label column in prediction file>\
          --output-path <path to tab-delimited output table>
```
