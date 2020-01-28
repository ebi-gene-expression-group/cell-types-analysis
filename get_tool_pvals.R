#!/usr/bin/env Rscript 

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))

### Using empirical CDFs, obtain p-values for statistics calculated in tool performance table

option_list = list(
    make_option(
    c("-i", "--input-table"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the table of tool statistics from get_tool_performance_table.R'
    ),
    make_option(
    c("-d", "--emp-dist-list"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the list of empirical distributions in .rds format'
    ),
    make_option(
    c("-o", "--output-table"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the modified output table in text format'
    )
)

opt = wsc_parse_args(option_list, mandatory = c("input_table", "emp_dist_list", "output_table"))

distr_list = readRDS(opt$emp_dist_list)
tools_table = read.delim(opt$input_table, stringsAsFactors=FALSE)
metrics = names(distr_list)
for(metric in metrics){
    metric_distr = distr_list[[metric]]
    metric_values = tools_table[, metric]
    tool_pvals = sapply(metric_values, function(x) round((1-metric_distr(x)), 3))
    col_name = paste(metric, "_pval", sep = '')
    tools_table[, col_name] = tool_pvals
}
write.table(tools_table, file = opt$output_table, sep="\t", row.names=FALSE)
