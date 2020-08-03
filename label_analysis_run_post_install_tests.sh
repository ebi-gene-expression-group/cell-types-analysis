#!/usr/bin/env bash
script_name=$0

function usage {
    echo "usage: run_post_install_tests.sh [action] [use_existing_outputs]"
    echo "  - action: what action to take, 'test' or 'clean'"
    echo "  - use_existing_outputs, 'true' or 'false'"
    exit 1
}

action=${1:-'test'}
use_existing_outputs=${2:-'false'}

if [ "$action" != 'test' ] && [ "$action" != 'clean' ]; then
    echo "Invalid action"
    usage
fi

if [ "$use_existing_outputs" != 'true' ] && [ "$use_existing_outputs" != 'false' ]; then
    echo "Invalid value ($use_existing_outputs) for 'use_existing_outputs'"
    usage
fi

test_working_dir=`pwd`/'post_install_tests'
output_dir=$test_working_dir/'outputs'
data_dir=$test_working_dir/'data'

# Clean up if specified
if [ "$action" = 'clean' ]; then
    echo "Cleaning up $test_working_dir ..."
    rm -rf $test_working_dir
    exit 0
elif [ "$action" != 'test' ]; then
    echo "Invalid action '$action' supplied"
    exit 1
fi 

# Initialise directories
mkdir -p $test_working_dir
mkdir -p $output_dir

################################################################################
# List tool outputs/inputs & parameters 
################################################################################
export eval_input_dir=$data_dir/'results_dir/'
export res_to_combine=$data_dir/'prod_outputs_scpred/'
export ref_labels_file=$data_dir/'reference_sdrf.tsv'

export sampling_test_10x_data=$data_dir/'10x_data/'
export sampling_test_sdrf=$data_dir/'E-GEOD-83139.sdrf.txt'
export sampling_cell_id_field='FactorValue..single.cell.identifier.'
export sampling_cell_type_field='Characteristics..inferred.cell.type.'
export sampling_arr_size_limit=2454500
export sampling_out_dir='10x_data_sampled'
export sampling_metadata_upd='E-GEOD-83139.sdrf_sampled.txt'

export SDRF_dir=$data_dir/'SDRFs'
export SDRF_cell_types="cell type"
export tmpdir="TMPDIR"
export ontology_graph="https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl-basic.obo"
export combined_tools_results=$data_dir/'prod_combined_tools'
export barcode_col_pred='cell_id'
export barcode_col_ref='Assay'
export label_column_ref='Sample.Characteristic.cell.type.'
export label_column_pred='predicted_label'

export label_cl_dict=$output_dir/'label_cl_dict.rds'
export label_cl_txt=$output_dir/'label_cl_txt.tsv'
export tool_perf_table=$output_dir/'tool_perf_table.tsv'
export empirical_dist=$output_dir/'empirical_dist_list.rds'
export tool_table_pvals=$output_dir/'tool_pvals.tsv'
export combined_results=$output_dir/'combined_results.tsv'

export raw_labels_table_path=$output_dir/'raw_labels_table.tsv'
export summary_table_path=$output_dir/'summary_output_table.tsv'

export num_iter=5
export sample_labs=50
export parallel='TRUE'
export num_cores=4
export top_labels_num=2
export use_existing_outputs

# retrieve test data 
test_data_url="ftp://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/atlas/cell-types-project-test-data/cell_types_analysis_test_data.tar.gz"
test_data_archive=$test_working_dir/`basename $test_data_url`

if [ ! -e "$test_data_archive" ]; then
    wget $test_data_url -P $test_working_dir
    tar -xzvf $test_data_archive -C $test_working_dir
fi

# Derive the tests file name from the script name
tests_file="${script_name%.*}".bats
# Execute the bats tests
$tests_file
