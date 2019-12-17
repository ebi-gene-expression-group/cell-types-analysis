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
output_dir=$test_working_dir/outputs

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
export input_dir='data/test/pred_labs/'
export ref_labels_file='data/test/metadata.tsv'
export ontology_graph='data/cl-basic.obo'
export tool_perf_table=$output_dir/'tool_perf_table.tsv'
export cell_anno_table=$output_dir/'cell_anno_table.tsv'
export empirical_dist=$output_dir/'empirical_dist_list.rds'
export tool_table_pvals=$output_dir/'tool_pvals.tsv'

export num_iter=5
export num_cores=4
export use_existing_outputs

# Derive the tests file name from the script name
tests_file="${script_name%.*}".bats
# Execute the bats tests
$tests_file




