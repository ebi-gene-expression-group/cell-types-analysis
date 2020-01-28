#!/usr/bin/env bats 

@test "build tool evaluation table" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$tool_perf_table" ]; then
        skip "$tool_perf_table exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $tool_perf_table && get_tool_performance_table.R\
                                    --input-dir $input_dir\
                                    --barcode-col-ref $barcode_col_ref\
                                    --barcode-col-pred $barcode_col_pred\
                                    --label-column-ref $label_column_ref\
                                    --label-column-pred $label_column_pred\
                                    --cell-ontology-col $cell_ontology_col\
                                    --ref-file $ref_labels_file\
                                    --ontology-graph $ontology_graph\
                                    --output-path $tool_perf_table
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$tool_perf_table" ]
  
}

@test "generate empirical CDF" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$empirical_dist" ]; then
        skip "$empirical_dist exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $empirical_dist && get_empirical_dist.R\
                                    --input-ref-file $ref_labels_file\
                                    --label-column-ref $label_column_ref\
                                    --cell-ontology-col $cell_ontology_col\
                                    --num-iterations $num_iter\
                                    --num-cores $num_cores\
                                    --ontology-graph $ontology_graph\
                                    --output-path $empirical_dist
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$empirical_dist" ]
  
}

@test "obtain p-values for calculated statistics" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$tool_table_pvals" ]; then
        skip "$tool_table_pvals exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $tool_table_pvals && get_tool_pvals.R\
                                    --input-table $tool_perf_table\
                                    --emp-dist-list $empirical_dist\
                                    --output-table $tool_table_pvals

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$tool_table_pvals" ]
}

@test "get per-cell statistics" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$cell_anno_table" ]; then
        skip "$cell_anno_table exists and use_existing_outputs is set to 'true'"
    fi
    
    run rm -f $cell_anno_table && get_cell_annotations_table.R\
                                    --input-dir $input_dir\
                                    --ref-file $ref_labels_file\
                                    --barcode-col-ref $barcode_col_ref\
                                    --barcode-col-pred $barcode_col_pred\
                                    --cell-ontology-col $cell_ontology_col\
                                    --label-column-ref $label_column_ref\
                                    --label-column-pred $label_column_pred\
                                    --ontology-graph $ontology_graph\
                                    --tool-table $tool_perf_table\
                                    --output-path $cell_anno_table
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$cell_anno_table" ]
   
}
