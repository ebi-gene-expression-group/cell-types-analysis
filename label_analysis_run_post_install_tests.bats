#!/usr/bin/env bats

@test "Build label - CL term mapping" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$label_cl_dict" ]; then
        skip "$label_cl_dict exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $label_cl_dict && build_cell_ontology_dict.R\
					--input-dir $SDRF_dir\
					--condensed-sdrf\
					--output-dict-path $label_cl_dict\
                    			--cell-label-col-name $SDRF_cell_types\
					--output-text-path $label_cl_txt

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$label_cl_dict" ]
}

@test "build tool evaluation table" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$tool_perf_table" ]; then
        skip "$tool_perf_table exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $tool_perf_table && get_tool_performance_table.R\
					--input-dir $eval_input_dir\
					--barcode-col-ref $barcode_col_ref\
					--parallel $parallel\
					--num-cores $num_cores\
					--barcode-col-pred $barcode_col_pred\
                    			--tmpdir $tmpdir\
                    			--ontology-graph $ontology_graph\
					--label-column-ref $label_column_ref\
					--label-column-pred $label_column_pred\
					--lab-cl-mapping $label_cl_dict\
					--ref-file $ref_labels_file\
                    			--include-sem-siml \
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
					--parallel $parallel\
					--num-iterations $num_iter\
					--sample-labs $sample_labs\
                    			--tmpdir $tmpdir\
					--lab-cl-mapping $label_cl_dict\
					--num-cores $num_cores\
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

@test "combine results" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$combined_results" ]; then
        skip "$combined_results exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $combined_results && combine_tool_outputs.R\
					--input-dir $res_to_combine\
					--top-labels-num $top_labels_num\
					--scores\
					--output-table $combined_results

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$combined_results" ]

}

@test "Get consensus output" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$summary_table_path" ]; then
        skip "$summary_table_path exists and use_existing_outputs is set to 'true'"
    fi

    run rm -rf $summary_table_path && get_consensus_output.R\
                     --input-dir $combined_tools_results\
                     --tool-table $tool_perf_table\
                     --cl-dictionary $lab_cl_mapping_consensus\
                     --tmpdir $tmpdir\
                     --parallel $parallel\
                     --num-cores $num_cores\
                     --include-sem-siml\
                     --summary-table-output-path $summary_table_path\
                     --raw-table-output-path $raw_labels_table_path\
		             --true-labels $true_labels


    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f "$summary_table_path" ]
}

@test "Run matrix down-sampling" {
    if [ "$use_existing_outputs" = 'true' ] && [ -d "$sampling_out_dir" ]; then
        skip "$sampling_out_dir exists and use_existing_outputs is set to 'true'"
    fi

    run rm -rf $sampling_out_dir && downsample_cells.R\
                    --expression-data $sampling_test_10x_data\
                    --metadata $sampling_test_sdrf\
                    --cell-id-field "$sampling_cell_id_field"\
                    --cell-count-threshold $cell_num_threshold\
                    --cell-type-field "$sampling_cell_type_field"\
                    --array-size-limit $sampling_arr_size_limit\
                    --output-dir $sampling_out_dir\
                    --metadata-upd $sampling_metadata_upd

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -d "$sampling_out_dir" ]
}

@test "Check labels" {
    if [ "$use_existing_outputs" = 'true' ] && [ -d "$SDRF_checked_output" ]; then
        skip "$SDRF_checked_output exists and use_existing_outputs is set to 'true'"
    fi

    run rm -rf $SDRF_checked_output && check_labels.R --input-file $SDRF_to_process\
                                                      --condensed\
                                                      --label-field $SDRF_cell_types\
                                                      --output-path $SDRF_checked_output


    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f "$SDRF_checked_output" ]
}





