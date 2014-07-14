function [vals_pd, labels_gt] = ...
      do_training(this, obj_fn, ...
                  labels_cl, valid_cl, edge_ids, ...
                  skip_training_acc, max_job_id, ...
                  dfeval_dir, dawmrlib_dir) 

  mlp_pe_batch = 1;
  use_gpu_copy = this.mlp_init.use_gpu;
  
  obj      = load(obj_fn);
  sn3_mn    = obj.this.sn3_mn{1}; 
  sn3_std   = obj.this.sn3_std{1}; 
  
  this.save_dir
  
  save_fn = sprintf('%s/saved_obj.mat', this.save_dir);
  save(save_fn, 'this', 'sn3_mn', 'sn3_std');

  dawmr_lib_exe = get_dawmr_lib_exe();
  
  q_ret = qsub_dist(@ec_run_mlp_joint, ...
                    mlp_pe_batch, ...
                    use_gpu_copy, 0, ...
                    [], dawmr_lib_exe, ...
                    {save_fn}, ...
                    {labels_cl}, {valid_cl}, ...
                    {edge_ids}, {max_job_id}, ...
                    {dfeval_dir}, {dawmrlib_dir});
  this.model = q_ret{1}{1};

  n_out = max(edge_ids); %length(this.model);
  
  vals_pd   = cell(1,n_out);
  labels_gt = cell(1,n_out);
  if(~skip_training_acc)    
    ids_per_worker = ceil(max_job_id/1000);
    ids_start = 1:ids_per_worker:max_job_id;
    ids_end   = min(ids_start + (ids_per_worker-1), max_job_id);
    ids_start = num2cell(ids_start);
    ids_end   = num2cell(ids_end);

    dawmr_lib_exe = get_dawmr_lib_exe(1);
    
    save(save_fn, 'this');
    obj_method_name = 'do_test_inference';
    q_ret = qsub_dist(...
      @run_obj_method_dist, ...
      [], [], [], [], dawmr_lib_exe, ...
      {save_fn}, {'this'}, {obj_method_name}, {2}, ...
      {sn3_mn}, {sn3_std}, {edge_ids}, ids_start, ids_end);
    
    for edge_id=edge_ids
      num_instances      = sum(valid_cl{edge_id});
      vals_pd{edge_id}   = zeros(num_instances, 1);
      labels_gt{edge_id} = zeros(num_instances, 1);
    end
    
    indices = zeros(1,n_out);
    
    for jj = 1:size(q_ret,1)
      job_vals_pd   = q_ret{jj}{1}{1};
      job_labels_gt = q_ret{jj}{1}{2};

      if(isempty(job_vals_pd)), continue, end
      
      for edge_id = edge_ids
        n_edge_labels = size(job_vals_pd{edge_id},1);
        
        vals_pd{edge_id}(  indices(edge_id)+(1:n_edge_labels)) = ...
            job_vals_pd{edge_id};
        labels_gt{edge_id}(indices(edge_id)+(1:n_edge_labels)) = ...
            job_labels_gt{edge_id};
        
        indices(edge_id) = indices(edge_id) + n_edge_labels;
      end
    end
  end
  
end
