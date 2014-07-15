function [labels_cl, valid_cl, ...
          db_offset_new] = ...
    set_features(this, cl_type, edge_ids, max_points_to_use, ...
                       data_size, save_feat_prefix, temp_dir, ...
                       pe_batch_size, db_offset, ...
                       num_outer_augmentations, augs, num_workers, set_norm_params)
  
  if(~exist('db_offset','var') || isempty(db_offset))
    db_offset = 0;
  end
  if(~exist('set_norm_params','var') || isempty(set_norm_params))
    set_norm_params = 1;
  end
  
  if(isempty(augs))
    n_augs = 1;
  else
    n_augs = length(augs);
  end
  
  
  
  n_out = max(edge_ids);%this.ds.num_outputs;
  bal_t = this.end_classifier.balanced_training;
  assert(isscalar(bal_t) || length(bal_t) == length(edge_ids), ...
         'DAWMRLIB:AssertionFailed', ...
         sprintf(...
           'balanced sampling invalid (%d outputs, %d fractions)', ...
           length(edge_ids), length(bal_t)));
  
  base_data_size = ...
      this.ds.get_data_size(cl_type,1);
  
  % TODO: clean up following code
  skip_feature_comp = 0;
  if(~isempty(save_feat_prefix) && ...
     ~exist(save_feat_prefix, 'dir'))
    system(sprintf('mkdir %s', save_feat_prefix));
  end
  if(~isempty(save_feat_prefix) && ...
     exist(sprintf('%s/000001.dat', save_feat_prefix), 'file'))
    skip_feature_comp = 1;
  end

  % TODO: fix this to work with saved, always need to initialize
  % TODO: fix estimate of num_points
  est_num_points = max_points_to_use;
  if(isempty(max_points_to_use) || isinf(max_points_to_use) ...
    || max_points_to_use <= 1)
    est_num_points = prod(base_data_size(1:3))/2;
    if(max_points_to_use <= 1)
      est_num_points = est_num_points * max_points_to_use;
    end
  end
  if( isempty(this.saved_table) && ...
      db_offset == 0 )
    this.end_classifier.initialize(...
      this.get_table_name(), ...
      feature_vec_length(this), ...
      est_num_points*num_outer_augmentations*n_augs, ...
      temp_dir);
  end

  % cube_size = 20 also in set_features_dist!
  inner_cube_size = 20;
  outer_cube_size = partition_space( ...
    base_data_size, ...
    (this.get_max_offset - this.get_min_offset) + 1, ...
    inner_cube_size, num_workers);

  num_inner_cubes = ceil(outer_cube_size / inner_cube_size);
  num_inner_cubes = prod(num_inner_cubes);
  index = 0;
  x_b = {}; y_b = {}; z_b = {};
  for ox=1:outer_cube_size(1):base_data_size(1)
    for oy=1:outer_cube_size(2):base_data_size(2)
      for oz=1:outer_cube_size(3):base_data_size(3)
        index       = index + 1;
        x_b{index}  = [ox, min(ox+outer_cube_size(1)-1, ...
                               base_data_size(1))]; %#ok<AGROW>
        y_b{index}  = [oy, min(oy+outer_cube_size(2)-1, ...
                               base_data_size(2))]; %#ok<AGROW>
        z_b{index}  = [oz, min(oz+outer_cube_size(3)-1, ...
                               base_data_size(3))]; %#ok<AGROW>
        
        offsets{index} = [x_b{index}(1), ...
                          y_b{index}(1), ...
                          z_b{index}(1)]; %#ok<AGROW>
        imsizes{index} = [x_b{index}(2)-x_b{index}(1)+1, ...
                          y_b{index}(2)-y_b{index}(1)+1, ...
                          z_b{index}(2)-z_b{index}(1)+1]; %#ok<AGROW>
      end
    end
  end

  
  % if(max_points_to_use > 1)
  %   max_points_to_use = ceil(max_points_to_use / index);
  % end
  dawmr_lib_exe = get_dawmr_lib_exe(1);
  if(max_points_to_use ~= Inf)
    the_labels_fn = this.ds.labels_fn;
    if(iscell(the_labels_fn))
      the_labels_fn = the_labels_fn{cl_type};
    end
    
    qret = qsub_dist(...
      @get_training_stats, ...
      pe_batch_size,[],[],[], dawmr_lib_exe, ...
      {the_labels_fn}, ...
      this.ds.mask_fn(  cl_type), ...
      offsets, imsizes);
    output_counts = sum(cell2mat([qret{:}]'),1);
    output_counts = reshape(output_counts, n_out, 2);
  else
    output_counts = zeros(n_out,2);
  end
  

  %% distributed feature computation
  dc_t = tic;
  this_fn  = sprintf('%s/this_%s.mat', temp_dir, ...
                     datestr(now,30));
  save(this_fn, 'this', '-v7.3');
  
  num_inner_cubes = num_inner_cubes*n_augs; %#ok<NASGU>
  % job_ids = num2cell(db_offset + ( ...
  %   0:num_inner_cubes:(index-1)*num_inner_cubes));
  job_ids = num2cell(db_offset + (1:index));
    
  if( ~skip_feature_comp && ...
      (isempty(this.saved_table) || cl_type == 3) )
        
    % TODO: add back in resource limited flag?
    qret = qsub_dist(...
      @run_obj_method_dist, ...
      pe_batch_size, [], [], [], dawmr_lib_exe, ...
      {this_fn}, {'this'}, {'set_features_dist'}, {5}, ...
      x_b, y_b, z_b, {cl_type}, {edge_ids}, {max_points_to_use}, ...
      {data_size}, job_ids, {augs}, {output_counts});

    n_qret     = size(qret,1);
    labels_cl  = cell(1, n_out);
    valid_cl   = cell(1, n_out);
    cl_indices = zeros(1, n_out);
    for edge_id = edge_ids
      num_entries = sum(cellfun(@(x) length(x{1}{1}{edge_id}), ...
                                qret));
      labels_cl{edge_id} = zeros(num_entries, 1);
      valid_cl{edge_id}  = zeros(num_entries, 1);
    end

    for qq = 1:n_qret
      for edge_id = edge_ids
        num_entries = length(qret{qq}{1}{1}{edge_id});
        
        labels_cl{edge_id}(cl_indices(edge_id) + ...
                           (1:num_entries)) = ...
          qret{qq}{1}{1}{edge_id};
        valid_cl{edge_id}( cl_indices(edge_id) + ...
                           (1:num_entries)) = ...
          qret{qq}{1}{2}{edge_id};
        
        cl_indices(edge_id) = cl_indices(edge_id) + num_entries;
      end
    end
    
    if(cl_type == 1 && this.svm_normalization > 0)
      qq_index = 0;
      for qq=1:n_qret
        if(isempty(qret{qq}{1}{3})), continue, end
        
        if(~exist('local_m','var'))
          n_feats = length(qret{qq}{1}{3});
          local_m = zeros(n_qret, n_feats);
          local_v = zeros(n_qret, n_feats);
          local_n = zeros(n_qret, 1);
        end
        
        qq_index = qq_index + 1;
        local_m(qq_index,:) = qret{qq}{1}{3};
        local_v(qq_index,:) = qret{qq}{1}{4};
        local_n(qq_index,:) = qret{qq}{1}{5};
      end
      
      local_m = local_m(1:qq_index,:);
      local_v = local_v(1:qq_index,:);
      local_n = local_n(1:qq_index,:);

      switch(this.svm_normalization)
        case 1
          assert(false, 'DAWMRLIB:AssertionFailed', ...
                 'not implemented');
          % this.sn1_max{1} = max(cell2mat(m(:,1)), [], 1);
          % this.sn1_max{1} ...
          %   (this.sn1_max{1} < 1e-5) = 1;
        case 2
          assert(false, 'DAWMRLIB:AssertionFailed', ...
                 'not implemented');
          % this.sn2_min{1} = min(cell2mat(m(:,1)), [], 1);
          % this.sn2_max{1} = ...
          %     max(cell2mat(m(:,2)), [], 1) - ...
          %     this.sn2_min{1};
        case 3
            if(set_norm_params)
                [this.sn3_mn{1}, this.sn3_std{1}] = ...
                    global_mean_var(local_m, local_v, local_n);
                std_eps = 1e-3;
                this.sn3_std{1} = sqrt(this.sn3_std{1}) + std_eps;
            end
      end
      
      if(this.end_classifier.do_feature_renormalization)
        % TODO: call 
        %   end_classifier.do_features_training_renormalization
        
        % save(this_fn, 'this', '-v7.3');
        % qsub_dist(@dawmr_load_write_features, ...
        %           pe_batch_size, ...
        %           0, 0, ...
        %           DFEVAL_DIR, dawmrlib_dir, ...
        %           objs_fn, labels_b, valid_b, feat_output_fns, ...
        %           edge_ids_b, db_offsets_b);      
      end
      
    end
  end
  
  system(sprintf('rm %s', this_fn));
    
  dc_t = toc(dc_t);
  fprintf('distributed job time: %g\n', dc_t);
  % feats{cli} = cell2mat(feats{cli}');
    
  db_offset_new = db_offset + index;%num_inner_cubes*index;
end


%% old feat_output_fns stuff
% feat_output_fns_prefix = cell(1,max(edge_id));
% feat_output_fns = cell(1,dc_size);

% if(~isempty(save_feat_prefix))
%   for edge_id = edge_ids
%     feat_output_fns_prefix{edge_id} = ...
%         sprintf('%s/%d', save_feat_prefix, edge_id);
%     system(sprintf('mkdir %s', feat_output_fns_prefix{edge_id}));
%   end
% else
%   if(db_offset == 0)
%     for edge_id = edge_ids
%       system(sprintf('mkdir %s/%d', temp_dir, edge_id));
%       feat_output_fns_prefix{edge_id} = ...
%           sprintf('%s/%d/%s_',temp_dir, edge_id, tm_stmp);
%     end
%   end
% end

% for i=1:dc_size
%   for edge_id = edge_ids
%     feat_output_fns{i}{edge_id} = ...
%         sprintf('%s%06d.mat', ...
%                 feat_output_fns_prefix{edge_id}, i);
%   end
% end

% % reverse indexing of feat_output_fns
% feat_output_fns_tmp = feat_output_fns;
% feat_output_fns = cell(1,max(edge_ids));
% for edge_id = edge_ids
%   for i = 1:dc_size
%     feat_output_fns{edge_id}{i} = ...
%         feat_output_fns_tmp{i}{edge_id};
%   end
% end
