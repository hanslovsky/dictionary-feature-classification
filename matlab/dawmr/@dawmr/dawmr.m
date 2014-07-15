classdef dawmr < handle
% dawmr class:
% Kmeans-based Affinity graph generator
  
  properties (Access = public)
    verbose       = 0;
    rand_sampling = 0;
    rand_ordering = [];
    balance_train = 0;
    balance_test  = 0;
    
    % 0 = no normalization
    % 1 = scale to [0,1], assume min is already 0
    % 2 = scale to [-1,1]
    % 3 = zero mean, unit std
    svm_normalization = 0;
    sn1_max = [];
    sn2_min = [];
    sn2_max = [];
    sn3_mn  = [];
    sn3_std = [];
    
    default_num_patches_kmeans = [];
    
    end_classifier
    
    saved_table = [];
    
    ds
    ds_additional
    dds = dawmr_data.empty;
    
  end
    
  properties %(GetAccess = public, SetAccess = private)    
    table_prefix
  end
  
  methods % defined externally
    learn_features(this, num_patches, learning_type, locs, ...
                   offset, imsize, cube_size)
    feats = get_features(this, inds, x, y, z, cl_type, ...
                               pe_batch_size, ...
                               aug, locs)
    [acc, labels_gt, labels_pd, vals_pd, auc, ...
     model_out, sn_out] = ...
        classifier(this, cl_type, max_points_to_use, ...
                   edge_ids, pe_batch_size, save_feat_prefix, ...
                   augs, num_workers)
               
    [acc, labels_gt, labels_pd, vals_pd ] = ...
        comp_features(this, cl_type, max_points_to_use, ...
                 edge_ids, pe_batch_size, save_feat_prefix, ...
                 augs, num_workers, set_norm_params)
    
    infer(this, offset, imdim, output_fn, ...
          val_shift, pe_batch_size, cube_size, outer_cube_size, ...
          new_data_fn, flag, use_sparse, n_recur)
    
    [tm, feats] = get_features_timing(this, xs, ys, zs, flag, ...
                                            aug)
    
    inds = set_inds(this, cl_type, mask, labels, ...
                          edge_ids, max_points_to_use, ...
                          output_counts)
    
    [labels_cl, valid_cl, ...
     db_offset_new] = ...
        set_features(this, cl_type, ...
                           edge_ids, max_points_to_use, ...
                           data_size, save_feat_prefix, ...
                           temp_dir, pe_batch_size, ...
                           db_offset, ...
                           num_outer_augmentations, ...
                           augs, num_workers, set_norm_params)
                       
    [labels_cl, valid_cl, m1, m2, m3] = ...
        set_features_dist( ...
          this, xs, ys, zs, cl_type, edge_ids, ...
          max_points_to_use, data_size, ...
          job_id, augs, output_counts)
      
     [labels_cl, valid_cl, feats ] = ... 
        dawmr_feat_comp_dist(this, cl_type, edge_ids, max_points_to_use, ... 
                       data_size, save_feat_prefix, temp_dir, ... 
                       pe_batch_size, db_offset, ... 
                       num_outer_augmentations, augs, num_workers)
                             
    [feats] = ...
        get_features_dist( ...
          this, xs, ys, zs, cl_type, edge_ids, ...
          max_points_to_use, data_size, ...
          job_id, augs, output_counts)
    
    min_offset = get_min_offset(this)
    max_offset = get_max_offset(this)
  end
  
  methods
    function obj = dawmr(ds, svm_normalization, end_classifier, ...
                      table_prefix)
      if(exist('ds','var') && isa(ds, 'dawmr_set'))
        obj.ds = ds;
      end
      
      if(exist('svm_normalization', 'var') && ...
         ~isempty(svm_normalization))
        obj.svm_normalization = svm_normalization;
      end
      
      assert(isa(end_classifier, 'end_classifier_abs'), ...
             'DAWMRLIB:AssertionFailed', ...
             'invalid end_classifier');
      obj.end_classifier = end_classifier;
      
      if(exist('table_prefix', 'var') && ...
         ~isempty(table_prefix))
        set_table_prefix(obj, table_prefix);
      end
    end
    
    function t = get_table_name(this)
      t = sprintf('%s', this.table_prefix);
    end
    function set_table_prefix(this, t)
      usr = getenv('USER');
      this.table_prefix = sprintf('%s_%s_%s', ...
                                  usr, t, ...
                                  datestr(now,30));
      
      if(length(this.table_prefix) > 64)
        table_short = [this.table_prefix(1:49) ...
                       datestr(now,30)];
        fprintf('warning: table length too long\n');
        fprintf('  shortening to: %s\n', table_short);
        this.table_prefix = table_short;
      end
      
      %, ...
      %get_random_id([],1));
    end
    
    function add_dd(this, dd)
      this.dds(end+1) = dd;
      dd.dawmr_p         = this;
    end
    
    function copy_dds(this, k)
      for dd = k.dds
        this.add_dd(dd);
      end
    end
    
    function n = num_clustering_cache(this)
      n = 0;
      for dd=this.dds
        n = n + length(dd.dcs);
      end
    end
    
    function n = num_pooling_cache(this)
      n = 0;
      for dd=this.dds
        for dc=dd.dcs
          n = n + length(dc.dps);
        end
      end
    end
    
    function n = feature_vec_length(this)
      n = 0;
      for dd=this.dds
        n = n + dd.feature_vec_length();
      end
    end
    
    function l = learned(this)
      l = false;
      for dd = this.dds
        if(~dd.learned())
          return
        end
      end
      l = true;
    end
    
    % function move_models(this, new_dir)
    %   assert(this.end_classifier == 0, 'DAWMRLIB:AssertionFailed', ...
    %          'only can be used with svm saved model files');
    %   assert(exist(new_dir, 'dir')>0, ...
    %          'DAWMRLIB:AssertionFailed', ...
    %          sprintf('directory %s does not exist', new_dir));
      
    %   for i=1:length(this.model)
    %     if(~isempty(this.model{i}))
    %       assert(exist(this.model{i},'file')>0, ...
    %              'DAWMRLIB:AssertionFailed', ...
    %              'model file does not exist');
          
    %       [~,model_fn,model_ext] = fileparts(this.model{i});
          
    %       new_model_fn = sprintf('%s/%s%s', new_dir, ...
    %                              model_fn, model_ext);
          
    %       system(sprintf('mv %s %s', this.model{i}, new_model_fn));
          
    %       this.model{i} = new_model_fn;
    %     end
    %   end
    % end
    
  end
end
