classdef ec_mlp_h5 < ec_mlp_abs
% mlp using file system
  
  properties
    resource_limited = 0
    threaded         = 1
    
    outer_size = 100;
    feats_per_file = 10;
  end
  
  methods
    function obj = ec_mlp_h5(mlp_init)
      if(nargin == 0)
        mlp_init = [];
      end
      obj = obj@ec_mlp_abs(mlp_init);
    end
    
    function initialize(this, exp_name, feature_dim, ...
                        num_points, save_dir)
      initialize@end_classifier_abs(this, exp_name, feature_dim, ...
                                    num_points, save_dir);
    end

    do_features_training(this, job_id, inner_id, feats, valid, ...
                         edge_ids, ...
                         labels, weighting)

    [vals_pd, labels_gt] = ...
        do_training(this, obj_fn, ...
                    labels_cl, valid_cl, edge_ids, ...
                    skip_training_acc, max_job_id, ...
                    dfeval_dir, dawmrlib_dir)
    
    [train_feats, train_labels, train_data, got_all] = ...
        get_training_batch(this, num_subinstances, ...
                                 edge_ids, pos_ratio_copy, ...
                                 train_data)
    
    [vals_pd, labels_gt] = do_test_inference(...
      this, sn3_mn, sn3_std, edge_ids, job_id_start, job_id_end)
  end
end
