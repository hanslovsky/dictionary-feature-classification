function [acc, labels_gt, labels_pd, vals_pd ] = ...
      comp_features(this, cl_type, max_points_to_use, ...
                 edge_ids, pe_batch_size, save_feat_prefix, ...
                 augs, num_workers, set_norm_params)
  % cl_type:
  %   1 = train
  %   3 = test
  
  global DFEVAL_DIR; 
  global DAWMRLIBPATH;
  
  
  %% initialization                               
  if(~exist('edge_ids','var') || isempty(edge_ids))
    edge_ids = 1:this.ds.num_outputs;
  end
  
  skip_training_acc = 0;
  if(nargout == 0)
    skip_training_acc = 1;
  end
  
  if(~exist('save_feat_prefix','var'))
    save_feat_prefix = [];
  end  
  if(~exist('max_points_to_use','var'))
    max_points_to_use = [];
  end
  if(~exist('pe_batch_size','var') || isempty(pe_batch_size))
    pe_batch_size = 1;
  end
  if(~exist('num_workers','var') || isempty(num_workers))
    num_workers = 2000;
  end
  if(~exist('set_norm_params','var') || isempty(set_norm_params))
    set_norm_params = 1;
  end
  if(~exist('augs','var'))
    augs = [];
  end
  assert(cl_type==1 || isempty(augs), ...
         'DAWMRLIB:AssertionFailed', ...
         'augs should be empty when running testing');
  
  % TODO: change files_ filename to something with table name?
  while(true)
    temp_dir = sprintf('%s/files_%s_%s', DFEVAL_DIR, ...
                       datestr(now,30), get_random_id([],1));
    if(exist(temp_dir,'dir')==0)
      break
    end
    pause(10);
  end
  system(sprintf('mkdir %s', temp_dir));
  fprintf('[files dir=%s]\n', temp_dir);
  
  data_size = cell(1, length(this.dds));
  
  max_edge_id = max(edge_ids);
  acc         = zeros(1,max_edge_id);
  labels_gt   = cell (1,max_edge_id);
  labels_pd   = cell (1,max_edge_id);
  vals_pd     = cell (1,max_edge_id);
  
  labels_cl   = cell(1,max_edge_id);
  valid_cl    = cell(1,max_edge_id);

  %%% set features in files/db
  
  if(cl_type == 1 && ~isempty(this.ds_additional))
    ds_all = [this.ds this.ds_additional];
  else
    ds_all = this.ds;
  end
  if(cl_type == 3)
    ds_mean_noise_old  = this.ds.mean_noise;
    ds_std_noise_old   = this.ds.std_noise;
    
    this.ds.mean_noise = 0;
    this.ds.std_noise  = 0;
  end
  ds_orig = this.ds;
  
  % TODO: use/remove?
  max_id_nonaug = []; %#ok<NASGU>
  db_offset = 0;
  for ds_i = ds_all
    this.ds = ds_i;
    
    for i=1:length(this.dds)
      dd = this.dds(i);
      res = dd.scaling;
      if(dd.downsampling==1)
        res = [1 1 1];
      end
      data_size{i} = ...
          this.ds.get_data_size(cl_type, res);
    end
    
    [labels_cl_tmp, valid_cl_tmp, ...
     db_offset_new] = ...
      this.set_features(cl_type, edge_ids, max_points_to_use, ...
                                 data_size, save_feat_prefix, ...
                                 temp_dir, ...
                                 pe_batch_size, ...
                                 db_offset, ...
                                 length(ds_all), ...
                                 augs, num_workers, set_norm_params);
    db_offset = db_offset_new;
            
    for edge_id = edge_ids
      labels_cl{edge_id} = [labels_cl{edge_id}; ...
                            labels_cl_tmp{edge_id}];
      valid_cl{edge_id}  = ...
            logical([valid_cl{edge_id}; ...
                     valid_cl_tmp{edge_id}]);
    end

  end
  this.ds = ds_orig;
  if(cl_type == 3 && set_norm_params)
    this.ds.mean_noise = ds_mean_noise_old;
    this.ds.std_noise  = ds_std_noise_old;
  end
  
  %%% start training/testing
  
  table_name_orig = this.table_prefix;
  if(~isempty(this.saved_table) && cl_type==1)
    this.table_prefix = this.saved_table;
  end
  
  this_fn  = sprintf('%s/this_%s.mat', temp_dir, ...
                     datestr(now,30));
  save(this_fn, 'this', '-v7.3');

  