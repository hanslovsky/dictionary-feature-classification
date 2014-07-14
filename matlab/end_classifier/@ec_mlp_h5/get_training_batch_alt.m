function [train_feats, train_labels, train_data] = ...
    get_training_batch_alt(this, max_job_id, ...
                                 num_subinstances, ...
                                 edge_id, pos_ratio_copy, ...
                                 train_data)
  % alternative method of pulling features from files
  %   aggregating into one single pull
  %   does not seem to be any faster than original version
  % does not support arbitrary output size
  
  if(~exist('train_data','var') || isempty(train_data))
    train_data_fn = [this.save_dir '/train_data.mat'];
    if(exist(train_data_fn, 'file'))
      train_data = load(train_data_fn);
      train_data = train_data.train_data;
    else
      train_data.feat_loc  = zeros(0,4);
      train_data.labels    = zeros(0,3,'single');
      train_data.weighting = zeros(0,1);
      
      for jj=1:max_job_id
        outer_id   = floor(jj/100);
        inner_id   = jj - this.outer_size*outer_id;
        job_dir    = sprintf('%s/%d/%d', this.save_dir, ...
                             outer_id, inner_id);

        jj_info = load(sprintf('%s/info.mat', job_dir));
        train_data.feat_loc(end+(1:size(jj_info.feat_loc,1)),:) = ...
            jj_info.feat_loc;
        train_data.labels(end+(1:size(jj_info.labels{1},1)),:) = ...
            [jj_info.labels{1}, ...
             jj_info.labels{2}, ...
             jj_info.labels{3}];
        train_data.weighting(end+(1:size(jj_info.weighting,1))) = ...
            jj_info.weighting;
      end
      
      save(train_data_fn, 'train_data', '-v7.3');
    end
  end
  
  num_points = size(train_data.feat_loc,1);
  if(pos_ratio_copy < 0)
    if(num_subinstances < num_points)
      ids = randsample(num_points, num_subinstances, true, ...
                       train_data.weighting);
    else
      ids = (1:num_points)';
    end
    
    train_feats  = get_feats_from_file(this.save_dir, ...
                                       this.feature_dim, ...
                                       ids, train_data.feat_loc);
    train_labels = train_data.labels(ids,:);
  else
    if(isempty(edge_id))
      edge_id = 1:3;
    end
    batch_size = ceil(num_subinstances/3);
    
    ee_index = 0;
    label_vals     = [ -1 1 ];
    pos_ratio_vals = [ (1-pos_ratio_copy) (pos_ratio_copy) ];
    
    ids_all = [];
    
    for ee = edge_id
      for vv = 1:2
        valid_ids = ...
            find(train_data.labels(:,ee) == label_vals(vv));
        ee_batch_size = ...
            ceil(batch_size*pos_ratio_vals(vv));
        ids = randsample(...
          valid_ids, ee_batch_size, true, ...
          train_data.weighting(valid_ids));
        
        train_feats_count(ee_index+vv) = length(ids);
        ids_all = [ids_all; ids];
        train_labels{ee_index+vv} = train_data.labels(ids,:); %#ok<AGROW>
      end
      ee_index = ee_index + 2;
    end
    
    train_feats_cumsum = cumsum([0 train_feats_count]);
    train_feats_all = get_feats_from_file(...
      this.save_dir, this.feature_dim, ...
      ids_all, train_data.feat_loc);
    for ee_index = 1:length(train_labels)
      train_feats{ee_index} = ...
          train_feats_all(...
            train_feats_cumsum(ee_index) + ...
            (1:train_feats_count(ee_index)), :); %#ok<AGROW>
    end
  end
  
end

function feats = get_feats_from_file(base_dir, feat_dim, ...
                                               ids, feat_loc)
  n_feats = length(ids);
  feats = zeros([n_feats feat_dim], 'single');

  parfor i=1:n_feats
    feats(i,:) = CompressLib.decompress(...
      get_feat(base_dir, ids(i), feat_loc) );
  end
end

function ff = get_feat(base_dir, id, feat_loc)
  [fid, msg] = fopen(sprintf('%s/%d/%d/data_%d', ...
                             base_dir, ...
                             feat_loc(id,1), ...
                             feat_loc(id,2), ...
                             feat_loc(id,3))); %#ok<NASGU>
  if(id == 1 || ...
     feat_loc(id-1,2) ~= feat_loc(id,2) || ...
     feat_loc(id-1,3) ~= feat_loc(id,3))
    prev_loc = 0;
  else
    prev_loc = feat_loc(id-1,4);
  end
  
  fseek(fid, prev_loc, 'bof');
  ff = fread(fid, feat_loc(id,4)-prev_loc, '*uint8');
  fclose(fid);
end
