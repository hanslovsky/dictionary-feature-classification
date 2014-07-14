function [train_feats, train_labels, train_data, got_all] = ...
    get_training_batch(this, num_subinstances, ...
                             edge_ids, pos_ratio_copy, ...
                             train_data)

  got_all = 0;
  
  if(~exist('train_data','var') || isempty(train_data))
    train_data_fn = [this.save_dir '/train_data.mat'];
    if(exist(train_data_fn, 'file'))
      train_data = load(train_data_fn);
      train_data = train_data.train_data;
      n_out = size(train_data.labels,2);
    else
      train_data.feat_loc  = zeros(0,5);
      train_data.weighting = zeros(0,1);
      
      dir_outer     = dir(this.save_dir);
      for do_i = dir_outer(3:end)'
        if(~do_i.isdir || isnan(str2double(do_i.name)))
          continue
        end
        
        outer_id  = str2double(do_i.name);
        dir_inner = dir(sprintf('%s/%d', this.save_dir, ...
                                outer_id));
        for di_i = dir_inner(3:end)'
          if(~di_i.isdir || isnan(str2double(di_i.name)))
            continue
          end

          inner_id = str2double(di_i.name);
          % jj       = outer_id * this.outer_size + inner_id;
      
          job_dir    = sprintf('%s/%d/%d', this.save_dir, ...
                               outer_id, inner_id);

          jj_fn = sprintf('%s/info.mat', job_dir);
          if(~exist(jj_fn, 'file')), continue, end
          jj_info = load(jj_fn);
        
          if(isempty(jj_info.labels) || ...
             isempty(jj_info.labels{1}))
            continue
          end
        
          if(~isfield(train_data, 'labels'))
            if(isempty(edge_ids))
              n_out = length(jj_info.labels);
              edge_ids = 1:n_out;
            else
              n_out = length(edge_ids);
            end
            train_data.labels = zeros(0,n_out,'single');
          end
          
          train_data.feat_loc(end+(1:size(jj_info.feat_loc,1)),:) ...
            = jj_info.feat_loc;
          curr_index = size(train_data.labels,1);
          for ii=1:length(edge_ids)
            ee = edge_ids(ii);
            train_data.labels( ...
              curr_index+(1:size(jj_info.labels{ee},1)),ii) = ...
              jj_info.labels{ee};
          end
          train_data.weighting(end+(1:size(jj_info.weighting,1))) ...
            = jj_info.weighting;
        end
      end
      
      save(train_data_fn, 'train_data', '-v7.3');
    end
  else
    n_out = size(train_data.labels,2);
  end
  
  num_points          = size(train_data.feat_loc,1);
  sample_fractions    = -1*ones(1,n_out);
  sample_fractions(:) = pos_ratio_copy;
  n_sampled           = sum(sample_fractions > 0);
  if(n_sampled == 0)
    if(num_subinstances < num_points)
      ids = randsample(num_points, num_subinstances, true, ...
                       train_data.weighting);
    else
      ids = (1:num_points)';
    end
    
    ids_prev  = max(ids-1, 1);
    ids_match = (ids_prev == ids);
    iter_feat_loc      = train_data.feat_loc(ids,:);
    iter_feat_loc_prev = train_data.feat_loc(ids_prev,:);
    iter_feat_loc_prev(ids_match,:) = 0;
    train_feats  = get_feats_from_file(this.save_dir, ...
                                       this.feature_dim, ...
                                       iter_feat_loc, ...
                                       iter_feat_loc_prev);
    train_labels = train_data.labels(ids,:);
  else
    if(isempty(edge_ids))
      edge_ids = 1:n_out;
    end
    batch_size = ceil(num_subinstances/n_sampled);
    
    label_vals = [ -1 1 ];
    ee_index   = 0;
    for ii = 1:length(edge_ids)
      sf_val = sample_fractions(ii);
      if(sf_val <= 0), continue, end
      
      pos_ratio_vals = [ (1-sf_val) (sf_val) ];
      for vv = 1:2
        valid_ids = ...
            find(train_data.labels(:,ii) == label_vals(vv));
        ee_batch_size = ...
            ceil(batch_size*pos_ratio_vals(vv));
        ids = randsample(...
          valid_ids, ee_batch_size, true, ...
          train_data.weighting(valid_ids));
        
        ids_prev  = max(ids-1, 1);
        ids_match = (ids_prev == ids);
        iter_feat_loc      = train_data.feat_loc(ids,:);
        iter_feat_loc_prev = train_data.feat_loc(ids_prev,:);
        iter_feat_loc_prev(ids_match,:) = 0;
        train_feats{ee_index+vv}  = get_feats_from_file(...
          this.save_dir, this.feature_dim, ...
          iter_feat_loc, iter_feat_loc_prev); %#ok<AGROW>
        train_labels{ee_index+vv} = train_data.labels(ids,:); %#ok<AGROW>
      end
      ee_index = ee_index + 2;
    end
  end
  
end

function feats = get_feats_from_file(base_dir, feat_dim, ...
                                               feat_loc, ...
                                               feat_loc_prev)
  n_feats = size(feat_loc,1);
  feats = zeros([n_feats feat_dim], 'single');

  parfor i=1:n_feats
    feats(i,:) = CompressLib.decompress(...
      get_feat(base_dir, feat_loc(i,:), ...
                         feat_loc_prev(i,:)) );
  end
end

function ff = get_feat(base_dir, ifeat_loc, ifeat_loc_prev)
  [fid, msg] = fopen(sprintf('%s/%d/%d/%d/data_%d', ...
                             base_dir, ...
                             ifeat_loc(1), ...
                             ifeat_loc(2), ...
                             ifeat_loc(3), ...
                             ifeat_loc(4))); %#ok<NASGU>
  if(...%id == 1 || ...
     ifeat_loc_prev(1) ~= ifeat_loc(1) || ...
     ifeat_loc_prev(2) ~= ifeat_loc(2) || ...
     ifeat_loc_prev(3) ~= ifeat_loc(3) || ...
     ifeat_loc_prev(4) ~= ifeat_loc(4))
    prev_loc = 0;
  else
    prev_loc = ifeat_loc_prev(5);
  end
  
  fseek(fid, prev_loc, 'bof');
  ff = fread(fid, ifeat_loc(5)-prev_loc, '*uint8');
  fclose(fid);
end
