function [vals_pd, labels_gt] = do_test_inference(...
  this, sn3_mn, sn3_std, edge_ids, job_id_start, job_id_end)
  
  vals_index = 0;
  for job_id = job_id_start:job_id_end
    outer_id   = floor(job_id/this.outer_size);
    inner_id   = job_id - this.outer_size*outer_id;
    job_dir    = sprintf('%s/%d/%d', this.save_dir, ...
                         outer_id, inner_id);
    
    jj_fn   = sprintf('%s/info.mat', job_dir);
    if(~exist(jj_fn,'file')), continue, end
    jj_info = load(jj_fn);
    
    n_locs = size(jj_info.feat_loc,1);
    ids_batch = 1000;    
    for ids_index = 1:ids_batch:n_locs    
      ids       = ids_index:min(ids_index+ids_batch-1, n_locs);
      ids_prev  = max(ids-1, 1);
      ids_match = (ids_prev == ids);
      iter_feat_loc      = jj_info.feat_loc(ids,:);
      iter_feat_loc_prev = jj_info.feat_loc(ids_prev,:);
      iter_feat_loc_prev(ids_match,:) = 0;
        
      feats  = get_feats_from_file(this.save_dir, ...
                                   this.feature_dim, ...
                                   iter_feat_loc, ...
                                   iter_feat_loc_prev);
      feats  = bsxfun(@rdivide, ...
                      bsxfun(@minus, feats, sn3_mn), sn3_std);
      
      n_ids  = length(ids);
      n_out  = max(edge_ids); %length(jj_info.labels);
      labels = zeros(n_ids, n_out);
      
      if(~exist('vals_pd','var'))
        vals_pd   = cell(1,n_out);
        labels_gt = cell(1,n_out);
      end
      
      if(isempty(ids)), continue, end
      
      for ii=edge_ids
        labels(:,ii) = jj_info.labels{ii}(ids,:);
      end
      
      vals_index = vals_index + 1;
      
      for edge_id = edge_ids
        valid_ids = find(~isnan(labels(:,edge_id)));
        
        vals_pd{edge_id}{vals_index,1}   = ...
            this.model{edge_id}.mlp_test( ...
              feats(valid_ids, :), length(valid_ids));
        labels_gt{edge_id}{vals_index,1} = ...
            labels(valid_ids, edge_id);
      end
    end
  end
  
  if(exist('vals_pd','var'))
    for edge_id = edge_ids
      if(~isempty(vals_pd{edge_id}))
        vals_pd{edge_id}   = cell2mat(vals_pd{edge_id});
        labels_gt{edge_id} = cell2mat(labels_gt{edge_id});
      end
    end
  else
    vals_pd = {};
    labels_gt = {};
  end
end

% TODO: remove duplicate in get_training_batch
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
