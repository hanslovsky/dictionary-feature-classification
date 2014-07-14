function do_features_training(this, job_id, inner_id, ...
                              feats, valid, ...
                              edge_ids, ...
                              labels, weighting) 

  if(isempty(feats)), return, end
  
  % feat is ( n_feats x n_samples )
                          
  job_outer_id = floor(job_id/this.outer_size);
  
  job_dir    = sprintf('%s/%d', this.save_dir, ...
                       job_outer_id);
                   
  n_feats       = size(feats,2);
  
  labPos = ( labels{1} > 0 );
  labNeg = ~labPos;
  
  fn = fullfile(job_dir,sprintf('features-%d.h5', job_id))
  N = [ nnz(labNeg) nnz(labPos) ];
    
  if(~exist(job_dir,'dir'))
    mkdir_st   = system(['mkdir -p ' job_dir]); %#ok<NASGU>
  end
  
  if( ~exist( fn, 'file'))
      if(N(1) > 0 )
        h5create( fn, '/labels/1', [ n_feats, N(1)  ], ...
            'datatype', 'single', ...
            'chunksize', max(ceil([ 0.1 .* n_feats, 0.1 .* N(1) ]), 1), ...
            'deflate', 9, ...
            'fillvalue', single(0));
      end
      if( N(2) > 0 )
        h5create( fn, '/labels/2', [ n_feats, N(2)  ],...
            'datatype', 'single', ...
            'chunksize', max(ceil([ 0.1 .* n_feats, 0.1 .* N(2) ]), 1), ...
            'deflate', 9, ...
            'fillvalue', single(0));
      end
  end
    
%   info = h5info( fn, '/labels/1');
%   sz1 = info.Dataspace.Size;
%   
%   info = h5info( fn, '/labels/2');
%   sz2 = info.Dataspace.Size;
  
  if(N(1) > 0 )
    negdat = feats( labNeg, :)';
    h5write( fn, '/labels/1', negdat );
  end
  
  if( N(2) > 0)
    posdat = feats( labPos, :)';
    h5write( fn, '/labels/2', posdat );
  end
%   h5write( fn, '/labels/1', negdat, [1 (sz1(2)+1)], size(negdat));
%   h5write( fn, '/labels/2', posdat, [1 (sz2(2)+1)], size(posdat) );
  
end