function success = writeFeaturesH5( feats, labels, fn )
% success = writeFeaturesH5( feats, labels, fn )
%
% feats  - the features (N x P vector (single))
%          or a cell array, with one entry per label
% labels - the labels( N x 1 vector)
%

success = false;

% list of labels in the dataset
if( iscell(feats))
    labelList = 1:length(feats);
    
else
    if( ~isa( feats, 'single') )
        fprintf('converting to float32\n');
        feats = single(feats);
    end
    sz = size(feats);
    labelList = unique(labels);
end


try
  if( ~exist(fn,'file'))
    
    % write the features for this label
    for l = 1:length(labelList)
        
        if( iscell(feats))
            writeme = feats{l};
            if( ~isa( writeme, 'single') )
                fprintf('converting to float32\n');
                writeme = single(writeme);
            end
        else
            idxs = labels == labelList(l);
            N = nnz( idxs );
            szLabel = [ N sz(2) ];
            writeme = feats( idxs, :);
        end

        datasetString = sprintf('/labels/%d', labelList(l));
        h5create(fn, datasetString, szLabel, ...
                    'datatype', 'single', ...
                    'fillvalue', single(0));
        
        h5write( fn, datasetString, writeme' );
    end
    
  else
     error('file exists ... exiting');
  end
  
catch e
  error( e.message );
  return;
end % trycatch


success = true;
