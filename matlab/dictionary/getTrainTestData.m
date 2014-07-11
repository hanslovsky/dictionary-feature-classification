function feats = getTrainTestData( datimfn, labimfn, mskimfn, N, D, patchSize, param )
% Usage:
%   data_set = getTrainTestData( datim, labim, N, param )
%
% N = vector containing number of training examples per class

%%

% datimfn = ds.data_fn{1};
% labimfn = ds.labels_fn{1};
% mskimfn = ds.mask_fn{1};

%%
fprintf('loading images...\n');
datim = read_image_stack( datimfn );

if( ~isempty(labimfn))
    labim = read_image_stack( labimfn );
    
    if( length( size(labim)) == 4)
        labim_bnd = uint8(round(mean(labim,4))) + 1;
    else
        labim_bnd = labim;
    end
    
    labelList = unique(labim_bnd);
    clear labim;
end

if( ~isempty(mskimfn))
    mskim = read_image_stack( mskimfn );
    if( length( size(mskim)) == 4)
        mskim = mskim(:,:,:,1) & mskim(:,:,:,2) & mskim(:,:,:,3);
    end
   mskim_bnd = false( size(datim));
   bnd_rad = (patchSize-1)./2 + 1;
   mskim_bnd(   bnd_rad(1) : end-bnd_rad(1), ...
            bnd_rad(2) : end-bnd_rad(2), ...
            bnd_rad(3) : end-bnd_rad(3)) = true;
        
	mskim = mskim & mskim_bnd;
        
else
   mskim = false( size(datim));
   bnd_rad = (patchSize-1)./2 + 1;
   mskim(   bnd_rad(1) : end-bnd_rad(1), ...
            bnd_rad(2) : end-bnd_rad(2), ...
            bnd_rad(3) : end-bnd_rad(3)) = true;

%    mskim(1:2:end, 1:2:end, 1:2:end) = false;
end

if( ~isempty(N))
    numClasses = length(N);
elseif( exist('labelList','var') ) % if label image is given, but N unspecified
    numClasses = length(labelList);
end

if(~exist('numClasses','var'))
    numClasses = 1;
    labelList = 1;
end
feats = cell( numClasses, 1);

for c = 1:numClasses
    
    label = labelList(c);
    fprintf('generating data for class (%d)\n', label);
    
    % ignore label
    if( ~isempty(labimfn))
        inds = find( mskim & (labim_bnd == label));
    else
        inds = find( mskim);
    end
    numClass = length( inds );
    
    % ignore N, use all point locations
    if( ~isempty(N))
        r = randperm( numClass, N(c));
        inds = inds(r);
    end
    
    [x,y,z] = ind2sub( size(datim), inds);
%     size(x)
%     pause;
    fprintf('  grabbing patches\n');
    patches = grabPatchesSimple( double(datim), patchSize, [], {x,y,z} );
    fprintf('  encoding\n');
    feats{c} = encoding(  patches', D, 'st', 0, param );
    fprintf('  done!\n');
end
