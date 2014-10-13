function [feats, patches, xyz ] = getTrainTestData( datimfn, labimfn, mskimfn, N, D, patchSize, param, dsFactorZ )
% Usage:
%  [feats, patches, xyz ] = getTrainTestData( datimfn, labimfn, mskimfn, N, D, patchSize, param, dsFactorZ )
%
% N - vector containing number of training examples per class
% D - the dictionary

%%

% datimfn = ds.data_fn{1};
% labimfn = ds.labels_fn{1};
% mskimfn = ds.mask_fn{1};

%%
isD = 0;
if( exist( 'D', 'var' ) && ~isempty(D))
    isD = 1;
end

doDownsample = 0;
if( exist( 'dsFactorZ', 'var' ) && ~isempty(dsFactorZ))
	doDownsample = 1;
    downsampler = Tid.getDownsampler3dzRs();
end

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

if( isD )
    feats = cell( numClasses, 1);
else
    feats = [];
end

patches = cell( numClasses, 1);

xyz = cell( numClasses, 3 );

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
    xyz(c, 1:3) = { x, y, z };

    fprintf('  grabbing patches\n');
    [ patchesThis ] = grabPatchesSimple( double(datim), patchSize, [], {x,y,z}, mskim );
    
    if( doDownsample )
        patchesThis = (downsampler( patchesThis', patchSize, dsFactorZ ))';
    end
    
    patches{c}  = patchesThis;
    
    fprintf('  encoding\n');
    if( isD )
        feats{c} = encoding(  patches{c}', D, 'sc', [], param );
    end
    
    fprintf('  done!\n');
end
