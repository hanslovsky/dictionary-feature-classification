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
labim = read_image_stack( labimfn );
mskim = read_image_stack( mskimfn );
mskim = mskim(:,:,:,1) & mskim(:,:,:,2) & mskim(:,:,:,3);

labim_bnd = uint8(round(mean(labim,4))) + 1;

clear labim;

numClasses = length(N);
labelList = unique(labim_bnd);

feats = cell( numClasses, 1);


for c = 1:numClasses
    
    label = labelList(c);
    fprintf('generating data for class (%d)\n', label);
    
    inds = find( mskim & (labim_bnd == label));
    numClass = length( inds );
    
    r = randperm( numClass, N(c));
    inds = inds(r);
    
    [x,y,z] = ind2sub( size(labim_bnd), inds);
    fprintf('  grabbing patches\n');
    patches = grabPatchesSimple( double(datim), patchSize, [], {x,y,z} );
    fprintf('  encoding\n');
    feats{c} = encoding(  patches', D, 'sc', 0, param );
    fprintf('  done!\n');
end