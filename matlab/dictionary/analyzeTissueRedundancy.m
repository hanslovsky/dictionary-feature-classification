% visDictionaryH5
%%

downsamplerRs = Tid.getDownsampler3dzRs();
normRows = @(x)( x./repmat( sqrt(sum( x.*x, 2 )), 1, size(x,2)) );

%%

N = 4000;
% patchSize  = [9 9 9];
patchSize  = [15 15 15];
dsfactor = 3;
patchSizeD = patchSize ./ [1 1 dsfactor];

donorm = 1;

%% grab some patches

training = 18;
test = 3;
ds = ds_medulla(training, test);


%% load the image
im = read_image_stack( ds.data_fn{1} );
lb = read_image_stack( ds.labels_fn{1} );
lb = (mean( lb , 4 ) > 0.5);

msk = read_image_stack( ds.mask_fn{1} );
msk = max( msk, [], 4 );

%% grab patches

% patches
p  = grabPatchesSimple( im, patchSize, N, [], msk );

% downsample the patches
pd =  (downsamplerRs(p',patchSize,dsfactor))';

if( donorm )
   disp('normalizing');
   p = normRows( p );
   pd = normRows( pd );
end

%% pairwise distances 

disp('HR pairwise distances');
d  = pdist( p );

disp('LR pairwise distances');
dd = pdist( pd );

disp('done');

%% analyze pairwise distances 

% joint histogram
jointHist( [d' dd'], 128 );

%% grab patches based on boundary prediction

% patches
ppos  = grabPatchesSimple( im, patchSize, N, [], (msk &  lb) );
pneg  = grabPatchesSimple( im, patchSize, N, [], (msk & ~lb) );

% downsample the patches
pposd =  (downsamplerRs( ppos', patchSize, dsfactor))';
pnegd =  (downsamplerRs( pneg', patchSize, dsfactor))';

if( donorm )
   disp('normalizing');
   ppos = normRows( ppos );
   pneg = normRows( pneg );
   
   pposd = normRows( pposd );
   pnegd = normRows( pnegd );   
end

%% pairwise distances between boundary and non-boundary patches

disp('HR pairwise boundary prediction distances');
dseg  = pdist2(  ppos,  pneg );

disp('LR pairwise boundary prediction distances');
dsegd = pdist2( pposd, pnegd );

%% joint histogram

jointHist( [dseg(:) dsegd(:)], 128 );

%% visualize most similar patches with different labels

[~,k]=sort(dsegd(:));
[i,j] = ind2sub( size(dsegd), k);

%%

for n = 1:50
    n
    
    norm( ppos(i(n),:) - pneg(j(n),:))
    i(n)
    j(n)
    figure;
    imdisp( permute(reshape( ppos(i(n),:), patchSize ), [1 2 4 3]), 'border', 0.1 );
    
    figure;
    imdisp( permute(reshape( pneg(j(n),:), patchSize ), [1 2 4 3]), 'border', 0.1 );
    
    pause;
    close all;

end

