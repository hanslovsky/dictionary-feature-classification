% analyzeTissueRedundancy
%%

downsamplerRs = Tid.getDownsampler3dzRs();
normRows = @(x)( x./repmat( sqrt(sum( x.*x, 2 )), 1, size(x,2)) );

%%

N = 8000;
patchSize  = [9 9 9];
% patchSize  = [15 15 15];

dsfactor = 3;
% dsfactor = 5;

patchSizeD = patchSize ./ [1 1 dsfactor];

donorm = 0;
distMeasure = 'correlation';

destdir='/data-ssd1/john/projects/imageAnalysis/interestingPatchPairs/9-9-9';
patchPointFile = '/groups/saalfeld/home/bogovicj/projects/dictionary/imageTissueAnalysis/patchCentersSimilarAtLowRes_9-9-9.csv';

%% grab some patches

training = 18;
test = 3;
ds = ds_medulla(training, test);


%% load the image
im = read_image_stack( ds.data_fn{1} );
lb = read_image_stack( ds.labels_fn{1} );
lb = (mean( lb , 4 ) > 0.5);
lb = uint8 ( lb );

msk = read_image_stack( ds.mask_fn{1} );
msk = max( msk, [], 4 );

%% grab patches
% 
% % patches
% p  = grabPatchesSimple( im, patchSize, N, [], msk );
% 
% % downsample the patches
% pd =  (downsamplerRs(p',patchSize,dsfactor))';
% 
% if( donorm )
%    disp('normalizing');
%    p = normRows( p );
%    pd = normRows( pd );
% end
% 
% %% pairwise distances 
% 
% disp('HR pairwise distances');
% d  = pdist( p, distMeasure );
% 
% disp('LR pairwise distances');
% dd = pdist( pd, distMeasure );
% 
% disp('done');
% 
% %% analyze pairwise distances 
% 
% % joint histogram
% jointHist( [d' dd'], 128 );

%% grab patches based on boundary prediction

% patches
[ppos,coordsPos]  = grabPatchesSimple( im, patchSize, N, [], (msk &  lb) );
[pneg,coordsNeg]  = grabPatchesSimple( im, patchSize, N, [], (msk & ~lb) );

% grab specific patches
% pnew = retrieveMatchPatches( patchPointFile, im, patchSize );

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

%% intensity variation within a patch

pposvar= var(ppos,0,2);
pnegvar= var(pneg,0,2);

% figure; hist( prod(patchSize).*pposvar, 128 );
% figure; hist( prod(patchSize).*pnegvar, 128 );

%% pairwise distances between boundary and non-boundary patches

disp('HR pairwise boundary prediction distances');
dseg  = pdist2(  ppos,  pneg, distMeasure );

disp('LR pairwise boundary prediction distances');
dsegd = pdist2( pposd, pnegd, distMeasure );
disp('done');

%% visualize patches that are similar at low res but high patch variance

num = 500;
[~,k]=sort(dsegd(:));
[i,j] = ind2sub( size(dsegd), k);
k = k(1:num);
i = i(1:num);
j = j(1:num);

% figure; hist( dsegd(k) );
vdown = prod(patchSize).*( pposvar(i) + pnegvar(j) );

[vdownSort,l] = sort( vdown, 'descend' );

for n = 1:num
    m  = l(n);
    
    cpos = { coordsPos{1}(i(m)), coordsPos{2}(i(m)), coordsPos{3}(i(m)) };
    mskpatchpos = grabPatchesSimple( msk, patchSize, [], cpos );
    lbpatchpos  = grabPatchesSimple(  lb, patchSize, [], cpos );
    
    cneg = { coordsNeg{1}(j(m)), coordsNeg{2}(j(m)), coordsNeg{3}(j(m)) };
    mskpatchneg = grabPatchesSimple( msk, patchSize, [], cneg );
    lbpatchneg  = grabPatchesSimple(  lb, patchSize, [], cneg );

    
    hppos = figure;
    imdisp( permute(reshape( ppos(i(m),:), patchSize ), [1 2 4 3]), 'border', 0.1 );
%     figure;
%     imdisp( permute(reshape( mskpatchpos, patchSize ), [1 2 4 3]), 'border', 0.1 );
    lbpos = figure;
    imdisp( permute(reshape( lbpatchpos, patchSize ), [1 2 4 3]), 'border', 0.1 );
    
    hpneg = figure;
    imdisp( permute(reshape( pneg(j(m),:), patchSize ), [1 2 4 3]), 'border', 0.1 );
%     figure;
%     imdisp( permute(reshape( mskpatchneg, patchSize ), [1 2 4 3]), 'border', 0.1 );
    lbneg = figure;
    imdisp( permute(reshape( lbpatchneg, patchSize ), [1 2 4 3]), 'border', 0.1 );
    
    hpposd = figure;
    imdisp( permute(reshape( pposd(i(m),:), patchSizeD ), [1 2 4 3]), 'border', 0.1 );
    
    hpnegd = figure; 
    imdisp( permute(reshape( pnegd(j(m),:), patchSizeD ), [1 2 4 3]), 'border', 0.1 );
    
    [ cpos{:} cneg{:} ]
    
%     pause;
    dosave = input( 'save to file? ' );
    if(  dosave > 0 )
        fprintf('saving...')
       figure(hppos);
       export_fig(sprintf('%s%sposPatch_id%d_pt_%d-%d-%d', destdir, filesep, n, cpos{1}, cpos{2}, cpos{3} )); 
       figure(hpneg);
       export_fig(sprintf('%s%snegPatch_id%d_pt_%d-%d-%d', destdir, filesep, n, cneg{1}, cneg{2}, cneg{3} )); 
       fprintf('done!\n');
    end
    close all;
    
end

%%
% n = 1;
% dsegd(k(n))
% norm( pposd(i(n),:) - pnegd(j(n),:))

%% joint histogram

jointHist( [dseg(:) dsegd(:)], 128 );

%% visualize most similar patches with different labels

[~,k]=sort(dsegd(:));
[i,j] = ind2sub( size(dsegd), k);

for n = 1:50
%     n
%     norm( ppos(i(n),:) - pneg(j(n),:))
%     i(n)
%     j(n)
    
    figure;
    imdisp( permute(reshape( ppos(i(n),:), patchSize ), [1 2 4 3]), 'border', 0.1 );
    
    figure;
    imdisp( permute(reshape( pneg(j(n),:), patchSize ), [1 2 4 3]), 'border', 0.1 );
    
    pause;
    close all;

end

