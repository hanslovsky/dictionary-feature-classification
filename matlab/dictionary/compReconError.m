% visDictionaryH5

%%

% fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0068_learnDictSpams_20140813133126/im_normalized_0mean_patches.h5'
% destdir = '/data-ssd1/john/projects/dictionary/vis/exp0044_patches';

% fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0069_learnDictSpams_20140813163922/im_normalized_0mean_patches.h5'
% destdir = '/data-ssd1/john/projects/dictionary/vis/exp0044_patches_downavg';

fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0071_learnDictSpams_20140814105705/im_normalized_0mean_patches.h5'
destdir = '/data-ssd1/john/projects/dictionary/vis/exp0044_patches_downavg2';

num = 25;
patchSize  = [9 9 9];
patchSizeD = [3 9 9];
ds=3;

set='testing'
suffixes ={'orig','down','recon'};

%% read the dictionary

pOrig = h5read( fn, sprintf('/%s/patches',set));
[numPatches, patchElem] = size( pOrig );

pDown  = h5read( fn, sprintf('/%s/patchesDown',  set));
pRecon = h5read( fn, sprintf('/%s/patchesRecon', set));

pList =  {pOrig, pDown, pRecon};
szList = { patchSize, patchSizeD, patchSize };
% permList = {[2 3 4 1]; [3 2 4 1]; [2 3 4 1]};
permList = {[2 3 4 1]; [2 3 4 1]; [2 3 4 1]};
downPerm = [2 3 1];


%%

if( ~exist( destdir, 'dir' ))
    mkdir( destdir );
end

writeImages = 0;

ssdRecon = zeros(numPatches,1);
ssdIntrp = zeros(numPatches,1);

% for i = randperm( numPatches, num )
for i = 816
% for i = 1:numPatches
%     i
    fprintf('patch %d of %d\n',i,numPatches);
    
    % estimate error for NN interpolation
    step = 1/ds; 
    half = (ds-1)/2;
    [xq,yq,zq] = meshgrid( 1:patchSize(1), 1:patchSize(2), ...
                            1-(half*step) : step : patchSize(3)/ds + (half*step));

    patchDown = permute( reshape(pList{2}(i,:),  szList{2}), downPerm );
    patchUp = interp3(patchDown, xq, yq, round(zq), 'nearest');
    patchUpFlat = permute(patchUp, [3 1 2]);

    % output SSD for upsampled and dictionary recon'd patches
    ssdRecon(i) = sum((pList{3}(i,:) - pList{1}(i,:)).^2);
    ssdIntrp(i) = sum((pList{3}(i,:) - patchUpFlat(:)').^2);
    
end
ssdRecon(i)
ssdIntrp(i)

meanSsdRecon = mean(ssdRecon);
meanSsdIntrp = mean(ssdIntrp);

fprintf('mean recon SSD: %f\n', meanSsdRecon );
fprintf('mean intrp SSD: %f\n', meanSsdIntrp );

fprintf('median recon SSD: %f\n', median(ssdRecon) );
fprintf('median intrp SSD: %f\n', median(ssdIntrp) );
