% visDictionaryH5

%%

fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0068_learnDictSpams_20140813133126/im_normalized_0mean_patches.h5'
destdir = '/data-ssd1/john/projects/dictionary/vis/exp0044_patches';

fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0069_learnDictSpams_20140813163922/im_normalized_0mean_patches.h5'
destdir = '/data-ssd1/john/projects/dictionary/vis/exp0044_patches_downavg';


num = 25;
patchSize  = [9 9 9];
patchSizeD = [3 9 9];

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

%%

if( ~exist( destdir, 'dir' ))
    mkdir( destdir );
end

for i = randperm( numPatches, num )
% for i = 750
    i
    
    for j = 1:length(suffixes)
        j
        patch = reshape(pList{j}(i,:),  szList{j} );
        
        sc( permute( patch, permList{j}));
        pause(0.2);
        export_fig( fullfile(destdir,sprintf('patch_%04d_%s.png',i, suffixes{j})) );
        close all;
    end
    
    % do difference 
    diff = reshape( pList{1}(i,:) - pList{3}(i,:),  szList{j} );
    sc( permute( diff, [2 3 4 1]), [-1 1]);
    pause(0.2);
    export_fig( fullfile(destdir,sprintf('patch_%04d_diff.png',i )) );
    close all;
end

%% testing stuff out


% patchOrig = reshape( pOrig(1,:), patchSize );
% figure; sc( permute( patchOrig, [2 3 4 1]));


%%
% 
% fn = '/groups/saalfeld/home/bogovicj/tmp/testPatches.h5';
% 
% pOrig = h5read( fn, '/patchesFlat' );
% [numPatches, patchElem] = size( pOrig );
% pDown  = h5read( fn, '/patchesDownFlat' );
% pRecon = h5read( fn, '/patchesUpFlat' );
% 
% pOrigR = h5read( fn, '/patches' );
% pDownR  = h5read( fn, '/patchesDown' );
% pReconR = h5read( fn, '/patchesUp' );

%%

% p = permute( reshape( pOrig, [9 9 9]), [2 3 1])
% permute( pOrigR, [2 3 1])
% 
% p = permute( reshape( pDown, [3 9 9]), [2 3 1])
% permute( pDownR, [2 3 1])

%%

% patchDown = reshape( pDown(1,:), [9 9 3]);
% figure; sc( permute( patchDown, [1 2 4 3]));

% patchDown = reshape( pDown(1,:), [9 3 9]);
% size( permute( patchDown, [1 3 4 2]) )
% figure; sc( permute( patchDown, [1 3 4 2]));

% patchDown = reshape( pDown(1,:), [3 9 9]);
% size( permute( patchDown,[2 3 4 1]) )
% figure; sc( permute( patchDown, [2 3 4 1]));