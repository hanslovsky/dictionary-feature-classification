% visDictionaryH5

%%
% fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0043_learnDictSpams_20140804162831/im_normalized_0mean_dict.h5';
% destdir = '/data-ssd1/john/projects/dictionary/vis/exp0043';

% fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0038_learnDictSpams_20140804162826/im_normalized_0mean_dict.h5'
% destdir = '/data-ssd1/john/projects/dictionary/vis/exp0038';

fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0068_learnDictSpams_20140813133126/im_normalized_0mean_patches.h5'
destdir = '/data-ssd1/john/projects/dictionary/vis/exp0044_patches';

num = 50;
patchSize = [9 9 9];

%% read the dictionary

D = h5read( fn, '/dict');
[dictElem, patchElem] = size( D );

%%

if( ~exist( destdir, 'dir' ))
    mkdir( destdir );
end

for i = randperm( dictElem, num )
    i
    patch = reshape(D(i,:), patchSize);
    sc( permute( patch, [1 2 4 3]));
%     col orbar;
    pause(0.2);
    export_fig( fullfile(destdir,sprintf('dict_%04d.png',i)) );
    close all;
end

%%

