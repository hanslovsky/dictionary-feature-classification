% visDictionaryH5

%%
fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0043_learnDictSpams_20140804162831/im_normalized_0mean_dict.h5';
% destdir = '/data-ssd1/john/projects/dictionary/vis/exp0043';

% fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0038_learnDictSpams_20140804162826/im_normalized_0mean_dict.h5'
% destdir = '/data-ssd1/john/projects/dictionary/vis/exp0038';

% fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0068_learnDictSpams_20140813133126/im_normalized_0mean_patches.h5'
% destdir = '/data-ssd1/john/projects/dictionary/vis/exp0044_patches';

num = 50;
patchSize  = [9 9 9];
patchSizeD = [9 9 3];
dsfactor = 3;


%% read the dictionary

D = h5read( fn, '/dict');
[dictElem, patchElem] = size( D );

%%

tmp  = mean(reshapePy( D, [], dsfactor, dictElem ), 2);
Dlow = reshape( tmp, [], dictElem )';

%%

distLow = pdist(Dlow);
distHi  = pdist(D);

[N,C] = hist3([distLow', distHi'], [128,128]);

cmap = diverging_map(linspace(0,1,256), [0.2 0.2 0.7], [0.7 0.2 0.2]);
sc(N, cmap);

distLowSq = squareform( distLow );
distHiSq = squareform( distHi );

%% visualize something

% i = find( (distLowSq < 0.25) & (distLowSq > 0) );
i = find( (distLowSq > 0.25) & (distLowSq < 0.4) & (distLowSq > 0) );

jj = find( distHiSq(i) < 0.75 & distHiSq(i) > 0.5); 
size(jj)

% [maxval,j] = max( distHiSq(i) );
% maxval

for t = 1:length(jj)
    t
    j = jj(t);
    [n,m] = ind2sub( size( distHiSq ), i(j) );
    
    % patch n
    patch_n = reshapePy(D(n,:), patchSize);
    patch_m = reshapePy(D(m,:), patchSize);
    
    patchLo_n = reshapePy(Dlow(n,:), patchSizeD);
    patchLo_m = reshapePy(Dlow(m,:), patchSizeD);
    
%     figure; 
%     sc( permute( patch_n, [1 2 4 3]));
%     figure; 
%     sc( permute( patch_m, [1 2 4 3]));
    
    fprintf( 'dist at low res:\t%f\ndist at hi res: \t%f\n', distLowSq(n,m), distHiSq(n,m) );
    figure;
    imdisp( permute( cat(3,patch_m,patch_n), [1 2 4 3]), ...
            'Size', [6 3], 'Border', 0.01);
    
    figure;
    imdisp( permute( cat(3, patchLo_m, patchLo_n ), [1 2 4 3]), ...
            'Size', [2 3], 'Border', 0.01);
    pause;
    close all;

end

%% test

% D = [ (0:35)', (36:71)' ]'
% D = reshape( 0:71, 2, []);
% dsfactor = 3;
% 
% [dictElem, patchElem] = size( D );
% tmp = reshapePy( D, [], 3, dictElem )

