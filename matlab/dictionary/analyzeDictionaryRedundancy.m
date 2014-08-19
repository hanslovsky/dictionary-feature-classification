% visDictionaryH5


%%
% fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0043_learnDictSpams_20140804162831/im_normalized_0mean_dict.h5';
% destdir = '/data-ssd1/john/projects/dictionary/vis/exp0043';

% fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0038_learnDictSpams_20140804162826/im_normalized_0mean_dict.h5'
% destdir = '/data-ssd1/john/projects/dictionary/vis/exp0038';

fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_sh_exp/exp0075_learnDictSpams_20140814142620/im_normalized_0mean_dict.h5';
destdirbase = '/groups/saalfeld/home/bogovicj/projects/dictionary/analysis/exp0075';

if( ~exist( destdirbase, 'dir' ))
    mkdir( destdirbase );
end

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

dat = [distLow', distHi'];
maxval = max(dat(:));
minval = min(dat(:));
% minmax = prctile(  dat(:), [1 99]);
% minval = minmax(1);
% maxval = minmax(2);

rng = linspace( minval, maxval, 128 );

% [Nlo, Clo ] = hist( distLow, rng );
% [Nhi, Chi ] = hist(  distHi, rng );


[N,C] = hist3( [distHi', distLow'], { rng, rng } );

cmap = diverging_map(linspace(0,1,256), [0.2 0.2 0.7], [0.7 0.2 0.2]);

distLowSq = squareform( distLow );
distHiSq = squareform( distHi );

make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

%%

figure('color', 'w', 'position', [1280 563 871 890] );
% subplot(2,2,1);
% sc(N, cmap);

subplot(2,2,2);
area( rng, sum(N,1) ); axis square;
axtmp = axis;
axis( [minval maxval axtmp(3) axtmp(4)]);
xlabel( 'low res patch pair distance');
ylabel( 'Count');

subplot(2,2,3);
area( sum(N,2), rng); axis square;
set( gca, 'xdir', 'reverse');
axtmp = axis;
axis( [axtmp(1) axtmp(2) minval maxval  ]);
xlabel( 'Count');


subplot(2,2,4);
sc(N, cmap); axis square;
axis xy;
ylabel( 'high res patch pair distance');


% export_fig( fullfile( destdirbase, 'patchDistanceHiLoRes_jointHist.png') );

%% visualize something

% i = find( (distLowSq < 0.25) & (distLowSq > 0) );

% i = find( (distLowSq > 0.25) & (distLowSq < 0.4) & (distLowSq > 0) );
% jj = find( distHiSq(i) < 0.75 & distHiSq(i) > 0.5); 

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


%% visualize something

i = find( (distHiSq > 1.4) & (distHiSq < 1.43) );
size(i)


prctileVals = prctile( distLowSq(i), [2 98 ]);
% jj = find( distLowSq(i) > prctileVals(2));
jj = find( distLowSq(i) < prctileVals(1));

size(jj)

%%

% t = 250;
% t = 556;
% t = 566;
% t = 760;

% for t = 1:length(jj)
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
    
    fprintf( 'patch match( %d, %d )\n', n, m );
    fprintf( 'dist at low res:\t%f\ndist at hi res: \t%f\n', distLowSq(n,m), distHiSq(n,m) );
    figure;
    imdisp( permute( cat(3,patch_m,patch_n), [1 2 4 3]), ...
            'Size', [6 3], 'Border', 0.01);
    
%     figure;
%     imdisp( permute( cat(3, patchLo_m, patchLo_n ), [1 2 4 3]), ...
%             'Size', [2 3], 'Border', 0.01);
    pause;
    close all;

%     export_fig( ...
%         fullfile('/groups/saalfeld/home/bogovicj/projects/dictionary/analysis', ...
%                 sprintf('dictElems_hiFarLoClose_%d_%d-%d.png',t,n,m)));

end

%% test

% D = [ (0:35)', (36:71)' ]'
% D = reshape( 0:71, 2, []);
% dsfactor = 3;
% 
% [dictElem, patchElem] = size( D );
% tmp = reshapePy( D, [], 3, dictElem )

