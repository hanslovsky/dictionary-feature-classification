% hrlr_concept3

%% helper functions
% downsampler = @(x,f) reshape( mean(reshape(permute(x,[3 1 2]),3,[])), size(x)./[1 1 f]);
% downsamplerRs = @(X,sz,f) reshape(mean(reshape( X,[ sz size(X,2)]),3),[],size(X,2));
            
downsampler = @(X,f)(permute(reshape(mean(reshape( permute(x, [3 2 1]), f, [] )), sz([3 2 1])./[f 1 1]), [3 2 1]));
downsamplerRs = @(X,sz,f) (reshape(permute(reshape(mean(reshape(permute( reshape(X, [sz size(X,2)] ), [ 3 2 1 4]), f, [])), [sz(3)./f sz(1) sz(2) size(X,2)]), [ 3 2 1 4]), [], size(X,2) ));

%% 3d proof-of-concept ( segmentation )

f = 3; % downsampling factor

patchSize = [9 9 9];
segProbs = [ 0.5 0.5 ];
numDict   = 20;
patchSizeDs = patchSize ./ [1 1 f];
Dtrue = DictionarySymSampler.dctDictionary3d( patchSize, numDict );

% deal with a binary segmentation for now
seg = mnrnd( 1, segProbs, numDict ); 
seg = seg(:,1);

rotRng = [1 2];
dss = DictionarySymSampler( Dtrue, [], patchSize, rotRng, [], [], [] );
[X,di] = dss.sample( 1000 );

Xds = downsamplerRs( X, patchSize, f ); 

%%
tid = Tid( X, patchSize, Tid.defaultParams() );
tid.buildDictionary();
tid.makeDictRotInv();

tid.getFeatures();

%% visualize the dictionary (for fun)

% D = tid.D;
% N = size( D, 2);
% for i = 1:N
%    imdisp( permute( reshape( D(:,i), patchSize ), [ 1 2 4 3] ), ...
%        'border', 0.1 );
%    pause;
% end

%% do a low res dictionary

tidLR = Tid( Xds, patchSizeDs, Tid.defaultParams() );
tidLR.buildDictionary();

%%
% hypothesis - building a dictionary as described above should
% be better for segmentationt that a dictionary build directly from 
% the low-res data

% we have samples




