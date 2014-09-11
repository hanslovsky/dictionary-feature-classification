% hrlr_concept

%% helper functions
downsampler = @(x,f) reshape( mean(reshape(permute(x,[3 1 2]),3,[])), size(x)./[1 1 f]);
downsamplerRs = @(X,sz,f) reshape(mean(reshape( X,[ sz size(X,2)]),3),[],size(X,2));

%% 3d proof-of-concept

f = 3; % downsampling factor

patchSize = [3 3 3];
patchSizeDs = patchSize ./ [1 1 f];
Dtrue = DictionarySymSampler.simpleDictionary3d();

rotRng = [1 2];
dss = DictionarySymSampler( Dtrue, [], patchSize, rotRng, [], [], [] );
X = dss.sample( 1000 );

Xds = downsamplerRs( X, patchSize, f ); 

tid = Tid( X, patchSize, Tid.defaultParams() );
tid.buildDictionary();
tid.makeDictRotInv();


