% hrlr_concept

%% helper functions
downsampler = @(x,f) reshape( mean(reshape(permute(x,[3 1 2]),3,[])), size(x)./[1 1 f]);
downsamplerRs = @(X,sz,f) reshape(mean(reshape( X,[ sz size(X,2)]),3),[],size(X,2));
upsamplerNnRsZ = @(X, f) repmat( X, [f 1]) ;

%% 3d proof-of-concept

f = 3; % downsampling factor

N = 500; % number of random samples

[Dtrue, patchSize] = DictionarySymSampler.simpleDictionary3d();
patchSizeDs = patchSize ./ [1 1 f];

% rotRng = [1 2];
rotRng = [ ];
dss = DictionarySymSampler( Dtrue, [], patchSize, rotRng, [], [], [] );
X = dss.sample( 500 );

Xds = downsamplerRs( X, patchSize, f );     % downsample by averaging
Xdsus = upsamplerNnRsZ( Xds, f );           % upsample with NN


%% make a 2d dictionary to start with

tidDs = Tid( Xds, patchSizeDs, Tid.defaultParams() );
tidDs.buildDictionary();
tidDs.makeDictRotInv();

%% vis the 2d dictionary

dsize = size(tidDs.D,2);

for i = 1:dsize
    imdisp( reshape(tidDs.D(:,i), patchSizeDs ));
    pause;
    close all;
end


%% make a low res dictionary

tid = Tid( Xdsus, patchSize, Tid.defaultParams() );
tid.buildDictionary();
tid.makeDictRotInv();

%% visualize patches

dsize = size(tid.D,2);
reshapeSz = [ patchSize 1 ];

for i = 1:dsize
    imdisp(reshape(tid.D(:,i), reshapeSz([1 2 4 3]) ), ...
            'size', [1 patchSize(3)], ...
            'border', 0.05);
    pause;
    close all;
end
