% hrlr_concept

%% helper functions

downsampler = @(a,f) reshape( mean(reshape(a,f,[])), [], size(a,2));


%% 2d proof-of-concept

f = 3; % downsampling factor

sml = (2/12); % a small value
lrg = (5/12); % a large value

% two 9 x 1 basis elements ( excluding background )
elems = [ sml lrg sml lrg sml lrg sml lrg sml; ...
          lrg sml lrg sml lrg sml lrg sml lrg ]';

% probability of each element
elemProb = [ 0.5, 0.5 ];

% a test
% x = downsampler(elems,f);

%% the combinatorics of all possible hi-res patches 

hrPatchList = cell( 2*2*9, 1);
lrPatchList = cell( 2*2*9, 1);

k = 1;
for e = 1:2
    for dorot = 0:1
        for row  = 1:9
            hrPatchList{k} = toHiPatch( elems(:,e), row, dorot);
            lrPatchList{k} = downsampler(hrPatchList{k},f);
            k = k + 1;
        end
    end
end

%%

N = 3*length(lrPatchList);
X =  zeros( 9, N );
i = 1;
for k = 1:length(lrPatchList)
    for r = 1:3
        patch = lrPatchList{k}(r,:);
        X(:,i) = patch(:);
        i = i + 1;
    end
end

% params
param.K = 11;  % dictionary size
param.lambda=0.1;
param.numThreads=4; % number of threads
param.batchsize = N;
param.verbose=0;
param.iter = 250;  % num iterations.

D = mexTrainDL(X,param);
% this dictionary looks the way we'd expect - hooray

%% get the dictionary to the low-res space
Dd = downsampler(D,f);

%% dictionary learning round 2 (on 9x3 size patches )

M = ( 9 ) * length(lrPatchList ); 

Xdown = zeros( 3, M );
for i = 1:length(lrPatchList)
    lrPatch = lrPatchList{i};
    for j = 1:9
        Xdown(:,k) = lrPatch(:,j);
    end
end

alpha = mexLasso( Xdown, Dd, param);

% hallucinate high-res
Xhr = D * alpha 


%%
% imagesc( elems ); axis equal;, colormap gray;
% caxis([0 1]);


