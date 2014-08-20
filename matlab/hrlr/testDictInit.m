% test initialization

% see if the initialization is smart enough to initialize with dis-similar
% patches

N = 3*length(lrPatchList);
X =  zeros( 9, N ); % num features x num observations
i = 1;
for k = 1:length(lrPatchList)
    for r = 1:3
        patch = lrPatchList{k}(r,:);
        X(:,i) = patch(:);
        i = i + 1;
    end
end

%%
X = zeros( 9, 20 );
N1 = 10
X(:, 1:N1) = repmat([0 0 0 1 1 1 0 0 0]', 1, N1 );
X(:, N1+1:end) = repmat([0 1 0 1 0 1 0 1 0]', 1, 20-N1 );


% params
param.K = 2;  % dictionary size
param.lambda=0.15;
param.numThreads=4; % number of threads
param.batchsize = N;
param.verbose=1;

% run one iteration so that things dont change much from the initialization
param.iter = 1;  

D = mexTrainDL(X,param);

D

%% conclusion:
% the dictionary is initialized randomly 
% and could perhaps benefit from a smater initialization