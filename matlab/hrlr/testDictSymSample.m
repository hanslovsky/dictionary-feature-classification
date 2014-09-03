% testDictSymSample


patchSize = [3 3 3];

%D(:,1) = reshape(eye(3),[],1);
%D(:,2) = reshape(ones(3,3),[],1);
%D(:,3) = reshape([1 1 1; 0 0 0; 0 0 0],[],1);

%D = DictionarySymSampler.standardDictionary2d();
Dtrue = DictionarySymSampler.simpleDictionary3d();

[xs,ys,zs] = ndgrid( -1:1, -1:1, -1:1 );
shiftRng = [ xs(:), ys(:), zs(:) ];

rotRng = [1 2];

dss = DictionarySymSampler( Dtrue, [], patchSize, rotRng, [], [], [] );

%dictProbs = zeros( size(Dtrue,2), 1);
%dss = DictionarySymSampler( Dtrue, dictProbs, patchSize, rotRng, [], [], [] );

%dss = DictionarySymSampler( Dtrue, [], patchSize, [], [], shiftRng );

X = dss.sample( 1000 );

 
clear param
param.K = 25;  % dictionary size
param.lambda=0.15;
param.numThreads=4; % number of threads
param.batchsize = 100;
param.verbose=1;
param.iter = 100;  

%D = mexTrainDL( X, param );

%% try Tid

tid = Tid( X, patchSize, param ); 
tid.workSomeMagic();
Dx = tid.xfmDict();


%% older tests of shift / rotation xfm's
%
%a = D(:,1);
%b = reshape( a, patchSize );
%c = padarray( b, [1,1,1] );
%
%si = shiftToIdx( [-1 0 0], patchSize );
%T = [ -1 0 0; 0 1 0; 0 0 1 ];
%ti = xfmToIdx( T, patchSize );
%
%c( si(ti))

%%

K = param.K;
distTol = 0.1;

numInvalid = 0;
dxsz  = size(Dx,2); % size of transformed dictionary
rxi= cell2mat(tid.rotXfmIdx);

invalid =  zeros(1,K);

allIdxs = 1:dxsz;
invariantIdxs = true(1,dxsz);

for i = allIdxs
    % check if this index is still valid
    % it may have been removed at a previous iteration
    if( ~invariantIdxs(i) )
        fprintf('skipping %d\n',i);
        continue;
    end
    
    j = setdiff( allIdxs, i );
    
    Dxfm = Dx( :, i );
    Dxfm = Dxfm( rxi );
    Draw = Dx( :, j );
    
    pd = min( pdist2( Dxfm', Draw' ), [], 1);
    similar = ( pd < distTol );
    
    % update indices of dict elems that are 
    % transformations of this dict element
    invariantIdxs( allIdxs(j(similar)) ) = false;
    
end

nnz(invariantIdxs)

%%

% Dx = Dx(:, invariantIdxs);
% size(Dx)


