% testDictSamplerSub

%% params

dsFactor  = 3;
patchSize = [ 9 9 ];
patchSizeHR = [ patchSize patchSize(1)];
patchSizeLR = [ patchSize patchSize(1)./dsFactor];

overlappingPatches = 1;

numDict = 20;

%% random dictionary
% 
% clear d23;
% 
% % D = rand( numDict, prod(patchSize));
% 
% % normalize rows of D
% D = bsxfun( @rdivide, D, sqrt(sum( D.^2, 2 ) ));
% 
% d23 = Dict2dTo3dSamplerSub( D, patchSize(1), dsFactor, overlappingPatches, 0 );
% d23.pc
% d23.pc.dimXyzList
% 
% x = rand( prod(patchSizeLR), 1);
% x = x ./ norm(x);
% 
% % [ patchParams, modelList, pv, patch ] = d23.solveHR( x );

%% gabor dictionary

Dcell = gaborFilterBank( 5, 8, patchSize(1), patchSize(1) );
D = zeros( numel(Dcell), prod(patchSize) ); 
for i = 1:numel(Dcell)
    D(i,:) = real(Dcell{i}(:));
end

D = bsxfun( @rdivide, D, sqrt(sum( D.^2, 2 ) ));

clear d23;
d23 = Dict2dTo3dSamplerSub( D, patchSize(1), dsFactor, overlappingPatches, 0 );

x = D(5,:) + 0.01.*randn( size(D(7,:)));
x = repmat( reshape( x, 9, 9 ), [1 1 3]);
% sc(reshape( x, 9, 9));

%%

[ patchParams, modelList, pv, patch ] = d23.solveHR( x(:), K );

% imdisp3d( reshape( pv, d23.pc.sz3d))

%%

i = 1;
K = 4;

% for i = 1 : d23.pc.numLocs
for i = 1 : 3
    [ dictIdxs, dictCosts, models ] = d23.bestKdicts(  x(:), i, K );

    dictIdxs 
    % dictCosts
    
end

%% test solveHR


%% test bestKdicts on HR observation

% i = 7;
% x = rand( 9*9*9, 1 );
% 
% % this, i, x, doBest, returnList )
% [ dictIdxs, dictCosts, models ] = d23.fitIdxAndModel( i, x(:), 0, 1 );
% dictIdxs

%%
% pc = d23.pc
% 
% % msk = d23.pc.planeMaskI( 1 );
% msk = d23.pc.planeMaskLRI( 1 );
% 
% % this = d23;
% % centered = false;
% % 
% % xyz = 1;
% % n   = 3;
% % f   = 3;

%%

