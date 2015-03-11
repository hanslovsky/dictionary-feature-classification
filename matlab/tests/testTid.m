% testTid

%% 
patchSize = [ 7 7 ];
dictPatchSize = [ 5 5 ];

X = rand( 49, 100 );

% X = rand( 49, 100 );

%%
patchSize = [ 7 7 ];
dictPatchSize = [ 5 5 ];

[tx,ty] = Tid.inferTranslations( patchSize, dictPatchSize )
[sp, mi] = Tid.translationSubPatches( patchSize, dictPatchSize );
sp 
mi

%% test dict helper

clear tid;
tid = Tid( X, patchSize, dictPatchSize );
tid.setComparator( 'ncc' );
tid.distTol = 5.0;

im1 = zeros( patchSize ); im1( 1:2, : ) = 1;
im2 = zeros( patchSize ); im2( 1:2, : ) = 1; im2( 3, : ) = 0.7; im2( 4, : ) = 0.3;
im3 = zeros( patchSize ); im3( 1:5, : ) = 1;

im1 = im1 + 0.1.*randn(size(im1));
im2 = im2 + 0.1.*randn(size(im2));
im3 = im3 + 0.1.*randn(size(im3));

% imdisp3d( cat( 3, im1, im2, im3 ));

D = [ reshape( im1, 1, [] ); reshape( im2, 1, [] ); reshape( im3, 1, [] ) ];

[ newD, sim ] = tid.translationInvariantHelper( D );
% sim;
sim(sim == 500) = 0;
aa = sum( sim, 3 )
sum(aa,1)

%% test get default dictionary function for debug
% clc;
% fun = Tid.getDictionaryBuildingFunction()