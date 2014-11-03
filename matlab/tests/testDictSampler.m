% testDictSampler
%%

patch1 = zeros( 9, 9 );
patch1(5,:) = 1;

patch2 = zeros( 9, 9 );
patch2(5,:) = repmat([0.5 2 0.5], 1, 3);

patch3 = zeros( 9, 9 );
patch3(:,5) = 1;

patch4 = zeros( 9, 9 );
patch4(:,5) = repmat([0.5 2 0.5]', 3, 1);

patch5 = (1./3).*ones( 9, 9 );

patch6 = zeros( 9, 9 );
patch6(1,:) = 1;

patch7 = zeros( 9, 9 );
patch7(:,1) = 1;

patch8 = zeros( 9, 9 );
patch8(end,:) = 1;

patch9 = zeros( 9, 9 );
patch9(:,end) = 1;

patch10 = zeros( 9, 9 );
patch10(1:3,1) = 1;

patch11 = zeros( 9, 9 );
patch11(1:3,end) = 1;

patch12 = zeros( 9, 9 );

patch13 = zeros( 9, 9 );
patch13([32 41 50]) = 1;

patch14 = zeros( 9, 9 );
patch14([1 10 19]) = 1;
 
patch15 = zeros( 9, 9 );
patch15([5 14 23]) = 1;

patch16 = zeros( 9, 9 );
patch16(1:27) = 1/3;
%patch16(1:27) = 1;

patch17 = rand( 9, 9 );


D = [   patch1(:), patch2(:), patch3(:), patch4(:), patch5(:), patch6(:), patch7(:), ...
        patch8(:), patch9(:), patch10(:), patch11(:), patch12(:), patch13(:), patch14(:),...
        patch15(:), patch16(:), patch17(:)]';
    

%%
dsFactor = 3;
patchSize = [ 9 9 ] ;

patchParams = [ 3 1  1; ... %set
                3 4 12; ... %set
                3 7 12; ... %set
                2 1 15; ... %set
                2 4 15; ... %set
                2 7 15; ... %set
                1 1 12; ... %
                1 4 16; ... % questionable
                1 7 12 ];   %

patchParamsLeftOut = patchParams(2:end,:);

%%

[nodeTest, rootTest] = Dict2dTo3d.patchParamArrayToNode( patchParamsLeftOut );

%% 
clear d23;
d23 = Dict2dTo3d( D, patchSize(1), dsFactor );


dim = patchParams( 1, 1);
xyz = patchParams( 1, 2);
idxTrue = patchParams( 1, 3);

costs = d23.allPatchConfigCosts( dim, xyz, nodeTest );

%%

clear d23s;
d23s = Dict2dTo3dSampler( D, patchSize(1), dsFactor );

b = d23s.pc.constraintValueList( d23s.D2d, patchParams(:,1:2) );

%%

i = 1;
[ bestidx, sims, cmtx2, b2] = d23s.bestPatchConfig( b, i );


