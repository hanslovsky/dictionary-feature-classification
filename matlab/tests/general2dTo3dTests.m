% general2dTo3dTests

%% patches for testing

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
patch16(1:27) = 1;

D = [   patch1(:), patch2(:), patch3(:), patch4(:), patch5(:), patch6(:), patch7(:), ...
        patch8(:), patch9(:), patch10(:), patch11(:), patch12(:), patch13(:), patch14(:),...
        patch15(:), patch16(:)]';

%% generate a known patch from consistent constraints

patchParams = [ 3 1  1; ... %set
                3 4 12; ... %set
                3 7 12; ... %set
                2 1 15; ... %set
                2 4 15; ... %set
                2 7 15; ... %set
                1 1 12; ... %
                1 4 16; ... %
                1 7 12 ];   %

[rc, lc] = rootConstraintFromArray( patchParams ); 

rc.printAsTree()

dsFactor = 3;
patchSize = [ 9 9 ] ;

clear d23;
% d23 = Dict2dTo3d( D, patchSize(1), dsFactor );
d23 = Dict2dTo3dConstr( D, patchSize(1), dsFactor );

[ pv, patch, cmtx, b ] = d23.patchFromParams( lc );
imdisp3d( patch )



%%
f = 3;

dim = 1;
xyz = 4;
patch = patch16;

% xyz dim
pm1 = Dict2dTo3d.planeMaskF( [9 9 9], xyz, dim, f ); 
imdisp3d( pm1 );

% dim xyz
I2 = Dict2dTo3d.fill3dWith2d( [9 9 9], dim, xyz, f, patch );
figure; imdisp3d( I2 );

%% what is the cost of this patch? 

% expect it to be correct?

p1 = squeeze(patch(5,:,:));
dim1 = 1;
xyz1 = 4;

% p2 = patch5;
% dim2 = 1;
% xyz2 = 4;

p2 = patch5;
dim2 = 1;
xyz2 = 4;

[ sim, x, cmtx, b, pm1, pm2, overlap ] = d23.patchSimilarity( p1, xyz1, dim1, p2, xyz2, dim2 );


