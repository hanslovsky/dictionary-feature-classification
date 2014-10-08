%% hrlr_samplingOpt_tests

import net.imglib2.algorithms.patch.SubPatch2dLocation;
import net.imglib2.algorithms.opt.astar.AStarMax;
import java.util.*;

%%  helper functions

downsamplerRs = Tid.getDownsampler3dzRs();
summer2Dx = Tid.sum2dx();
summer2Dy = Tid.sum2dy();
summer2Dxy = Tid.sum2dxy();

normRows = @(x)( x./repmat( sqrt(sum( x.*x, 2 )), 1, size(x,2)) );


%%

% patch2d = [ 0 0 0; 1 1 1; 0 0 0 ]
patch1 = zeros( 9, 9 );
patch1(5,:) = 1;

patch2 = zeros( 9, 9 );
patch2(5,:) = repmat([0.5 2 0.5], 1, 3);

patch3 = zeros( 9, 9 );
patch3(:,5) = 1;

patch4 = zeros( 9, 9 );
patch4(:,5) = repmat([0.5 2 0.5]', 3, 1);

patch5 = (1./3).*ones( 9, 9 );

X = [patch1(:), patch2(:), patch3(:), patch4(:), patch5(:)]';

%%

% % peasy = ones( 9, 9 );
% 
% clear d23;
% sz = 9;
% f  = 3;
% 
% xyz = 5;
% n1 = 1;
% n2 = 2;
% 
% d23 = Dict2dTo3d( X, sz, f );
% d23
% 
% d23.buildConstraintMatrix();
% [psim, x, cmtx, b, pm1, pm2, overlap] = d23.patchSimilarity( patch5, xyz, n1, patch2, xyz, n2 );
% 
% psim
% 
% values = zeros( size( overlap ));
% 
% values( overlap ) = x;
% 
% imdisp( permute(values, [1 2 4 3]), 'border', 0.02 );

%%
% clear d23
%dbstatus
% sz = 9;
% f  = 3;
% d23 = Dict2dTo3d( X, sz, f );
% d23
% 
% d23.sample2dTo3d()

%%

clear d23;

sz = 9;
f  = 3;
d23 = Dict2dTo3d( X, sz, f );

best = d23.testOutThings();

root = d23.p2dFill3d.getRoot();

%% test data for constraints 
% sz = [ 9 9 9 ];
% 
% xyz1 = [ 5 5 5 ];
% n1   = 1;
% 
% f = 3;
% szSm = [ 9 9 ] ./ f;
% 
% patchIdxs = zeros( [ 3 max(szSm) ] );
% 
% patchIdxs(1,1) = 4;
% patchIdxs(3,2) = 3
% 
% constraints = d23.getSumConstraints( patchIdxs )
% 
% d   = 2;
% xyz = 5;
% cc = d23.collectConstraints( constraints, 2, d )
% 
% sim = Dict2dTo3d.patchConsistencyConstraints( cc, X, [ 9 9 ], f )


%%
% sz3 = [ 9 9 9 ];
% xyz1 = 5;
% xyz2 = 5;
% n1 = 1;
% n2 = 3;
% 
% pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, n1, f );
% pm2 = Dict2dTo3d.planeMaskF( sz3, xyz2, n2, f );
% 
% overlap = (pm1 > 0)  & (pm2 > 0 );


% overlap = Dict2dTo3d.planeIntersections( sz3, xyz1, n1, xyz2, n2 );
% 
% p1Line = Dict2dTo3d.patchIntersectionFromMask( patch1(:), xyz1, n1, overlap );
% p2Line = Dict2dTo3d.patchIntersectionFromMask( patch2(:), xyz2, n2, overlap );

%%

% pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, n1, f );
% pm1(overlap)
% 
% subpatch = patch3( unique(pm1( overlap )))
% 
% nnz( subpatch )

%% test data for constraints

% sz = [ 9 9 9 ];
% f = 3;
% 
% patchIdxs = zeros( [ 3 max(sz) ] );
% patchIdxs(1,1) = 4;
% patchIdxs(3,2) = 3


%%

% d23.buildConstraintMatrix();
% psim = d23.patchSimilarity( patch1, xyz, d, patch2, xyz, 1 );
% 
% psim

%%


% spl1 = SubPatch2dLocation( 1, 1, 10.0);
% spl2 = SubPatch2dLocation( 2, 1, 5.0);
% 
% splCandidateList = ArrayList( );
% splCandidateList.add( spl2 );
%     
% as = AStarMax( spl1 );
% as.addCandidates( splCandidateList );

%%

% d = 1;
% xyz = 1;
% switch d
%     case 1
%         patchSums( Dict2dTo3d.planeMask( [3 3 3], xyz, d )) = summer2Dxy( patch1, 3 );
%     case 2
%         val = xyz(2);
%     case 3
%         val = xyz(3);
%     otherwise
%         error('invalid normal direction');
% end


%%
% 
% patch1 = zeros( 3, 3 );
% patch1(2,:) = 1;
% 
% patch2 = zeros( 3, 3 );
% patch2(2,:) = [0.5 2 0.5];
% 
% %%
% sz = [3 3 3];
% 
% xyz1 = 2
% n1   = 1;
% 
% xyz2 = 2
% n2   = 2;
% 
% [msk,d] = Dict2dTo3d.planeIntersections( sz, xyz1, n1, xyz2, n2 );
% line = Dict2dTo3d.patchIntersectionFromMask( patch1, xyz1, n1, msk );
% fprintf('line should be non-zero: ( %d %d %d )\n', line);
% 
% xyz2 = 1;
% [msk,d] = Dict2dTo3d.planeIntersections( sz, xyz1, n1, xyz2, n2 );
% line2 = Dict2dTo3d.patchIntersectionFromMask( patch1, xyz1, n1, msk );
% fprintf('line should be zero    : ( %d %d %d )\n', line2);

