%% hrlr_samplingOpt_tests

import net.imglib2.algorithms.patch.SubPatch2dLocation;
import net.imglib2.algorithms.opt.astar.AStarMax;
import net.imglib2.algorithms.opt.astar.SortedTreeNode;
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

% patch5 = (1./3).*ones( 9, 9 );

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

% X = [   patch1(:), patch2(:), patch3(:), patch4(:), patch5(:), patch6(:), patch7(:), ...
%        patch8(:), patch9(:), patch10(:), patch11(:)]';

X = [   patch1(:), patch2(:), patch3(:), patch4(:), patch6(:), patch7(:), ...
        patch8(:), patch9(:), patch10(:), patch11(:)]';
    
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

%% test patch similarity / cost

% clear d23;
% sz = 9;
% f  = 3;
% d23 = Dict2dTo3d( X, sz, f );

% %%
% 
% dim1 = 3;
% xyz1 = 2;
% 
% dim2 = 1;
% xyz2 = 5;
% 
% p1 = patch1;
% p2 = patch10;
% 
% I1 = Dict2dTo3d.fill3dWith2d( [9 9 9], dim1, xyz1, f, p1 );
% I2 = Dict2dTo3d.fill3dWith2d( [9 9 9], dim2, xyz2, f, p2 );
% 
% [ sim1, x, cmtx, b, pm1, pm2, overlap ] = d23.patchSimilarity( p1, xyz1, dim1, p2, xyz2, dim2 );
% sim1
% 
% %% compare the above to the computation with allSims
% 
% allSims = d23.allSims;
% 
% idx1 = d23.locXyzDim2Idx( dim1, xyz1 )
% idx2 = d23.locXyzDim2Idx( dim2, xyz2 )
% 
% allCosts = allSims( 1, :, idx1, idx2 );
% allCosts
% 
%%
% numDifferent = 0;
% 
% for i = 1:size(X,1)
% %     p1 = reshape( X(i,:), [9 9]);
% %     for j = 1:size(X,2)
% %         p2 = reshape( X(j,:), [9 9 ]);
% %         
% %         for n = 1:size(d23.dimXyzList,1)
% %             dim1 = d23.dimXyzList(n,1);
% %             xyz1 = d23.dimXyzList(n,2);
% %             for m = 1:size(d23.dimXyzList,1)
% %                 dim2 = d23.dimXyzList(m,1);
% %                 xyz2 = d23.dimXyzList(m,2);
% %                 
% %                 [ sim1, x, cmtx, b, pm1, pm2, overlap ] = d23.patchSimilarity( p1, xyz1, dim1, p2, xyz2, dim2 );
% %                 [ sim2, x, cmtx, b, pm1, pm2, overlap ] = d23.patchSimilarity( p1, xyz2, dim2, p2, xyz1, dim1 );
% %                 
% %                 if ( sim1 ~= sim2 )
% %                     fprintf('different\n');
% %                     numDifferent = numDifferent + 1;
% %                 end
% %                 
% %             end
% %         end
% %     end
% end

%% test filling of 3d patches with 2d patches

% I1 = Dict2dTo3d.fill3dWith2d( [9 9 9], 1, 5, f, patch1 );
% I2 = Dict2dTo3d.fill3dWith2d( [9 9 9], 2, 5, f, patch1 );

%% test consistency of patches given constraints from other patches


% sz = 9;
% f  = 3;
% d23 = Dict2dTo3d( X, sz, f );
% 
% rootNode = SortedTreeNode(  ...
%     SubPatch2dLocation( 1, 5, 1, 0 ));
% 
% tmpNode = SortedTreeNode(  ...
%     SubPatch2dLocation( 2, 5, 6, -1 ), ...
%     rootNode );
% 
% [ xsectList ] = Dict2dTo3d.intersectingParents( tmpNode );
% size( xsectList )
% 
% iiThis = 5
% costs = d23.patchCosts( iiThis, xsectList );

%% build the 3d patch from 2d patches!!

clear d23;

sz = 9;
f  = 3;
d23 = Dict2dTo3d( X, sz, f );

%%

best = d23.build3dPatch();

root = d23.p2dFill3d.getRoot();

%% print out full parameters

printme = best
params = [  best.getData().dim, ...
            best.getData().xyz, ...
            best.getData().idx ];

i = 1;
while( ~printme.isRoot())
    printme = printme.getParent()
    params = [ params; ...
                printme.getData().dim, ...
                printme.getData().xyz, ...
                printme.getData().idx ];
                
end
params

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

