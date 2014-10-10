classdef Dict2dTo3d < handle
    % Dict2dTo3d 
    %
    % This class contains logic that builds a high-res 3d dictionary from 
    % 2d patches at a lower resolution.
    %
    % John Bogovic
    % HHMI
    % September 2014
    
    properties ( SetAccess = private )
        
        D2d;      % The dictionary elements
        numDict;  % Number of 2d dictionary elements
        sums2d;   % the dictionary element sums
        
        sz2d;     % the size of 2d patches
        sz3d;     % the size of 3d patches
        
        f;        % the downsampling factor ( resolution in z of patches )
        
        summer2Dxy; % downsampling by sum function
        
        % Constraint matrices
        constraintMtxX;
        constraintMtxY;
        constraintMtxZ;
        
        % 
        p2dFill3d;
        pairLocRng;
        dimXyzList;
        allSims;
        
        
        % optimization options
        Nbest = 3;
        maxItersPerPatch = 5000;
        
        ndims = 3; 
    end
    
    methods
        
        % Constructor
        % D2d - 
        % sz  - size of 2d patches
        % f   - downsampling factor
        function this = Dict2dTo3d( D2d, sz, f )
            
            
%             import net.imglib2.algorithms.opt.astar.AStarMax;
            import net.imglib2.algorithms.patch.Patch2dFill3d;
            import java.util.*;

            this.D2d     = D2d;
            this.numDict = size( this.D2d, 1 );
           
            this.f    = f;
            
            if( length(sz ) > 2 )
               error('sz must be a 2-vector or scalar');
            end
            if( length(sz ) == 2 )
                if( sz(1) ~= sz(2) )
                    error('sz must be a 2-vector');
                end
                this.sz2d = sz;
                this.sz3d = [ sz sz(1) ];
            end
            if( isscalar( sz ))
                this.sz2d = [ sz sz ];
                this.sz3d = [ sz sz sz ];
            end
                        
%             this.sums2d = Dict2dTo3d.allSums( this.D2d, this.sz2d, this.f);
%             this.summer2Dxy = Tid.sum2dxy();

            half = (this.f - 1)./2;
            this.pairLocRng = (1+half) : this.f : this.sz3d(1);
            
            this.allSims = this.allSimilarities();

            [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
            this.dimXyzList = [ dList(:), xyzList(:) ];
            
        end
        
        function idx = locXyzDim2Idx( this, dim, xyz )
            idx = ( this.dimXyzList(:,1) == dim & ...
                    this.dimXyzList(:,2) == xyz );
        end
        
        function [dim,xyz] = locIdx2XyzDim( this, idx )
            dim = this.dimXyzList( idx, 1 );
            xyz = this.dimXyzList( idx, 2 );
        end
        
        function patch = sample2dTo3d( this )
            import net.imglib2.algorithms.patch.SubPatch2dLocation;
            
            szSm = this.sz2d ./ this.f;
            
            % stores the order in which patches are added to particular bins 
            % (hard-coded '3' corresponds to 3 spatial-dimensions)
            param = zeros( [ 3 max(szSm) ] ); 
            
            % stores which observed 2d patch goes into which 3d position / orientation "bin"
            patchIdxs = zeros( [ 3 max(szSm) ] );
            
            % patch Loc and Idx
            % [ dim xyz idx ]
            patchLocIdx = zeros( 3.*max(szSm), 3 );
            
            % how many bins are there
            mskN = numel( param );
            
            % randomize the order in which the bins will be filled
            xyzn = randperm( mskN );
            
            for ii = 1:mskN
                fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
                
                dim = this.dimXyzList(ii,1);
                xyz = this.dimXyzList(ii,2);
                dxyzi1 = xyzn(ii);
                
                if( ii > 1 )
                    
%                     xsectList = Dict2dTo3d.findIntersections( param, [dim xyz]);
%                     xsectList

                    [ xsectCoords, ~ ] = Dict2dTo3d.findIntersectionsList( patchLocIdx, dim,  xyz);
                    % xsectCoords
                    
                    % get the next best candidate for patch 
                    
                    
                    for nn = 1:size( xsectCoords, 1 )
                        
                        fixedDim = xsectCoords(nn,1);
                        fixedXyz = xsectCoords(nn,2);
                        fixedIdx = xsectCoords(nn,3);
                        
                        dxyzi2 = ( this.dimXyzList(:,1) == fixedDim & ...
                                   this.dimXyzList(:,2) == fixedXyz);
                                
%                         dxyzi2 = find( this.dimXyzList(:,1) == fixedDim & ...
%                                        this.dimXyzList(:,2) == fixedXyz);
%                         tsim = this.allSims(  fixedIdx, :, dxyzi1, dxyzi2 );

                        % get similarities for all patches and the given orientations
                        tsim = this.allSims(  fixedIdx, :, dxyzi1, dxyzi2 )

                        [tsimSort, sortIdx ] = sort( tsim );
                      
                        for mm = 1:this.Nbest
                            aa=1;
                        end
                        
                    end
                    
                else
                    rootSpl = SubPatch2dLocation( dim, xyz, 0 );
                    this.astar.setRoot( rootSpl );
                    patchIdxs(ii) = randi( this.numDict, 1 ); 
                    patchLocIdx(ii,:) = [ dim, xyz, patchIdxs(ii) ];
                end
                
                param( xyzn(ii) ) = ii;
                patchLocIdx(ii,:) = [ dim, xyz, randi( this.numDict, 1 ) ];
                
            end
            
            patch = [];
        end
        
        function [patchParams,iteration] = build3dPatch( this )
            
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            N = size( this.dimXyzList, 1 );
            
            % randomize the order in which the patch bins will be filled
            xyzn = randperm( N );
            
            done = false;
            iteration = 1;
            while( ~done )
%                 iteration 

                if ( iteration > this.maxItersPerPatch )
                    
                    patchParams = this.currentBestPatchParams();
                    return;
                    
                elseif( iteration == 1 )
                   
                    ii = xyzn( 1 );
                    dimIni = this.dimXyzList( ii, 1 );
                    xyzIni = this.dimXyzList( ii, 2 );
                    
                    % dxyzi1 = xyzn(ii);
                    
                    initialPatchIdx = randi( this.numDict );
                    rootSpl = SubPatch2dLocation( dimIni, xyzIni, initialPatchIdx, 0 );
                    this.p2dFill3d = Patch2dFill3d( rootSpl );
                    
                else
                    
                    prevNode = this.p2dFill3d.next();
                    depth = prevNode.getDepth();
                    if( depth == N-1 )
                        
                        % the next node is best and describes the 
                        % parameters for the locally optimal patch
                        patchParams = prevNode;
                        
                        break;
                    end
                    
                    iiPrev = xyzn( depth + 1 );
                    iiThis = xyzn( depth + 2 );
                    
                    % dimPrev = this.dimXyzList( iiPrev, 1 );
                    % xyzPrev = this.dimXyzList( iiPrev, 2 );
                    
                    dimThis = this.dimXyzList( iiThis, 1 );
                    xyzThis = this.dimXyzList( iiThis, 2 );
                    
                    tmpNode = SortedTreeNode(  ...
                                SubPatch2dLocation( dimThis, xyzThis, -1, -1 ), ...
                                prevNode );
                    [ xsectList ] = Dict2dTo3d.intersectingParents( tmpNode );
%                     size( xsectList )
                            
                    % TODO - CHECK IF PARENT LIST IS EMPTY - done!
                    %        INCLUDE CONSTRAINTS OF nextNode in list - done
                    % % This is old and may not be correct
%                     [ xsectParentList ] = Dict2dTo3d.intersectingParents( prevNode );
%                     size( xsectParentList )
%                     
%                     xsectList = [ prevNode.getData().dim prevNode.getData().xyz prevNode.getData().idx; ...
%                                   xsectParentList ];
%                     size( xsectList )
                              

                    if( isempty( xsectList ))
%                         fprintf('no intersections...picking randomly\n');
                        % pick 'Nbest' patches randomly and give them equal cost 
                        randomPatchIndexes = randi( this.numDict, 1, this.Nbest);
                        costs = ones( this.numDict, 1 );
                        % TODO - what cost should be given?
                        costs( randomPatchIndexes ) = 0;
                    else
%                         fprintf('computing costs\n');
                        costs = this.patchCosts( iiThis, xsectList );
                    end

                    candList = java.util.ArrayList( this.Nbest );
                    [ sortedCosts, sortedCostIdxs ] = sort( costs );
                    for nn = 1:length( sortedCosts )
                       
                       val = sortedCosts(nn);
                       spl = SubPatch2dLocation( dimThis, xyzThis, sortedCostIdxs(nn), val );
                       candList.add( spl );
                    end
                    
                    % % for debug purposes
%                     candList = java.util.ArrayList();
%                     for nn = 1:this.Nbest
%                        % pick random value for now
%                        val = rand;
%                        spl = SubPatch2dLocation( dim, xyz, nextPatchIdx(nn), val );
%                        candList.add( spl );
%                     end
                    
                    this.p2dFill3d.addFromNode( prevNode, candList );
                    
                    fprintf('iteration %d\n', iteration);
                    fprintf('set has %d elements\n',this.p2dFill3d.getSet().size());
                    % this.p2dFill3d
                    
%                     if( depth == N - 1 )
%                         done = true;
%                     end
                end
                iteration = iteration + 1;
                %pause;
            end
        end
        
        function currentBest = currentBestPatchParams( this )
            N = size( this.dimXyzList, 1 );
            it = this.p2dFill3d.getSet().iterator();
            while( it.hasNext())
               thisone = it.next(); 
               if( thisone.getDepth() == (N-1))
                  currentBest = thisout;
                  return;
               end
            end
            currentBest = [];
        end
        
        function constraints = getSumConstraints( this, patchIdxs )
            
            szdown = this.sz3d ./ this.f;
            constraints = nan( szdown );
            
            [d,xyz] = ind2sub( size( patchIdxs ), find( patchIdxs > 0 ));
            for i = 1:length(xyz)
                patch = reshape( this.D2d( patchIdxs(d(i),xyz(i)), : ), this.sz2d );
                constraints( Dict2dTo3d.planeMask( szdown, xyz(i), d(i) )) = this.summer2Dxy( patch, this.f );
            end
        end
        
        function c = collectConstraints( this, constraints, xyz, d )
            szSm3 = this.sz3d ./ this.f;
            szSm2 = this.sz2d ./ this.f;
            planeMask = Dict2dTo3d.planeMask( szSm3, xyz, d );
            c = reshape( constraints( planeMask ), szSm2 );
        end
        
        function buildConstraintMatrix( this )
           
            [x,y,z] = ndgrid( 1:this.f, 1:this.f, 1:this.f );
            N = numel(x);
            
            this.constraintMtxX = zeros( this.f, N );
            this.constraintMtxY = zeros( this.f, N );
            this.constraintMtxZ = zeros( this.f, N );
            
            for i = 1:this.f
                % rows
                this.constraintMtxX(i,:) = reshape((x==i),1,[]);
                
                % colums
                this.constraintMtxY(i,:) = reshape((y==i),1,[]);
                
                % slices
                this.constraintMtxZ(i,:) = reshape((z==i),1,[]);
                
            end
        end

        % returns 
        function [ costs ] = patchCosts( this, dxyzi1, intersectionList )
            
            N = size(intersectionList,1);
            constraintCosts = zeros( N, this.numDict );
            
            for nn = 1:N
                
                fixedDim = intersectionList(nn,1);
                fixedXyz = intersectionList(nn,2);
                fixedIdx = intersectionList(nn,3);
                
                dxyzi2 = ( this.dimXyzList(:,1) == fixedDim & ...
                           this.dimXyzList(:,2) == fixedXyz);
                               
                constraintCosts( nn,:) = this.allSims(  fixedIdx, :, dxyzi1, dxyzi2 );
            end
            
%             fprintf('constraint costs\n');
%             size( constraintCosts )
            
            % take max over constraint costs for final cost
            costs = max( constraintCosts, [], 1 );
            
%             fprintf('new patch costs\n');
%             size(costs)
        end
        
        function [ cost ] = patchCostSlow( this, p1, xyz1, n1, intersectionList )
            N = size(intersectionList,1);
            similarityList = zeros( N, 1);
            for i = 1:N
                p2 = reshape( this.X( intersectionList(i,3), :), this.sz2d );    
                similarityList(i) = this.patchSimilarity( p1, xyz1, n1, ...
                                        p2, intersectionList(i,1), intersectionList(i,2));
                                    
            end
            cost = max( similarityList );
        end
        
        % make this methods abstract when more possibilities are available
        % this version will use unconstrained linear least squares 
        function [ sim, x, cmtx, b, pm1, pm2, overlap ] = patchSimilarity( this, p1, xyz1, n1, p2, xyz2, n2 )
            
            % get the value of the sum-contraint each patch describes over
            % the overlapping area
            sz3 = this.sz3d;

            pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, n1, this.f );
            pm2 = Dict2dTo3d.planeMaskF( sz3, xyz2, n2, this.f );
            overlap = (pm1 > 0) & (pm2 > 0);
           
            [cmtx1, b1] = Dict2dTo3d.contraintsFromMask( p1, overlap, pm1 );
            [cmtx2, b2] = Dict2dTo3d.contraintsFromMask( p2, overlap, pm2 );

            cmtx = [ cmtx1; cmtx2 ] ;
            b    = [ b1; b2 ];
            x = pinv(cmtx) * b;
            
            sim = sum((cmtx*x - b).^2);
        end
        
        % returns a 4d array describing the similarities between all pairs
        % of patches in all positions / orientations
        % 
        % dimensions are : ( patch1Index patch2Index patch1XyzDim patch2XyzDim )
        % the this.dimXyzList property stores the patch orientations.
        %
        % As an example:
        % allSims( i, j, k, l ) stores the similarities for patches i and j
        %   for orientations 
        %   this.dimXyzList( k ) and
        %   this.dimXyzList( l ) 
        %
        function allSims = allSimilarities( this )

            [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
            N = numel(dList);
            
            this.dimXyzList = [ dList(:), xyzList(:) ];
            
            allSims = zeros( this.numDict, this.numDict, N, N );
            for i = 1:this.numDict
                p1 = reshape( this.D2d(i,:), this.sz2d );
                for j = i:this.numDict
                    p2 = reshape( this.D2d(j,:), this.sz2d );
                    for k = 1:N
                        xyz1 = xyzList(k);
                        n1   = dList(k);
                        for l = 1:N
                            xyz2 = xyzList(l);
                            n2   = dList(l);
                            sim = this.patchSimilarity( p1, xyz1, n1, p2, xyz2, n2 );
                            
                            allSims( i, j, k, l ) = sim;
                            allSims( j, i, k, l ) = sim;
                        end
                    end
                end
            end
            
        end
        
    end
    
    methods( Static )
        
        % isRow - if true, sample ith row, 
        %         else sample ith column
        function slc = slice( patch, i, isRow )
            if( isRow )
                slc = patch(i,:);
            else
                slc = patch(:,i);
            end
        end
        
        % isRow - if true, sample ith row, 
        %         else sample ith column
        function slc = slicePadded( patch, i, pad, isRow )
            if( isRow )
                slc = patch(i-pad:i+pad,:);
            else
                slc = patch(:,i-pad:i+pad);
            end
        end
        
        % downsample vector 
        function ds = downsampleVec( vec, f )
            ds = vecToRowCol( vec, 'row' );
            ds = mean(reshape( ds, f, [] ));
        end
        
        function similarity = BADpatchConsistency2d( p1, i1, d1, p2, i2, d2, f )
            
            slc1 = Dict2dTo3d.slice( p1, i1, d1 );
            slc2 = Dict2dTo3d.slice( p2, i2, d2 );
            
            slc1ds = Dict2dTo3d.downsampleVec( slc1, f );
            slc2ds = Dict2dTo3d.downsampleVec( slc2, f );
            
            similarity = norm( slc1ds - slc2ds );
        end
        
        function similarity = patchConsistency2d( p1, i1, d1, p2, i2, d2, f )
            
            half = (f-1)./2;
            slc1 = Dict2dTo3d.slicePadded( p1, i1, half, d1 );
            slc2 = Dict2dTo3d.slicePadded( p2, i2, half, d2 );
            
            dsFun = Tid.sum2dxy();
            
            slc1ds = sum(dsFun( slc1, f ));
            slc2ds = sum(dsFun( slc2, f ));
            
            similarity = norm( slc1ds - slc2ds );
        end
        
        function similarity = patchConsistencyConstraints( constraints, X, sz, f )
           
            N = size( X, 1 );
            similarity = zeros( N, 1 );
            dsFun = Tid.sum2dxy();
            
            for i = 1:N
                
                p = reshape( X(i,:), sz );
                pd = dsFun( p, f );
                notNanInds = ~isnan(constraints);
%                 constraints(notNanInds)
%                 pd(notNanInds)
%                 constraints(notNanInds) - pd(notNanInds)
%                 fprintf('\n\n')
                similarity(i) = norm(constraints(notNanInds) - pd(notNanInds));
                
            end
        end
        
        function psv = patchSumVecs( p, i, isRow, downsampleFactor )
            half = (downsampleFactor-1)./2;
            slc = Dict2dTo3d.slicePadded( p, i, half, isRow );
            dsFun = Tid.sum2dxy();
            psv = dsFun( slc, downsampleFactor );
        end
        
        function sumMtx = allSums( X, sz, f )
            [N,M] = size( X );
            if( prod(sz) ~= M )
                error('size inconsistent with patch matrix X')
            end
            
            half = (f-1)./2;
            numValid = sz - ( 2 * half );
            L = sum( numValid );
            
            sumSz = max(sz./f);
            
            sumMtx = zeros( N*L, sumSz );
            
            k = 1;
            for n = 1:N
                for d = 0:1
                    for ri = 2 : sz(d+1)-1;
                        sumMtx( k, : ) = Dict2dTo3d.patchSumVecs( reshape( X(n,:), sz ), ri, d, f );
                        k = k + 1;
                    end
                end
            end
        end
        
        % 
        function intersection = findIntersections( msk, v )
            j = setdiff( 1:3, v(1));
            m = msk( j, : );
            intersection = m( m > 0 );
        end
        
        function [coords, logicalIdx] = findIntersectionsList( xyzDimList, dim, xyz )
            logicalIdx = (xyzDimList(:,2) > 0 ) & (xyzDimList(:,1) ~= dim);
            coords = xyzDimList( logicalIdx, : );
        end
                
        % xyz{1,2} are 3-vectors giving a point planes 1 and 2
        % n{1,2} are 3-vectors giving the normal vectors 
        function [msk,d] = planeIntersections( sz, xyz1, n1, xyz2, n2 )
        %function [msk,d,ix,iy,iz] = planeIntersections( sz, xyz1, n1, xyz2, n2 )
            
            msk1 = Dict2dTo3d.planeMask( sz, xyz1, n1);
            msk2 = Dict2dTo3d.planeMask( sz, xyz2, n2);
            msk = msk1 & msk2;

            if( nargout == 1)
                return;
            end
            
            if( n1 == n2 )
                d = 0;
            else
                d = setdiff( 1:3, [n1 n2] );
            end
            
            %if( nargout > 2 )
            %    [x,y,z] = meshgrid(1:sz(1), 1:sz(2), 1:sz(3));
            %    ix = x(msk);
            %    iy = y(msk);
            %    iz = z(msk);
            %end
            
        end

        function line = patchIntersectionFromMask( patch, xyz, n, msk )
            sz = size( msk );
            pmsk = Dict2dTo3d.planeMask( sz, xyz, n );
            line = patch( msk( pmsk )); 
        end
        
        function I = fill3dWith2d( sz3d, dim, xyz, f, patch)
            I = zeros( sz3d );
            msk = Dict2dTo3d.planeMaskF( sz3d, xyz, dim, f);
            I( msk > 0 ) = patch( msk(msk>0) );
        end
        
        % n is {1,2,3}
        function msk = planeMask( sz, xyz, n )
            msk = false( sz );
            
            if( isscalar(xyz) )
                val = xyz;
            else
                switch n
                    case 1
                        val = xyz(1);
                    case 2
                        val = xyz(2);
                    case 3
                        val = xyz(3);
                    otherwise
                        error('invalid normal direction');
                end
            end
            
            switch n
               case 1
                   msk( val, :, : ) = true;
               case 2
                   msk( :, val, : ) = true;
               case 3
                   msk( :, :, val ) = true;
               otherwise 
                   error('invalid normal direction');
            end
        end
        
        % n is {1,2,3}
        function msk = planeMaskF( sz, xyz, n, f, centered )
            
            msk = zeros( sz );
            
            if( ~exist('centered','var') || isempty( centered ))
               centered = false; 
            end
            
            if( isscalar(xyz) )
                val = xyz;
            else
                switch n
                    case 1
                        val = xyz(1);
                    case 2
                        val = xyz(2);
                    case 3
                        val = xyz(3);
                    otherwise
                        error('invalid normal direction');
                end
            end
            
            if( centered )
                half = (f-1)./2;
                rng = val-half : val+half;
            else
                rng = val : val + f - 1;
            end
            N = prod(sz(1:2));
            
            v = repmat( reshape(1:N, sz(1:2)), [1 1 f]);
            
            switch n
                case 1
                    msk( rng, :, : ) = permute(v, [3 1 2]);
                case 2
                    msk( :, rng, : ) = permute(v, [1 3 2]);
                case 3
                    msk( :, :, rng ) = v;
                otherwise
                    error('invalid normal direction');
            end
        end
        
        % Deprecated
        function intersection = findIntersectionsBig( msk, v )
           n = v(4);
           switch n
               case 1
                   m = msk( v(1), :, :, [2 3] );
               case 2
                   m = msk( :, v(2), :, [1 3] );
               case 3
                   m = msk( :, :, v(3), [1 2] );
               otherwise 
                   error('invalid normal vector');
           end
           intersection = m( m > 0 );
        end
        
        % X must be (N x M)
        % where N is the number of patches and
        % M is prod( sz )
        % and sz is a vector describing the patch size
        function simMtx = allSimilaritiesSums( X, sz ) 
            [N,M] = size( X );
            if( prod(sz) ~= M )
                error('size inconsistent with patch matrix X')
            end
            half = (f-1)./2;
            numValid = sz - ( 2 * half );
            L = sum( numValid );
            
            simMtx = zeros( N , L );
        end
        
        function [cmtx, b] = contraintsFromMask( patch, constraintPtMask, pm )
            
            pm = pm( constraintPtMask );
             
            % remove zero 
            % since it indicates points not in the patch
            patchIdxsInMask = setdiff( unique( pm ), 0 );
            
            N = length( patchIdxsInMask );  % the number of constraints
            M = nnz( constraintPtMask ); % the number of elements in the HR patch
            cmtx = zeros( N , M );
            
            for i = 1:N
                cmtx( i,  (pm == patchIdxsInMask(i)) ) = 1;
            end
            
            b = patch(  pm( patchIdxsInMask ));
        end
        
        function [ xsectParentList ] = intersectingParents( dimXyzNode )
            
            depth = dimXyzNode.getDepth();
            xsectParentList = zeros( depth, 3 );
            parent = dimXyzNode.getParent();
            n = 0;
            for i = 1:depth
                
                if( dimXyzNode.getData().dim ~= parent.getData().dim )
                    n = n + 1;
                    xsectParentList( n, : ) = [ parent.getData().dim, ...
                                                parent.getData().xyz,...
                                                parent.getData().idx ];
                end
                
                parent = parent.getParent();
            end
            
            xsectParentList = xsectParentList(1:n,:);
            
        end
        
    end % static methods
    
end % class
