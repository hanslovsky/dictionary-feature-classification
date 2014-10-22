classdef Dict2dTo3dConstr < Dict2dTo3d
    % Dict2dTo3d 
    %
    % This class contains logic that builds a high-res 3d dictionary from 
    % 2d patches at a lower resolution.
    %
    % John Bogovic
    % HHMI
    % September 2014
    

    properties ( SetAccess = private )
       qpOps; % quadratic program options
    end
    
    methods
        
        % Constructor
        % D2d - 
        % sz  - size of 2d patches
        % f   - downsampling factor
        function this = Dict2dTo3dConstr( D2d, sz, f, opts )
            
            this = this@Dict2dTo3d( D2d, sz, f );
            
            if( ~exist( 'opts', 'var')  || isempty(opts))
                this.qpOps = optimoptions('quadprog');
                this.qpOps.Display = 'off';
                this.qpOps.Algorithm = 'interior-point-convex';
            else
                this.qpOps = opts;
            end
            this.qpOps
            
        end

        function [patchParams,iteration] = build3dPatch( this, iniPatchFill, num2exclude )
            
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            if( ~exist('num2exclude','var') || isempty(num2exclude))
                num2exclude = 0;
            end
            
            N = size( this.dimXyzList, 1 );
            
            % randomize the order in which the patch bins will be filled
            xyzn = randperm( N );
            
            haveIni = 0;
            if( exist( 'iniPatchFill', 'var' ) && ~isempty(iniPatchFill))
                
                fprintf('initializing with the input\n');
                haveIni = 1;
                
                this.p2dFill3d = Patch2dFill3d( iniPatchFill );
            
                iniParams = Dict2dTo3d.patchParamNodeToArray( iniPatchFill );
                
                exclude = find( ismember(this.dimXyzList(:,1), iniParams(:,1)) & ...
                            ismember(this.dimXyzList(:,2), iniParams(:,2)));
                
                xyzn = [exclude', setdiff( xyzn, exclude, 'stable')];
                
            end
            
            done = false;
            iteration = 1;
            while( ~done )

                if ( iteration > this.maxItersPerPatch )
                    
                    patchParams = this.currentBestPatchParams();
                    return;
                    
                elseif( iteration == 1 && ~haveIni )
                   
                    fprintf('initializing with random 2d patch\n');
                    ii = xyzn( 1 );
                    dimIni = this.dimXyzList( ii, 1 );
                    xyzIni = this.dimXyzList( ii, 2 );
                    
                    % dxyzi1 = xyzn(ii);
                    
                    initialPatchIdx = randi( this.numDict );
                    rootSpl = SubPatch2dLocation( dimIni, xyzIni, initialPatchIdx, 0 );
                    rootNode = SortedTreeNode( rootSpl );
                    this.p2dFill3d = Patch2dFill3d( rootNode );
                    
                else
                    
                    prevNode = this.p2dFill3d.next();
                    depth = prevNode.getDepth();
                    if( depth == N-1 )
                        
                        % the next node is best and describes the 
                        % parameters for the locally optimal patch
                        patchParams = prevNode;
                        
                        break;
                    end
                    
                    % iiPrev = xyzn( depth + 1 );
                    iiThis = xyzn( depth + 2 );
                    
                    dimThis = this.dimXyzList( iiThis, 1 );
                    xyzThis = this.dimXyzList( iiThis, 2 );
                    
                    tmpNode = SortedTreeNode(  ...
                                SubPatch2dLocation( dimThis, xyzThis, -1, -1 ), ...
                                prevNode );
                    [ xsectList ] = Dict2dTo3d.intersectingParents( tmpNode );
                   
                    if( isempty( xsectList ))
                        
%                         fprintf('no intersections...picking randomly\n');
                        % pick 'Nbest' patches randomly and give them equal cost 
                        randomPatchIndexes = randi( this.numDict, 1, this.Nbest);
                        costs = ones( this.numDict, 1 );
                        
                        % TODO - what cost should be given?
                        costs( randomPatchIndexes ) = 0;
                    else
%                         fprintf('computing costs\n');
                        if( isempty( this.allSims ))
                            costs = this.patchCostsAllSlow( xyzThis, dimThis, xsectList );
                        else
                            costs = this.patchCosts( iiThis, xsectList );
                        end
                        
                    end

                    costs = costs( 1: length(costs)-num2exclude);
                    
                    [ sortedCosts, sortedCostIdxs ] = sort( costs );
                    
                    
                    candList = java.util.ArrayList( this.Nbest );
                    for nn = 1:this.Nbest
                       val = sortedCosts(nn);
                       spl = SubPatch2dLocation( dimThis, xyzThis, sortedCostIdxs(nn), val );
                       candList.add( spl );
                    end
                    
                    this.p2dFill3d.addFromNode( prevNode, candList );
                    
                end
                iteration = iteration + 1;
            end
        end
        
            
        % make this methods abstract when more possibilities are available
        % this version will use unconstrained linear least squares 
        function [ sim, x, Aeq, b, pm1, pm2, overlap ] = patchSimilarity( this, p1, xyz1, n1, p2, xyz2, n2 )
            
            % get the value of the sum-contraint each patch describes over
            % the overlapping area
            sz3 = this.sz3d;

            pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, n1, this.f );
            pm2 = Dict2dTo3d.planeMaskF( sz3, xyz2, n2, this.f );
            overlap = (pm1 > 0) & (pm2 > 0);
            
            numConstraintsPerPatch = length(unique( pm1(overlap) ));
            N = nnz( overlap );
            
            H = [eye(N) -eye(N); -eye(N) eye(N)];
            f = zeros( 2*N, 1 );
            
            [cmtx1, b1] = Dict2dTo3d.contraintsFromMask( p1, overlap, pm1 );
            [cmtx2, b2] = Dict2dTo3d.contraintsFromMask( p2, overlap, pm2 );

            Aeq = [  cmtx1, zeros( numConstraintsPerPatch, N ); ...
                     zeros( numConstraintsPerPatch, N ), cmtx2 ];
            b    = [ b1; b2 ];
%             lb = repmat( min( b ), 1, length(b)); % lower bound
%             ub = repmat( max( b ), 1, length(b)); % upper bound
            lb = [];
            ub = [];
            [ x ] = quadprog(H, f, [], [], Aeq, b, lb, ub, [], this.qpOps );
            
            sim = norm( x(1:N) - x(N+1:end) );
            
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
        function allSims = allSimilaritiesFast( this )
            
            numP = size( this.D2d, 1 );
            
            sz3 = this.sz3d;
            [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
            numOrient = numel(dList);
            
            % this.dimXyzList = [ dList(:), xyzList(:) ];
            
            allSims = 9999.*ones( numP, numP, numOrient, numOrient );
            for k = 1:numOrient
                xyz1 = xyzList(k);
                n1   = dList(k);
                
                for l = 1:numOrient
                    
                    % don't bother computing similarities
                    % for identical orientations - they'll
                    % never be used
                    if( l == k )
                        continue;
                    end
                    
                    xyz2 = xyzList(l);
                    n2   = dList(l);
                    
                    % skip when the constraint planes are 'parallel'
                    % because this will never be used.
                    if( n1 == n2 )
                       continue; 
                    end
                    
                    pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, n1, this.f );
                    pm2 = Dict2dTo3d.planeMaskF( sz3, xyz2, n2, this.f );
                    overlap = (pm1 > 0) & (pm2 > 0);
                    
                    numConstraintsPerPatch = length(unique( pm1(overlap) ));
                    N = nnz( overlap );
                    
                    [ cmtx1, idxs1 ] = Dict2dTo3d.contraintsMtx( overlap, pm1 );
                    [ cmtx2, idxs2 ] = Dict2dTo3d.contraintsMtx( overlap, pm2 );
                    
                    H = [eye(N) -eye(N); -eye(N) eye(N)];
                    f = zeros( 2*N, 1 );
                    Aeq = [  cmtx1, zeros( numConstraintsPerPatch, N ); ...
                             zeros( numConstraintsPerPatch, N ), cmtx2 ];
                    
                    for i = 1:numP
                        b1 = this.D2d( i, idxs1 );
                        
                        for j = i:numP
                            
                            b2 = this.D2d( j, idxs2 );
                            
                            b    = [ b1'; b2' ];
%                             lb = repmat( min( b ), length(f), 1 );
%                             ub = repmat( max( b ), length(f),1 );
                            lb = [];
                            ub = [];
                            [ x ] = quadprog(H, f, [], [], Aeq, b, lb, ub, [], this.qpOps );
                            
                            if( isempty(x))
                                sim = 8888;
                            else
                                sim = norm( x(1:N) - x(N+1:end) );
                            end
                            

                            allSims( i, j, k, l ) = sim;
                            allSims( j, i, k, l ) = sim;

                        end
                    end
                end
%                 endtime = toc;
%                 fprintf('time: %fs\n', endtime );
            end
            this.allSims = allSims;
        end
        
        
        function allSims = specSimilaritiesFast( this, iprange, jprange )
            
            asym = 0;
            if( ~exist( 'iprange', 'var') || isempty( iprange ) )
               fprintf('default iprange\n');
               iprange =  1 : size( this.D2d, 1 );
               asym = 1;
            end
            
            if( ~exist( 'jprange', 'var') || isempty( jprange ) )
               fprintf('default jprange\n');
               jprange =  1 : size( this.D2d, 1 );
               asym = 1;
            end
            
            numIp = length( iprange );
            numJp = length( jprange );
            
            sz3 = this.sz3d;
            [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
            numOrient = numel(dList);
            
            % this.dimXyzList = [ dList(:), xyzList(:) ];
            
            allSims = 9999.*ones( numIp, numJp, numOrient, numOrient );
            for k = 1:numOrient
                xyz1 = xyzList(k);
                n1   = dList(k);
                
                for l = 1:numOrient
                    
                    % don't bother computing similarities
                    % for identical orientations - they'll
                    % never be used
                    if( l == k )
                        continue;
                    end
                    
                    xyz2 = xyzList(l);
                    n2   = dList(l);
                    
                    % skip when the constraint planes are 'parallel'
                    % because this will never be used.
                    if( n1 == n2 )
                       continue; 
                    end
                    
                    pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, n1, this.f );
                    pm2 = Dict2dTo3d.planeMaskF( sz3, xyz2, n2, this.f );
                    overlap = (pm1 > 0) & (pm2 > 0);
                    
                    N = nnz( overlap );
                    numConstraintsPerPatch = length(unique( pm1(overlap) ));
                    
                    H = [eye(N) -eye(N); -eye(N) eye(N)];
                    f = zeros( 2*N, 1 );
                    
                    [ cmtx1, idxs1 ] = Dict2dTo3d.contraintsMtx( overlap, pm1 );
                    [ cmtx2, idxs2 ] = Dict2dTo3d.contraintsMtx( overlap, pm2 );
                    
                    Aeq = [  cmtx1, zeros( numConstraintsPerPatch, N ); ...
                             zeros( numConstraintsPerPatch, N ), cmtx2 ];
                    
                    for i = 1:numIp
                        b1 = this.D2d( iprange(i), idxs1 );
                        
                        for j = 1:numJp
                            
                            b2 = this.D2d( jprange(j), idxs2 );
                            
                            b    = [ b1'; b2' ];
%                             lb = repmat( min( b ), 1, length(b));
%                             ub = repmat( max( b ), 1, length(b));
                            lb = [];
                            ub = [];
                            [ x ] = quadprog(H, f, [], [], Aeq, b, lb, ub, [], this.qpOps );
                            sim = norm( x(1:N) - x(N+1:end) );
                           
                            allSims( i, j, k, l ) = sim;
                            
                        end
                    end
                end
%                 endtime = toc;
%                 fprintf('time: %fs\n', endtime );
            end
            
        end
        
        
    end
    
    
end % class
