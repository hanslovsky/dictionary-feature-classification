classdef PatchConstraints < handle
    
    properties 
        f;  % downsampling factor
        sz2d;
        sz3d;
        
        pairLocRng;
        dimXyzList;
        xsectList;
        
        cmtx;     % the constraint matrix 
        cmtxInv;  % pseudo-inverse of the constraint matrix
        locToConstraint;
        
        subCmtxAndInvs;
        constraintVecXsectSubsets; % 
        constraintVecSubsets; % 
    end
    
    methods
        
        function this = PatchConstraints( sz, f )
        % Constructor
            this.f    = f;
            
            if( length( sz ) > 2 )
               error('sz must be a 2-vector or scalar');
            end
            if( length( sz ) == 2 )
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
            
            this.pairLocRng = (1) : this.f : this.sz3d(1);
            [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
            this.dimXyzList = [ dList(:), xyzList(:) ];
            
            this.buildIntersectionList();
        end
        
        function idx = locXyzDim2Idx( this, dim, xyz )
            idx = ( this.dimXyzList(:,1) == dim & ...
                    this.dimXyzList(:,2) == xyz );
        end
        
        function [dim,xyz] = locIdx2XyzDim( this, idx )
            dim = this.dimXyzList( idx, 1 );
            xyz = this.dimXyzList( idx, 2 );
        end
        
        function buildCmtx( this )
            patchNumElem = prod( this.sz2d );
            numVariables = prod( this.sz3d );
            
            numPatchLocs   = size( this.dimXyzList, 1 );
            numConstraints = numPatchLocs * patchNumElem;
            
            this.cmtx = zeros( numConstraints, numVariables );
            this.locToConstraint = false( numPatchLocs, numConstraints );
            this.constraintVecSubsets = false( numPatchLocs, numConstraints );
            
            k = 1;
            for i = 1 : numPatchLocs
                thisdim = this.dimXyzList(i,1);
                thisxyz = this.dimXyzList(i,2);
                
                rngStart = ( patchNumElem * (i-1) + 1);
                rng =  rngStart : rngStart + patchNumElem - 1;
                this.constraintVecSubsets( i , rng ) = true;
                
                % TODO - dont recompute this msk every time
                msk = Dict2dTo3d.planeMaskF( this.sz3d, thisxyz, thisdim, this.f );
                
                for j = 1:patchNumElem
                    
                    this.cmtx( k, (msk==j) ) = 1;
                    this.locToConstraint( i, k ) = 1;
                    
                    k = k + 1;
                end
            end
            
            this.cmtxInv = pinv( this.cmtx );
        end
        
        function compXsectInverses( this )
            N =  size(this.locToConstraint, 1);
            this.subCmtxAndInvs = cell( N, 2 );
            this.constraintVecXsectSubsets = cell( N ,2 );
            for i = 1:N
                inds = this.xsectList(i,:);
                
                [subcmtx, bxsectList] = ...
                        this.cmtxFromConstraintsList( inds );
                [cmtxThis, ~] = ...
                        this.cmtxFromConstraintsList( i );
                
                subcmtxWCurrent = [ cmtxThis; subcmtx ];
                cmtxInv = pinv( subcmtxWCurrent );
                this.subCmtxAndInvs{i,1} = subcmtxWCurrent;
                this.subCmtxAndInvs{i,2} = cmtxInv;
                
                % determine the 
                this.constraintVecXsectSubsets{i} = bxsectList;
            end
        end
        
        function constraints = cmtxFromConstraints( this, obj )
            if( isa( obj, 'net.imglib2.algorithms.opt.astar.SortedTreeNode'))
                constraints = this.cmtxFromConstraintsNode( obj );
            else
                constraints = this.cmtxFromConstraintsList( obj );
            end
        end
        
        function constraints = cmtxFromConstraintsNode( this, node )
            thisnode = node;
            depth = node.getDepth();
            
            cmtxIdxs = false( 1, size(this.locToConstraint,2));
            for d = 1:depth+1
                
                idx = this.locXyzDim2Idx( thisnode.getData().dim, ...
                                          thisnode.getData().xyz );
                
                cmtxIdxs = cmtxIdxs | this.locToConstraint(idx,:);
                thisnode = thisnode.getParent();
            end
              
            constraints = this.cmtx( cmtxIdxs, : );
        end
        
        function [constraints, binds] = cmtxFromConstraintsList( this, locXyz )
        % outputs a constraint matrix from a list of constraints

            if( islogical( locXyz ))
                locXyz = find( locXyz );
            end

            ndim  = nnz( size(locXyz) > 1);
            if( ndim == 1 )
                N = numel(locXyz);
            else
                N = size(locXyz,1);
            end
            
            cmtxIdxs = false( 1, size(locXyz,1));
            binds    = false( 1, prod( this.sz3d) );
            for d = 1:N 
                
                if( ndim <= 1)
                    idx =  locXyz(d);
                elseif( size( locXyz, 2) == 2)
                    idx = this.locXyzDim2Idx(   locXyz(d,1), ...
                                                locXyz(d,2));
                    
                else
                    error('locXyz must have 1 or 2 columns ');
                end
                
                cmtxIdxs = cmtxIdxs | this.locToConstraint(idx,:);
                binds = binds | this.constraintVecSubsets( idx,: );
            end
              
            constraints = this.cmtx( cmtxIdxs, : );
        end
        
        function xsectList = buildIntersectionList( this )
            % i^th row of output xsectList
            % is a logical vector such that the coordinates
            % dimXyzList( xsectList,:) intersect with dimXyzList( i, :)
            
            numPatchLocs   = size( this.dimXyzList, 1 );
            xsectList = false( numPatchLocs );
            
            for i = 1:numPatchLocs
                xsectList(i,:) = (this.dimXyzList(:,1) ~= this.dimXyzList(i,1));
            end
            
            this.xsectList = xsectList;
        end
           
        function b = constraintValue( this, patchMtx, obj )
            if( isa( obj, 'net.imglib2.algorithms.opt.astar.SortedTreeNode'))
                b = this.constraintValueNode( patchMtx, obj );
            else
                b = this.constraintValueList( patchMtx, obj );
            end
        end
        
        function b = constraintValueNode( this, patchMtx, node )
            b = [];
            error('not yet implemented');
        end
        
        function [paramsOut, rowPermutation] = reorderPatchParams( this, patchParams )
            % patchParams is an N x 3 matrix where the columns represent
            % [ dim, xyz, idx ].
            %
            % This function re-orders the rows of patchParams so that
            % the rows are in the same order as this.dimXyzList
            
            nr = size( this.dimXyzList, 1 );
            rowPermutation = zeros( nr, 1);
            
            for r = 1:nr
                dim = this.dimXyzList(r,1);
                xyz = this.dimXyzList(r,2);
                
                i = find( patchParams(:,1) == dim &  ...
                          patchParams(:,2) == xyz );
                rowPermutation(r) = i;
            end
            
            paramsOut = patchParams( rowPermutation, : );
            
        end
            
        function b = constraintValueList( this, patchMtx, idxList )
        % idxList must be in the same order
        
            if( islogical( idxList ))
                idxList = find( idxList );
            end
            
            ndim  = nnz( size(idxList) > 1);
            if( ndim == 1 )
                N = numel(idxList);
            else
                N = size(idxList,1);
            end
            
            patchNumElem = prod( this.sz2d );
            M = patchNumElem * N;
            
            b = zeros( M, 1 );
            brng = 1:patchNumElem;
            for i = 1:N
                if( ndim <= 1)
                    idx =  idxList( i );
                elseif( size( idxList, 2) == 2 )
                    idx = this.locXyzDim2Idx(   idxList(i,1), ...
                                                idxList(i,2));
                end
                b(brng) = patchMtx( idx, : );
                brng = brng + patchNumElem;
            end
        end
        
        function bnew = updateConstraints( this, patchMtx, b, j, idx )
        % bnew = updateConstraints( this, patchMtx, b, j, idx )
        %   patchMtx - N x M matrix where N is number of dictionary elements
        %   b        - current b vector
        %   j        - index of patch location (into dimXyzList)
        %   idx      - index into patchMtx of replacing patch
            bnew = b;
            patchNumElem = prod( this.sz2d );
            start = patchNumElem * (j - 1) + 1;
            rng = start : start + patchNumElem - 1;
            
            bnew( rng ) = patchMtx( idx, : );
        end
        
    end
end
