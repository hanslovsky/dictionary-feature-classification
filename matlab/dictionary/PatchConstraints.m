classdef PatchConstraints < handle
    
    properties 
        f;  % downsampling factor
        sz2d;
        sz3d;
        
        pairLocRng;
        dimXyzList;
        xsectList;
        
        cmtx;
        locToConstraint;
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
            
            k = 1;
            for i = 1 : numPatchLocs
                thisdim = this.dimXyzList(i,1);
                thisxyz = this.dimXyzList(i,2);
                
                % TODO - dont recompute this msk every time
                msk = Dict2dTo3d.planeMaskF( this.sz3d, thisxyz, thisdim, this.f );
                
                for j = 1:patchNumElem
                    
                    this.cmtx( k, (msk==j) ) = 1;
                    this.locToConstraint( i, k ) = 1;
                    
                    k = k + 1;
                end
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
        
        function constraints = cmtxFromConstraintsList( this, locXyz )
           
            ndim  = nnz( size(locXyz) > 1);
            if( ndim == 1 )
                N = numel(locXyz);
            else
                N = size(locXyz,1);
            end
            
            cmtxIdxs = false( 1, size(locXyz,1));
            for d = 1:N 
                
                if( ndim == 1)
                    idx =  locXyz(d);
                elseif( size( locXyz, 2) == 2)
                    idx = this.locXyzDim2Idx(   locXyz(d,1), ...
                                                locXyz(d,2));
                    
                else
                    error('locXyz must have 1 or 2 columns ');
                end
                
                cmtxIdxs = cmtxIdxs | this.locToConstraint(idx,:);
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
        
    end
end