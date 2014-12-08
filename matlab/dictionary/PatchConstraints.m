classdef PatchConstraints < handle
    
    properties ( SetAccess = protected )
        f;  % downsampling factor
        sz2d;
        sz3d;
        
        numConstraints;
        numLocs;
        
        pairLocRng;
        dimXyzList;
        xsectList;
        
        cmtx;     % the constraint matrix
        cmtxInv;  % pseudo-inverse of the constraint matrix
        locToConstraint;
        
        subCmtxAndInvs;
        constraintVecXsectSubsets; %
        constraintVecSubsets; %
        scaleByOverlap = 1;
        
        overlapFraction;    % fraction of this location that
                            % overlaps with field-of-view
    end
    
    properties
        overlappingPatches;
        overlappingFull = 1;
        % aBigNumber = 1000;
    end
   
    methods
        
        function this = PatchConstraints( sz, f, overlappingPatches, scaleByOverlap, diffsz2d )
        % Constructor
            this.f    = f;
            if (~exist('overlappingPatches','var') || isempty(overlappingPatches))
                this.overlappingPatches = 0;
            else
                this.overlappingPatches = overlappingPatches;
            end
            if (exist('scaleByOverlap','var') && ~isempty(scaleByOverlap))
                this.scaleByOverlap = scaleByOverlap;    
            end
            
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
            
            if (exist('diffsz2d','var') && ~isempty(diffsz2d))
                this.sz2d = diffsz2d;
            end
            
            this.iniLocs();

        end
        
        function setOverlappingFull( this, overlappingFull )
           this.overlappingFull = overlappingFull;
           this.iniLocs();
        end
        
        function idx = locXyzDim2Idx( this, dim, xyz )
            num = length( dim );
            idx = zeros( num, 1 );
            for i = 1:num
                idx(i) = find( this.dimXyzList(:,1) == dim(i) & ...
                               this.dimXyzList(:,2) == xyz(i) );
            end
        end
        
        function [dim,xyz] = locIdx2XyzDim( this, idx )
            dim = this.dimXyzList( idx, 1 );
            xyz = this.dimXyzList( idx, 2 );
        end
        
        function buildCmtx( this )
            patchNumElem = prod( this.sz2d );
            numVariables = prod( this.sz3d );
            
            numPatchLocs   = size( this.dimXyzList, 1 );
            this.overlapFraction = zeros( numPatchLocs, 1);
            
            this.cmtx = zeros( this.numConstraints, numVariables );
            this.locToConstraint = false( numPatchLocs, this.numConstraints );
            this.constraintVecSubsets = false( numPatchLocs, this.numConstraints );
             
            k = 1;
            for i = 1 : numPatchLocs
                
                thisdim = this.dimXyzList(i,1);
                thisxyz = this.dimXyzList(i,2);
                
                rngStart = ( patchNumElem * (i-1) + 1);
                rng =  rngStart : rngStart + patchNumElem - 1;
                this.constraintVecSubsets( i , rng ) = true;
                
                % TODO - dont recompute this msk every time
                msk = Dict2dTo3d.planeMaskF( this.sz3d, thisxyz, thisdim, this.f );
                this.overlapFraction( i ) = nnz( msk ) ./ ( patchNumElem .* this.f );
                for j = 1:patchNumElem
                    
                    if( this.scaleByOverlap )
                        this.cmtx( k, (msk==j) ) = 1./this.overlapFraction(i);
                    else
                        this.cmtx( k, (msk==j) ) = 1;
                    end
                    
                    this.locToConstraint( i, k ) = 1;
                    
                    k = k + 1;
                end
            end
            
            this.cmtxInv = pinv( this.cmtx );
        end
               
        function [ H, f, A, bineq, Aeq, beq, lb, ub ] = buildQuadProg( this, b, constrainScalesPos, ...
                                                    constrainHrMinMax )
            f = [];
            if( ~exist('constrainScalesPos','var') || isempty( constrainScalesPos ))
                constrainScalesPos = true;
            end
            
            if( ~exist('constrainHrMinMax','var') || isempty( constrainHrMinMax ))
                constrainScalesPos = false;
            end
            
            
            patchNumElem = prod( this.sz2d );   
            numVariables = prod( this.sz3d );
            numParams = numVariables + this.numLocs; 
            
            % put hard constraint on sum of HR pixel values
            Aeq = zeros( 1, numParams );
            Aeq( this.numLocs + 1 : end ) = 1;
            beq = 1;
            
            A = [];
            bineq = [];
            
            aBigNumber = 1000;
            
            if ( constrainScalesPos || constrainHrMinMax )
                lb = -aBigNumber .* ones( numParams, 1 );
                ub =  aBigNumber .* ones( numParams, 1 );
            else
                lb = [];
                ub = [];
            end
            
            % ensure scales are positive
            if( constrainScalesPos )
%                 A = zeros( this.numLocs, numParams);
%                 A( 1:this.numLocs, 1:this.numLocs ) = -eye( this.numLocs );
%                 bineq = zeros( this.numLocs, 1);

                lb( 1:this.numLocs )  = 0;
            end
            

%             min_b = min( b( this.numLocs + 1 : end) );
            min_b = 0;
            max_b = max( b( this.numLocs + 1 : end) );
            if( constrainHrMinMax )
                
%                 AHrMinMax = zeros( 2 .* numVariables, numParams );
%                 AHrMinMax( :, this.numLocs + 1 : end ) = [ -eye(numVariables); eye(numVariables) ];
% 
%                 A = [ A; AHrMinMax ];
%                 bHrMinMax = [ repmat( min_b, numVariables, 1); repmat( max_b, numVariables, 1) ];
%                 bineq = [ bineq; bHrMinMax ];
                
                lb( this.numLocs + 1 : end ) = min_b;
                ub( this.numLocs + 1 : end ) = max_b;
                
            end
            
            % build H matrix
            H = zeros( numParams, numParams );
            Htemplate = PatchConstraints.buildQuadConstraintTemplate( this.f );
            
            k = 1;
            for i = 1 : this.numLocs
                
                thisdim = this.dimXyzList(i,1);
                thisxyz = this.dimXyzList(i,2);
                
                msk = Dict2dTo3d.planeMaskF( this.sz3d, thisxyz, thisdim, this.f );
                
                for j = 1:patchNumElem
                    
                    mskIdx = find( msk == j ) + this.numLocs;
                    Hsub = PatchConstraints.fillInQuadConstraintTemplate( Htemplate, b(k));
                    HsubRng = [ i; mskIdx ];
                    M = length( HsubRng );
                    H( HsubRng, HsubRng ) = H( HsubRng, HsubRng ) + Hsub(1:M, 1:M);
                    
                    k = k + 1;
                end % loop over pixels per patch
            end % loop over patch locations
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
            binds    = false( 1, this.numConstraints );
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
            xsectList = logical(eye( numPatchLocs ));
            
            for i = 1:numPatchLocs
                xsectList(i,:) = (this.dimXyzList(:,1) ~= this.dimXyzList(i,1));
            end
            
            this.xsectList = xsectList;
        end
        
        function xsectList = buildIntersectionListOL( this )
            
            xsectList = logical(eye( this.numLocs ));
            for i = 1:this.numLocs
                for j = i+1:this.numLocs
                    if( this.doOverlap( this.dimXyzList(i,:), this.dimXyzList(j,:)))
                        xsectList(i,j) = true;
                        xsectList(j,i) = true;
                    end
                end
            end
        end
        
        function areOverlapping = doOverlap( this, dimxyz1, dimxyz2 )
        % dimxyz1 and dimxyz2 are 2-vectors where
        % [dim xyz ]
        % and dim is the dimension of the normal vector
        % and xyz is the offset coordinate in that dimension
        
            areOverlapping = ( dimxyz1(1) ~= dimxyz2(1) || ...
                               abs(dimxyz1(2) - dimxyz2(2)) < this.f );
                             
        end
           
        function b = constraintValue( this, patchMtx, obj )
            if( isa( obj, 'net.imglib2.algorithms.opt.astar.SortedTreeNode'))
                b = this.constraintValueNode( patchMtx, obj );
            else
                b = this.constraintValueList( patchMtx, obj );
            end
        end
        
        function b = constraintValueNode( this, patchMtx, node )
            
            patchNumElem = prod( this.sz2d ); % num constraints per patch
            depth = node.getDepth();
            b = zeros( patchNumElem * (depth+1), 1 );
             
            for i = 0:depth

                dat = node.getData();
                j = this.locXyzDim2Idx( dat.dim, dat.xyz );
                idx = dat.idx;

                start = patchNumElem * (j - 1) + 1;
                rng = start : start + patchNumElem - 1;

                b( rng ) = patchMtx( idx, : );

                node = node.getParent();
            end
             
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
                if( isempty( i ))
                    continue;
                end
                    
                rowPermutation(r) = i;
            end
            rowPermutation = rowPermutation(rowPermutation~=0);
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
                    % if idxList is a column vector, assume that
                    % the indices are given in the same order as 
                    % this.dimXyzList
                    idx =  idxList( i );
                elseif( size( idxList, 2) == 3 )
                    j = this.locXyzDim2Idx( idxList(i,1), ...
                                            idxList(i,2));
                                        
                    idx = idxList( i, 3 );
                    
                    start = patchNumElem * (j - 1) + 1;
                	brng = start : start + patchNumElem - 1;
                    
                end
                
                b(brng) = patchMtx( idx, : );
                brng = brng + patchNumElem;
            end
        end
        
        function bnew = updateConstraints( this, patchMtx, b, jList, idx )
        % bnew = updateConstraints( this, patchMtx, b, j, idx )
        %   patchMtx - N x M matrix where N is number of dictionary elements
        %   b        - current b vector
        %   j        - index of patch location (into dimXyzList)
        %   idx      - index into patchMtx of replacing patch
            bnew = b;
            patchNumElem = prod( this.sz2d );
            
            for i = 1:length(jList)
                j = jList( i );
                start = patchNumElem * (j - 1) + 1;
                rng = start : start + patchNumElem - 1;
            
                bnew( rng ) = patchMtx( idx(i), : );
            end
        end
        
    end
    
    methods( Access = protected )
        
        function iniLocs( this )
            % Initialize Location parameters 
            if(  this.overlappingPatches )
                if( this.overlappingFull )
                    this.pairLocRng = (-this.f + 2) : this.sz3d(1);
                else
                    this.pairLocRng = (1) : this.sz3d(1)- this.f + 1;
                end
            else
                this.pairLocRng = (1) : this.f : this.sz3d(1);
            end
            
            [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
            this.dimXyzList = [ dList(:), xyzList(:) ];
            this.numLocs = size( this.dimXyzList, 1 );
            
            this.numConstraints = prod( this.sz2d ) *  this.numLocs;
            this.buildIntersectionList();
        end
    end
    
    methods( Static )
        
        function [ Hsub ] = fillInQuadConstraintTemplate( Htemplate, b )
            % Fills in the specific values of a submatrix given
            %
            if( ~isscalar( b ))
                error( 'b must be scalar' );
            end
            
            Hsub = Htemplate;
            Hsub( Htemplate == -1 ) =  -b;
            Hsub( Htemplate == -2 ) =  b.*b;
        end
        
        function [ Htemplate, ri, ci ] = buildQuadConstraintTemplate( f )
            
            % each constraint effects one scale and 'f' HR pixel intensities
            numParamsPerConstraint = f + 1;
            
            Htemplate = zeros( numParamsPerConstraint, numParamsPerConstraint );
            
            % the scale goes in the first row/column
            Htemplate( 1, 2:end ) = -1; % the pixel-scale interaction
            Htemplate( 2:end, 1 ) = -1; % the scale-pixel interaction
            Htemplate( 1, 1)   = -2;  % the scale-scale element is different
            Htemplate( 2:end, 2:end ) = 1; % the pixel-pixel interaction
        end
    end
end
