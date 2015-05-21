classdef PatchConstraintsGeneralZ < PatchConstraints
    
    
    properties
        overlappingXY = 1;
        overlappingOutofBoundsXY = 1;
    end
    
    methods
        
        function this = PatchConstraintsGeneralZ( sz2d, sz3d, f, overlappingPatches, scaleByOverlap, ...
                    overlappingXY, overlappingOutofBoundsXY )
            % sz3d is the size of a low res patch
            % 
            
            % Constructor
            this = this@PatchConstraints( sz3d(1), f, overlappingPatches, scaleByOverlap, []  );
            
            if( ~exist('overlappingXY','var') || isempty(overlappingXY) )
                overlappingXY = 1;
            end
            if( ~exist('overlappingOutofBoundsXY','var') || isempty(overlappingOutofBoundsXY) )
                overlappingOutofBoundsXY = 1;
            end
            
            this.sz2d = sz2d;
            this.sz3d = sz3d;
            this.sz3dHR = sz3d .* [ 1 1 this.f ];
            
            this.overlappingXY = overlappingXY;
            this.overlappingOutofBoundsXY = overlappingOutofBoundsXY;
            
            if( ~isa( this, 'PatchConstraintsGeneralZ'))
                this.iniLocs();
            end
        end
        
        function setOverlapping( this, overlappingPatches, overlappingFull, overlappingXY )
            this.overlappingPatches = overlappingPatches;
            this.overlappingFull    = overlappingFull;
            this.overlappingXY      = overlappingXY;
            this.iniLocs();
        end
        
        function idx = locXyzDim2Idx( this, dim, xyz )
            num = length( dim );
            idx = zeros( num, 1 );
            for i = 1:num
                idx(i) = find( this.dimXyzList(:,1) == dim(i) & ...
                               this.dimXyzList(:,2) == xyz(i,1) & ...
                               this.dimXyzList(:,3) == xyz(i,2) & ...
                               this.dimXyzList(:,4) == xyz(i,3) );
            end
        end
        
        function [dim,xyz] = locIdx2XyzDim( this, idx )
            dim = this.dimXyzList( idx, 1 );
            xyz = this.dimXyzList( idx, 2:4 );
        end
            
        function buildCmtxWShifts( this )
            f = this.f;
            sz2d = this.sz2d;
            sz3d = this.sz3d;
            
            patchNumElem = prod( sz2d );
            numVariables = prod( sz3d );
            
            numPatchLocs   = size( this.dimXyzList, 1 );
            this.overlapFraction = zeros( numPatchLocs, 1);
            
            if( this.doSparseCmtx )
                this.cmtx = spalloc( this.numConstraints, numVariables, this.numConstraints .* this.f );
                this.constraintVecSubsets = logical( spalloc( numPatchLocs, this.numConstraints, numPatchLocs * patchNumElem ));
%                 this.locToConstraint      = logical( spalloc( numPatchLocs, this.numConstraints, numPatchLocs * patchNumElem ));
            else
                this.cmtx = zeros( this.numConstraints, numVariables );
%                 this.locToConstraint = false( numPatchLocs, this.numConstraints );
                this.constraintVecSubsets = false( numPatchLocs, this.numConstraints );
                this.constraintToHrLoc    = false( numVariables, this.numConstraints );
            end
            
            baseSz = sz3d + 2.*( f - 1);
            baseMask = zeros( baseSz );
            
            N = prod( sz2d );
            baseMask( 1:sz2d(1), 1:sz2d(2), 1:f ) = ...
                repmat(reshape( 1:N, sz2d ), [1 1 f]);
            
            locIdx = 1;
            for x = 0:sz3d(1)
                for y = 0:sz3d(2)
                    for z = 0:sz3d(3)
                        
                        tmpmsk = circshift( baseMask, [(x) (y) (z)] );
                        mskPt = tmpmsk( f:f+sz3d(1)-1, f:f+sz3d(2)-1, f:f+sz3d(3)-1);
                        
                        for d = 1:3 % dimensions
                            
                            if( d == 1 )
                                msk = shiftdim( mskPt, 2 );
                            elseif( d == 2)
                                msk = permute( mskPt, [1 3 2] );
                            else
                                msk = mskPt;
                            end
                            
                            this.cmtx( locIdx, (msk > 0) ) = 1;
                            this.overlapFraction( locIdx ) = nnz( msk );
                            
%                             this.constraintVecSubsets( locIdx, rng ) = true;
%                             this.locToConstraint( locIdx, 1 ) = true;
                            
                            locIdx = locIdx + 1;
                            
                        end %d
                    end % z
                end % y
            end % x

            
        end
        
        function buildCmtx2( this )
            patchNumElem = prod( this.sz2d );
            numVariables = prod( this.sz3d );
            
            numPatchLocs   = size( this.dimXyzList, 1 );
            this.overlapFraction = zeros( numPatchLocs, 1);
            
            if( this.doSparseCmtx )
                this.cmtx = spalloc( this.numConstraints, numVariables, this.numConstraints .* this.f );
                this.constraintVecSubsets = logical( spalloc( numPatchLocs, this.numConstraints, numPatchLocs * patchNumElem ));
                this.locToConstraint      = logical( spalloc( numPatchLocs, this.numConstraints, numPatchLocs * patchNumElem ));
            else
                this.cmtx = zeros( this.numConstraints, numVariables );
                this.locToConstraint = false( numPatchLocs, this.numConstraints );
                this.constraintVecSubsets = false( numPatchLocs, this.numConstraints );
                this.constraintToHrLoc    = false( numVariables, this.numConstraints );
            end
            
            k = 1;
            for i = 1 : numPatchLocs
                if( mod( i , 10 ) == 0 )
                    fprintf('location %d of %d\n', i, numPatchLocs);
                end
                thisdim = this.dimXyzList(i,1);
                thisxyz = this.dimXyzList( i, 2:4 );
                
                for j = 1:patchNumElem
                    msk = this.planeMaskFast( thisdim, thisxyz, j );
                    this.cmtx( k,  msk(:) ) = true;
                    k = k + 1;
                end
            end
        end
        
        function buildCmtxFastSparseIHope( this )
            patchNumElem = prod( this.sz2d );
            numVariables = prod( this.sz3d );
            
            numPatchLocs   = size( this.dimXyzList, 1 );
            this.overlapFraction = zeros( numPatchLocs, 1);
            
            offsets = [ 1 (sz3d(2)) prod(sz3d(2:3)) ];
            additions = (this.dimXyzList(:,2:4) - 1) * offsets';
             
        end
        
        function [ basei, basej ] = buildBaseCmtx( this, dim )
            
            patchNumElem = prod( this.sz2d );
            numVariables = prod( this.sz3d );
            numPatchLocs   = size( this.dimXyzList, 1 );
            
            offsets = [ 1 (this.sz3d(2)) prod(this.sz3d(2:3)) ];
            
            switch( dim )
                case 1
                    % x normal
                    gridargs = {0, 0:this.sz2d(1)-1, 0:this.sz2d(2)-1 };
                case 2
                    % y normal
                    gridargs = {0:this.sz2d(1)-1, 0, 0:this.sz2d(2)-1};
                case 3
                    % z normal
                    gridargs = {0:this.sz2d(1)-1, 0:this.sz2d(2)-1, 0 };
                otherwise
                    error( 'invalid dim' )
            end
            
            [ gridx, gridy, gridz] = ndgrid( gridargs{:} );
            additions = [ gridx(:) gridy(:) gridz(:) ] * offsets';
            msk = this.planeMaskFast( dim, [ 1 1 1 ], 1 );
            
            basej = bsxfun( @plus, repmat( find(msk)',patchNumElem,1), additions );
            basei = bsxfun( @times, ones( size( basej )), (1:size(basej,1))' );
        end
        
        function buildCmtx( this )
            patchNumElem = prod( this.sz2d );
            numVariables = prod( this.sz3d );
            
            numPatchLocs   = size( this.dimXyzList, 1 );
            this.overlapFraction = zeros( numPatchLocs, 1);
            
            if( this.doSparseCmtx )
%                 fprintf('allocating cmtx with %d elements\n', (this.numConstraints .* this.f));
%                 fprintf('allocating constraintVecSubsets with %d elements\n', (numPatchLocs * patchNumElem));
%                 fprintf('allocating locToConstraint with %d elements\n', (numPatchLocs * patchNumElem));
                this.cmtx = spalloc( this.numConstraints, numVariables, this.numConstraints .* this.f );
                this.constraintVecSubsets = logical( spalloc( numPatchLocs, this.numConstraints, numPatchLocs * patchNumElem ));
                this.locToConstraint      = logical( spalloc( numPatchLocs, this.numConstraints, numPatchLocs * patchNumElem ));
            else
                this.cmtx = zeros( this.numConstraints, numVariables );
                this.locToConstraint = false( numPatchLocs, this.numConstraints );
                this.constraintVecSubsets = false( numPatchLocs, this.numConstraints );
                this.constraintToHrLoc    = false( numVariables, this.numConstraints ); 
            end
            
            
            k = 1;
            for i = 1 : numPatchLocs
                
                if( mod( i , 10 ) == 0 )
                   fprintf('location %d of %d\n', i, numPatchLocs); 
                end
                thisdim = this.dimXyzList(i,1);
                thisxyz = this.dimXyzList( i, 2:4 );
                
                rngStart = ( patchNumElem * (i-1) + 1);
                rng =  rngStart : rngStart + patchNumElem - 1;
                this.constraintVecSubsets( i , rng ) = true;
                
                % TODO - dont recompute this msk every time
                msk = this.planeMask( thisdim, thisxyz, this.f );
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
        end
        
        function computeCmtxInverse( this )
            this.cmtxInv = pinv( this.cmtx );
        end
        
        function b = constraintValueList( this, patchMtx, idxList, model )
            % idxList must be in the same order as this.dimXyzList
            
            fprintf('PatchConstraintsGeneralZ: constraintValueList\n');
            if( ~exist( 'model','var') )
                model = [];
            end
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
            brng = 1:patchNumElem;
            
            b = zeros( prod(this.sz3d), 1 );
            for i = 1:N
                
                if( ndim <= 1)
                    % if idxList is a column vector, assume that
                    % the indices are given in the same order as
                    % this.dimXyzList
                    idx =  idxList( i );
                elseif( size( idxList, 2 ) == size( patchMtx, 1 ))
                    idx = idxList( i, : )';
                else
                    error( 'invalid index type' );
                end % checking idx
                
                imsk = this.planeMaskI( i );
                pmsk = imsk( imsk > 0 );
                bsubrng = brng( pmsk );
                
                % add the constraints
                if( length( idx ) > 1 )
                    b( bsubrng  ) = patchMtx( :, pmsk )' * idx;
                else
                    
                    b( bsubrng  ) = patchMtx( idx, pmsk );
                end
                
                if( ~isempty( model ) && ~isempty(model{i}))
                    b( bsubrng ) = feval( model{i},  b( bsubrng ));
                end
                
                brng = brng + patchNumElem;
            end
        end
        
        function b = constraintValueListMin( this, patchMtx, idxList, model )
            % idxList must be in the same order as this.dimXyzList
            
            fprintf('PatchConstraintsGeneralZ: constraintValueList\n');
            if( ~exist( 'model','var') )
                model = [];
            end
            if( islogical( idxList ))
                idxList = find( idxList );
            end
            
            ndim  = nnz( size(idxList) > 1);
            if( ndim == 1 )
                N = numel(idxList);
            else
                N = size(idxList,1);
            end
            
            b = zeros( prod(this.sz3d), 1 );
            bstart = 1;
            for i = 1:N
                
                if( ndim <= 1)
                    % if idxList is a column vector, assume that
                    % the indices are given in the same order as
                    % this.dimXyzList
                    idx =  idxList( i );
                elseif( size( idxList, 2 ) == size( patchMtx, 1 )) 
                    idx = idxList( i, : )';
                else
                    error( 'invalid index type' );
                end % checking idx
                
                imsk = this.planeMaskI( i );
                
                % see how many constraints are contributed by this sub-patch
                numConstraintsHere = nnz( imsk );
                bend = bstart + numConstraintsHere - 1;
                
                % add the constraints
                b( bstart:bend ) = patchMtx( idx, imsk( imsk > 0 ));
                
                if( ~isempty( model ) && ~isempty(model{i}))
                     b( bstart:bend ) = feval( model{i},  b( bstart:bend ));
                end
                
                bstart = bstart + numConstraintsHere;
            end
        end
        
        function msk = planeMaskI( this, i, centered )
            if( ~exist('centered','var'))
                centered = [];
            end
            msk = this.planeMask( this.dimXyzList(i,1), this.dimXyzList(i,2:4), this.f, centered );
        end
        
        function [ dsXmtx, dsYmtx ] = buildProjectionMatrix2d( this )
                 
            [dsxi,dsxj] = ndgrid( 1:this.f:this.sz2d(1), 1:this.sz2d(2));
            [dsyi,dsyj] = ndgrid( 1:this.sz2d(1), 1:this.f:this.sz2d(2));
            
            N = prod( this.sz2d );
            M = N ./ this.f;
            
            dsXmtx = zeros( M, N );
            dsYmtx = zeros( M, N );
            
            for i = 1:M %#ok<BDSCI>
                dsxk = sub2ind( this.sz2d, dsxi(i):dsxi(i)+this.f-1, repmat( dsxj(i), 1, this.f));
                dsXmtx( i, dsxk ) = 1; 
                
                dsyk = sub2ind( this.sz2d, repmat( dsyi(i), 1, this.f), dsyj(i):dsyj(i)+this.f-1 );
                dsYmtx( i, dsyk ) = 1;
            end
        end
        
        function [obs,msk,mtx,oinds] = getSubPatchI( this, im, loci, isLR )
            [obs,msk,mtx,oinds] = this.getSubPatch( im, this.dimXyzList(loci,1), ...
                         this.dimXyzList(loci,2:end), isLR );
        end
        
        function [obs, msk, mtx, omsk] = getSubPatch( this, im, dim, xyz, isLR, centered, permloc )
            % returns an observation at the given location
            % also returns a mask indicating which part of the 2d
            % dictionary this observation should be compared against
            
            if( ~exist( 'centered', 'var'))
                centered = [];
            end
            if( ~exist( 'permloc', 'var') || isempty(permloc))
                permloc = false;
            end
            
            % TODO - make centered option work ( low priority )
            sz2d = this.sz2d;
            
            start = xyz;
            obsSzRaw = [ sz2d, this.f ];
            
            % always end up downsampling along y
            [~, mtx ] = this.buildProjectionMatrix2d();
            
            zdim = 3;
            switch( dim )
                case 1
                    perm = [ 2 3 1 ];
%                     iperm = [3 1 2];
                    zdim = 2;
                case 2
                    perm = [ 1 3 2 ];
%                     iperm = [ 1 3 2 ];
                    zdim = 2;
                case 3
                    perm = [1 2 3];
%                     iperm = [1 2 3];
                    mtx = eye( prod( sz2d ));
                otherwise
                    error( 'invalid dimension')
            end
            
            % permute or ipermute ?
            imp = permute( im, perm );
            
            startInds = start( perm );
            
            if( isLR )
               obsSzRaw( zdim ) = obsSzRaw( zdim )./ this.f;
            end
            endInds = startInds + obsSzRaw - 1;
            
            sIndsOrig = startInds;
            eIndsOrig = endInds;
            startInds = max( startInds, [1 1 1]);
            endInds   = min( endInds, size( imp ));
            
            mskStrtInds =   1 + ( startInds(1:2) - sIndsOrig(1:2) );
            mskEndInds = sz2d - ( eIndsOrig(1:2) -   endInds(1:2) );
            msk = false( sz2d );
            msk( mskStrtInds(1):mskEndInds(1), ...
                 mskStrtInds(2):mskEndInds(2)) = true;
             
             if( ~isempty( mtx ))
                 mtx = mtx( :, msk );
                 
                 % remove rows with all zeros
                 mtx = mtx( any(mtx,2), : );
                 
                 % and old thing that seems not to work
                 % mtx( :, msk ) = 0;
             end
            
            obs = imp( startInds(1):endInds(1), ...
                       startInds(2):endInds(2), ...
                       startInds(3):endInds(3));
            obs = sum( obs, 3 );
            
            if( isLR )
                sz3d = this.sz3d;
            else
                sz3d = this.sz3dHR;
            end
            
            omsk = false( sz3d );
            omsk( startInds(1):endInds(1), ...
                  startInds(2):endInds(2), ...
                  startInds(3):endInds(3)) = true;
            omsk = ipermute( omsk, perm );
            oinds = find( omsk );
                   
        end
        
        function im = planeMask( this, d, xyz, f, centered )
            
            sz2d = this.sz2d;
            sz3d = this.sz3d;
            
            im = zeros( sz3d );
            
            if( ~exist('centered','var') || isempty( centered ))
                centered = false;
            end
            
            if( centered )
                half = (f-1)./2;
                zrng = xyz(3)-half : xyz(3)+half;
            else
                zrng = xyz(3) : xyz(3) + f - 1;
            end
            
            validZ = (zrng>0) & (zrng<=sz3d(3));
            zrng = zrng( validZ );
            
            N = prod(sz2d);
            v = repmat( reshape(1:N, sz2d), [1 1 f]);
            
            xd = xyz(1) : xyz(1) + sz2d(1) - 1;
            yd = xyz(2) : xyz(2) + sz2d(2) - 1;
            
            validX = (xd > 0) & (xd <= sz3d(1));
            validY = (yd > 0) & (yd <= sz3d(2));
            
            xd = xd( validX );
            yd = yd( validY );
            
            v = v( validX, validY, validZ );
            
            switch d
                case 1
                    %                     size( im( zrng, xd, yd ))
                    %                     size(permute( v, [3 1 2]))
                    im( zrng, xd, yd ) = permute( v, [3 1 2]);
                case 2
                    im( xd, zrng, yd ) = permute( v, [1 3 2]);
                case 3
                    im( xd, yd, zrng ) = v;
                otherwise
                    error('invalid normal direction');
            end
        end
        
        function im = planeMaskLR( this, d, xyz, f, centered )
            
            sz3d = this.sz3d;
            im = zeros( sz3d );
            
            if( ~exist('centered','var') || isempty( centered ))
                centered = false;
            end
            
            [ xi, yi, zi ] = this.getSubImageRangeLR( [d, xyz] );
            
            args = { xi, yi, zi };
            vv = true(3,1);
            vv( d ) = false;
            
            M = prod(cellfun( @length, args(vv)));
            
            N = length( args{d} );
            for i = 1:N
                tmpargs = args;
                tmpargs{ d } = args{ d }( i ); 
                sz = cellfun( @length, tmpargs );
                im( tmpargs{:} ) = reshape( 1:M, sz );
            end
        end
        
        function im = planeMaskLRI( this, i, centered )
            if( ~exist('centered','var') || isempty( centered ))
                centered = false;
            end
            [dim,xyz] = this.locIdx2XyzDim( i );
            im = this.planeMaskLR( dim, xyz, this.f, centered );
        end
        
        function msk = planeMaskFast( this, dim, xyz, k )
            msk = false( this.sz3d );
            xo = (1:this.sz2d(1)) - 1;
            yo = (1:this.sz2d(2)) - 1;
            zo = (1:this.f) - 1;
            
            if( exist('k','var') && ~isempty( k ))
               [ i, j ] = ind2sub( this.sz2d, k );
               xo = xo( i ); 
               yo = yo( j );
            end
            
            switch dim
                case 1
                    msk( zo + xyz(1), xo + xyz(2), yo + xyz(3) ) = true;
                case 2
                    msk( xo + xyz(1), zo + xyz(2), yo + xyz(3) ) = true;
                case 3
                    msk( xo + xyz(1), yo + xyz(2), zo + xyz(3) ) = true;
                otherwise
                    error( 'invalid normal dimension' );
            end
        end
        
        function K = numLocs( this )
            l = this.sz2d(1);
            n = this.sz3d(1);
            D = this.f;
            K = 3*( ((l-1) + n).^2 .* (n+D-1));
        end
        
        function iniLocs( this )
            
            xStart = 1;
            yStart = 1;
            
            if( this.overlappingXY )
                if( this.overlappingOutofBoundsXY )
                    xStart = -this.sz2d(1)+2;
                    yStart = -this.sz2d(2)+2;
                    
                    xEnd = this.sz3d(1);
                    yEnd = this.sz3d(2);
                else
                    xEnd = this.sz3d(1);
                    yEnd = this.sz3d(2);
                end
            else
                %error( 'not yet implemented');
                xEnd = this.sz3d(1) - this.f + 1;
                yEnd = this.sz3d(2) - this.f + 1;
            end
            
            % Initialize Location parameters
            if(  this.overlappingPatches )
                if( this.overlappingFull )
                    this.pairLocRng = (-this.f + 2) : this.sz3d(3);
                else
                    this.pairLocRng = (1) : this.sz3d(1)- this.f + 1;
                end
            else
                this.pairLocRng = (1) : this.f : this.sz3d(3);
            end
            
            [ dList, x, y, z ] = ndgrid( 1:3, xStart:xEnd, yStart:yEnd, this.pairLocRng );
            this.dimXyzList = [ dList(:), x(:), y(:), z(:) ];
            this.numLocs = size( this.dimXyzList, 1 );
            
            this.numConstraints = prod( this.sz2d ) *  this.numLocs;
        end
        
    end
    
    methods( Static )
        
        function msk = planeMaskFStatic( sz2d, sz3d, dim, xyz, f, centered )
            % xyz is a 3-vector
            % dim is {1,2,3}
            
            msk = zeros( sz3d );
            
            if( ~exist('centered','var') || isempty( centered ))
                centered = false;
            end
            
            if( centered )
                half = (f-1)./2;
                zrng = xyz(3)-half : xyz(3)+half;
            else
                zrng = xyz(3) : xyz(3) + f - 1;
            end
            
            validZ = (zrng>0) & (zrng<=sz3d(3));
            zrng = zrng( validZ );
            
            N = prod(sz2d);
            v = repmat( reshape(1:N, sz2d), [1 1 nnz(validZ)]);
            
            xd = xyz(1) : xyz(1) + sz2d(1) - 1;
            yd = xyz(2) : xyz(2) + sz2d(2) - 1;
            
            validX = (xd > 0) & (xd <= sz3d(1));
            validY = (yd > 0) & (yd <= sz3d(2));
            
            xd = xd( validX );
            yd = yd( validY );
            
            xr = 1:sz2d(1);
            yr = 1:sz2d(2);
            
            xr = xr( validX );
            yr = yr( validY );
            
            m3 = true( sz3d );
            m2 = true( sz2d );
            
            switch dim
                case 1
                    msk( ZR3, XR3, YR3 ) = permute( v(XR2,YR2), [3 1 2]);
                case 2
                    msk( XR3, ZR3, YR3 ) = permute( v(XR2,YR2), [1 3 2]);
                case 3
                    msk( XR3, YR3, ZR3 ) = v(XR2,YR2);
                otherwise
                    error('invalid normal direction');
            end
        end
        
        function msk = planeMaskFOld( sz2d, sz3d, dim, xyz, f, centered )
            % xyz is a 3-vector
            % dim is {1,2,3}
            
            msk = zeros( sz3d );
            
            if( ~exist('centered','var') || isempty( centered ))
                centered = false;
            end
            
            if( centered )
                half = (f-1)./2;
                zrng = xyz(3)-half : xyz(3)+half;
            else
                zrng = xyz(3) : xyz(3) + f - 1;
            end
            
            valid = (zrng>0) & (zrng<=sz3d(3));
            zrng = zrng( valid );
            
            N = prod(sz2d);
            v = repmat( reshape(1:N, sz2d), [1 1 f]);
            
            %             [ x, y ] = ndgrid( 1:sz2d(1), 1:sz2d(2) );
            xd = xyz(1) : xyz(1) + sz2d(1) - 1;
            yd = xyz(2) : xyz(2) + sz2d(2) - 1;
            
            validX = (xd > 0) & (xd <= sz3d(1));
            validY = (yd > 0) & (yd <= sz3d(2));
            
            xd = xd( validX );
            yd = yd( validY );
            
            v = v( repmat(validX',1,2) & repmat(validY,2,1) );
            
            switch dim
                case 1
                    msk( zrng, xd, yd ) = v; % permute(v, [3 1 2]);
                case 2
                    msk( xd, zrng, yd ) = v; % permute(v, [1 3 2]);
                case 3
                    msk( xd, yd, zrng ) = v;
                otherwise
                    error('invalid normal direction');
            end
        end
        
    end
end
