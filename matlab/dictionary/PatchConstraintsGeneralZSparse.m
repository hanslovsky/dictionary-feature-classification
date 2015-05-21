classdef PatchConstraintsGeneralZSparse < PatchConstraintsGeneralZ
    
    properties
        minOverlapRatio = 0.25
        maxShiftFrac = 0.25;
    end
    
    properties ( SetAccess = protected )
        overlapAmount;
        shiftAmount;
        xySpacing;
        
        patchLocAdjustments;
    end
    
    properties ( SetAccess = immutable )
        dimXyzList_original;
    end
    
    methods
        
        function this = PatchConstraintsGeneralZSparse( sz2d, sz3d, f, scaleByOverlap, ...
                M, minOverlapRatio, maxShiftFrac )
            % Constructor
            this = this@PatchConstraintsGeneralZ( sz2d, sz3d, f, true, scaleByOverlap, ...
                true, false );
            
            if( exist( 'minOverlapRatio','var') && ~isempty( minOverlapRatio ))
                this.minOverlapRatio = minOverlapRatio;
            end
            if( exist( 'maxShiftFrac','var') && ~isempty( maxShiftFrac ))
                this.maxShiftFrac = maxShiftFrac;
            end
            
            this.overlappingFull = false; % set the default value to false
            
            this.iniLocs( M );
            this.dimXyzList_original = this.dimXyzList;
            this.resetAdjustments();
        end
        
        function resetLocs( this )
            % DEPRECATED
            this.dimXyzList = this.dimXyzList_original;
        end

        function modifyLocs( this, newLocs, i )
            % DEPRECATED
            if( ~exist('i','var') || isempty(i))
               i = 1:this.numLocs; 
            end
            if( length(i) ~= size( newLocs, 1 ))
                error('length of i and newLocs not compatible');
            end
            this.dimXyzList( i, 2:3 ) = newLocs;
        end
        
        function nudgeLocs( this, locDelta, i )
            % DEPRECATED!
            if( ~exist('i','var') || isempty(i))
                i = 1:this.numLocs;
            end
            if( length(i) ~= size( locDelta, 1 ))
                error('length of i and newLocs not compatible');
            end
            this.dimXyzList( i, 2:3 ) = this.dimXyzList( i, 2:3 ) + locDelta;
        end
        
        function resetAdjustments( this )
            this.patchLocAdjustments = zeros( size( this.dimXyzList, 1), 3 );
        end
        
        function updateLocAdjustments( this, locDelta, i )
            if( ~exist('i','var') || isempty(i))
                i = 1:this.numLocs;
            end
            if( length(i) ~= size( locDelta, 1 ))
                error('length of i and newLocs not compatible');
            end
            this.patchLocAdjustments( i, : ) = this.patchLocAdjustments( i, : ) - locDelta;
        end
        
        function setOverlapRatio( this, minOverlapRatio, maxShiftFrac )
            this.minOverlapRatio = minOverlapRatio;
            this.maxShiftFrac = maxShiftFrac;
            this.iniLocs();
        end

        function buildCmtx( this )
            fprintf('PatchConstraintsGeneralZSparse - buildCmtx\n');
            
            patchNumElem = prod( this.sz2d );
            numVariables = prod( this.sz3dHR );
            
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
        
        function buildCmtxTmp( this )
            
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
        
        function b = constraintValueList( this, patchMtx, idxList, model )
            % idxList must be in the same order as this.dimXyzList
            
            fprintf('PatchConstraintsGeneralZSparse: constraintValueList\n');
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
                
                dim = this.dimXyzList( i, 1 );
                % incorporate the nudge
                xyz = this.dimXyzList( i, 2:end ) - this.patchLocAdjustments( i, : );
                
                imsk = this.planeMask( dim, xyz, this.f );
                
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
            
            fprintf('PatchConstraintsGeneralZSparse: constraintValueListMin\n');
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
                
                dim = this.dimXyzList( i, 1 );
                % incorporate the nudge
                xyz = this.dimXyzList( i, 2:end ) - this.patchLocAdjustments( i, : );
                
                imsk = this.planeMask( dim, xyz );
                
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
        
        function iniLocs( this, M )
            % M - the number of subpatches to fit in the large patch
            % Enforces that no part of a subpatch is outside the larger
            % image
            
            sz2d = this.sz2d;
            mxsz = max( sz2d );
            sz3dHR = repmat( M*mxsz, 1, 3 );
            this.sz3dHR = sz3dHR;
            this.sz3d   = vecToRowCol(this.sz3dHR, 'row')./ [1 1 this.f ];
            
            psz = [this.sz2d, this.f ];
            rad2d = ( sz2d - 1 ) ./ 2;
            
            permmtxs = { [3 1 2 ],[1 3 2],[1 2 3]};
            locList = [];
            
            for dim = 1:3
                sz2dperm = psz( permmtxs{dim} );
                xrng = 1:sz2dperm(1):sz3dHR(1);
                yrng = 1:sz2dperm(2):sz3dHR(2);
                zrng = 1:sz2dperm(3):sz3dHR(3);
                [x,y,z] = ndgrid( xrng, yrng, zrng );
                locList = [ locList; dim.*ones( numel(x),1), x(:), y(:), z(:) ]; %#ok<AGROW> (its only twice)
            end
            this.dimXyzList = locList;
            
            this.numLocs = size( this.dimXyzList, 1 );
            this.numConstraints = prod( this.sz2d ) *  this.numLocs;
        end
        
        function iniLocsOld( this )
            
            this.overlapAmount = floor( this.sz2d .* this.minOverlapRatio );
            this.shiftAmount   = floor( this.sz2d .* this.maxShiftFrac  );
            this.xySpacing = this.sz2d - (  this.shiftAmount.^2 ) - this.overlapAmount;
            
            % try to even the sizes of 
            rad2d = ( this.sz2d - 1 ) ./ 2;
            m = mod( this.sz3dHR(1:2), this.sz2d );
            
            if( m(1) > 0  && m(1) < rad2d(1) + 1 )
                xrng = m(1) - rad2d(1) : this.sz2d(1) : this.sz3dHR(1);
            else
                xrng = 1 : this.sz2d(1) : this.sz3dHR(1);
            end
            if( m(2) > 0 && m(2) < rad2d(2) + 1 )
                yrng = m(2) - rad2d(2) : this.sz2d(2) : this.sz3dHR(2);
            else
                yrng = 1 : this.sz2d(2) : this.sz3dHR(2);
            end
            
            % Initialize Location parameters
            this.pairLocRng = (1) : this.f : this.sz3dHR(3);
            
            [ dList, x, y, z ] = ...
                    ndgrid( 1:3, xrng, yrng, this.pairLocRng );
                                     
            this.dimXyzList = [ dList(:), x(:), y(:), z(:) ];
            
            this.numLocs = size( this.dimXyzList, 1 );
            this.numConstraints = prod( this.sz2d ) *  this.numLocs;
        end
        
        function im = planeMask( this, d, xyz, f, centered )
            fprintf('PatchConstraintsGeneralZSparse - planeMask\n');
            sz2d = this.sz2d;
            sz3d = this.sz3dHR;
            
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
            size( v )
            
            switch d
                case 1
                    % size( im( zrng, xd, yd ))
                    % size(permute( v, [3 1 2]))
                    im( zrng, xd, yd ) = permute( v, [3 1 2]);
                case 2
                    im( xd, zrng, yd ) = permute( v, [1 3 2]);
                case 3
                    im( xd, yd, zrng ) = v;
                otherwise
                    error('invalid normal direction');
            end
        end
        
        function [ pmsk ] = patchNudgeMask( this, nudge )
            
            if( ~isempty( nudge ) && sum(abs(nudge))~=0 )
                xstart = 1;
                xend   = this.sz2d(1);
                
                ystart = 1;
                yend   = this.sz2d(2);
                
                if( nudge(1) > 0 )
                    xstart = xstart + nudge(1);
                end
                if( nudge(2) > 0 )
                    ystart = ystart + nudge(2);
                end
                
                if( nudge(1) < 0 )
                    xend = xend + nudge(1);
                end
                if( nudge(2) < 0 )
                    yend = yend + nudge(2);
                end
                pmsk = false(this.sz2d );
                pmsk( xstart:xend, ystart:yend ) = true;
            else
                pmsk = true(this.sz2d );
            end
        end
        
        function [ D_shifts ] = buildShiftedDictionary( this, D, shifter )
            if( ~isa( shifter, 'PatchCompareShift' ))
                error( 'input must be of type ''PatchCompareShift''');
            end
            
            numDict = size( D, 1 );
            numShft = shifter.numShifts;
            D_shifts = cell( numShft, 1 );
            
            for j = 1:numShft
                for i = 1:numDict
                    
                    if( ~isa( shifter, 'PatchCompareShiftPad'))
                        d = reshape( D(i,:), this.sz2d );
                        [ ds, mc ] = shifter.shiftImage( d, j );
                    else
                        % UGLY HACK!
                        % but at least it avoids recoding this method
                        [ ds ] = shifter.shiftImage( D(i,:), j );
                        mc = true( size( ds )); 
                    end
                    
                    if( i == 1 )
                        dict_ds = zeros( numDict, nnz( mc ));
                    end
                    dict_ds( i, : ) = ds( mc );
                end
                D_shifts{ j } = dict_ds;
            end
            
        end
        
        function [ D_shifts ] = buildShiftedDownsampledDictionary( this, D, shifter )
            if( ~isa( shifter, 'PatchCompareShift' ))
                error( 'input must be of type ''PatchCompareShift''');
            end
            
            numDict = size( D, 1 );
            numShft = shifter.numShifts;
            D_shifts = cell( numShft, 1 );
            
            for j = 1:numShft
                for i = 1:numDict
                    
                    d = reshape( D(i,:), this.sz2d );
                    [ ds, mc ] = shifter.shiftImage( d, j );
                    ds = reshape(ds( mc ), this.sz2d );
                    ddown = PatchConstraints.downsamplePatch( ds, this.sz2d, this.f );
                    if( i == 1 )
                        dict_ds = zeros( numDict, numel( ddown ));
                    end
                    dict_ds( i, : ) = ddown(:);
                end
                D_shifts{ j } = dict_ds;
            end
            
        end
        
        function [ I, imsk ] = fill3dWith2d( this, locIdx, patch )
            
            dim = this.dimXyzList( locIdx, 1 );
            
            % incorporate the nudge
            xyz = this.dimXyzList( locIdx, 2:end ) - this.patchLocAdjustments( locIdx, : );
            %             fprintf( 'xyz: %d %d %d\n', xyz );
            
            I = zeros( this.sz3d );
            imsk = this.planeMask( dim, xyz, this.f );
            I( imsk > 0 ) = patch( imsk( imsk > 0 ));
            
            %             pmsk = this.patchNudgeMask( this.patchLocAdjustments( locIdx, : ) );
            %             ii = find( pmsk );
            %             I( imsk > 0 ) = patch( ii( imsk( imsk > 0 )));
            
        end
    end
    
end
