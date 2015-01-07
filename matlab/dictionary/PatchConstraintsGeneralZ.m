classdef PatchConstraintsGeneralZ < PatchConstraints
    
    properties
        overlappingXY = 1;
        overlappingOutofBoundsXY = 1;
    end
    
    methods
        
        function this = PatchConstraintsGeneralZ( sz2d, sz3d, f, overlappingPatches, scaleByOverlap  )
            % Constructor
            this = this@PatchConstraints( sz3d(1), f, overlappingPatches, scaleByOverlap, []  );
            
            this.sz2d = sz2d;
            this.sz3d = sz3d;
            
            this.iniLocs();
        end
        
        function setOverlapping( this, overlappingPatches, overlappingFull, overlappingXY )
            this.overlappingPatches = overlappingPatches;
            this.overlappingFull    = overlappingFull;
            this.overlappingXY      = overlappingXY;
            this.iniLocs();
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
        
        function im = planeMask( this, d, xyz, f, centered)
            
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
                xEnd = 1;
                yEnd = 1;
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
