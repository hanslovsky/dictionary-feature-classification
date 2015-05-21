classdef PatchConstraintsGeneralZSparsePadShift < PatchConstraintsGeneralZSparse
    
    methods
        
        function this = PatchConstraintsGeneralZSparsePadShift( sz2d, sz3d, f, scaleByOverlap, ...
                minOverlapRatio, maxShiftFrac )
            % Constructor
            this = this@PatchConstraintsGeneralZSparse( sz2d, sz3d, f, scaleByOverlap, ...
                minOverlapRatio, maxShiftFrac );
            
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

        function b = constraintValueList( this, patchMtx, idxList, model )
            % idxList must be in the same order as this.dimXyzList
            
            fprintf('PatchConstraintsGeneralZSparsePadShift: constraintValueList\n');
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
        
    end
    
end
