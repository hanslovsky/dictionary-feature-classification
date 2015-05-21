classdef PatchCompareShiftPad < PatchCompareShift
    % Only deals with 2d for now ...

    properties( SetAccess = protected )
        bigSz;
        smallSz;
        msk_center;
    end
    
    methods
        
        function this = PatchCompareShiftPad( method, smallSz, bigSz  )
            this = this@PatchCompareShift( method, [0 0] );
            
            this.bigSz   = bigSz;
            this.smallSz = smallSz;
            
            this.maxshift = ( bigSz - smallSz )./2;
            this.msk_center = PatchCompareShiftPad.centerMask( smallSz, bigSz );
            
            [ xs, ys ] = ndgrid( -this.maxshift(1) : this.maxshift(1), ...
                                 -this.maxshift(2) : this.maxshift(2));
            this.shifts = [xs(:) ys(:)];
            this.numShifts = size( this.shifts, 1 );
        end
        
        function [ im_shift ] = shiftImage( this, im, shiftIdx, sz )
            if( exist( 'sz', 'var') && ~isempty( sz ))
                im = reshape( im, sz );
            else
                im = reshape( im, this.bigSz );
            end
            msk_cen = this.msk_center;
            msk = circshift( msk_cen, this.shifts( shiftIdx, : ));
            im_shift = im( msk );
        end
        
        function [ dists, shifts ] = shiftDistances( this, x, Y )
            % Assumes vector input for x
            % and matrix input for Y [ numObservations x numVariables ]
            % As such, sz is a required input
            
            numComps = size( Y, 1 );
            dists  = zeros( numComps, 1 );
            shifts = zeros( numComps, 2 );
            for i = 1 : numComps
                [d,s] = this.shiftDistance( x, Y(i,:) );
                dists( i ) = d;
                shifts( i , : ) = s;
            end
        end
        
        function [ dist, shift ] = shiftDistance( this, x, y )
            % shift y and hold x constant
            
            % build a mask
            if( ~isequal( this.smallSz, size(x) ))
                x = reshape(x,this.smallSz);
            end
            
            if( ~isequal( this.bigSz, size(y) ))
                y = reshape(y,this.bigSz);
            end
            
            msk_cen = this.msk_center;
            
            dist = inf; % we'd better find something smaller than infinity
            for i = 1 : this.numShifts
                
                % shift the mask
                msk = circshift( msk_cen, this.shifts( i, : ));
                
                % compare values only within the shifted max
                tmpdist = this.distance( x(:), y( msk ));
                
%                 % % For debug
%                 if( this.shifts( i, 1 ) == 0 &&  this.shifts( i, 2 ) == 0 )
%                    face = 1;
%                 end
%                 fprintf( 'i=%d : shift( %d %d )   %f v %f\n', i, this.shifts(i,:), tmpdist, dist );

                if( tmpdist < dist ) 
                    % if we have a lower distance, set it and the shift
                    dist = tmpdist;
                    shift = this.shifts( i, : );
                elseif( tmpdist == dist && (norm(this.shifts(i,:)) < norm(shift)) )
                    % border case - 
                    % if the distance is the same, compare magnitude of
                    % the shift, preferring smaller magnitude shifts
                    dist = tmpdist;
                    shift = this.shifts( i, : );
                end
            end
            
        end % shiftDistance
    end % methods
    
    methods( Static )
        
        function [ msk_center ] = centerMask( sz, szBig )
            dim = length( szBig );
            
            rad = ( sz - 1 ) ./ 2;
            radBig = ( szBig - 1 ) ./ 2;
            
            msk_center = false( szBig );
            midptBig = 1 + radBig;
            
            rngStart = midptBig - rad;
            rngEnd   = midptBig + rad;
            if( dim == 2 )
                msk_center( rngStart(1):rngEnd(1), rngStart(2):rngEnd(2) ) = true;
            elseif( dim == 3 )
                msk_center( rngStart(1):rngEnd(1), rngStart(2):rngEnd(2), ...
                            rngStart(3):rngEnd(3) ) = true;
            else
               error( 'only accepts 2 or 3 dimensional inputs'); 
            end
            
        end % centerMask
    end % static methods
end