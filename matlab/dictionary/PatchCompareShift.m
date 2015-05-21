classdef PatchCompareShift < PatchCompare
    % Only deals with 2d for now
    
    properties
        
        compOverUnionOnly = 0;
        padmethod = 'replicate';
    end
    
    properties( SetAccess = protected )
        maxshift;
        shifts;
        numShifts;
        noShiftIndex;
    end
    
    methods
        
        function this = PatchCompareShift( method, maxshift )
            this = this@PatchCompare( method );
            this.maxshift = maxshift;
            
            [ xs, ys ] = ndgrid( -maxshift(1) : maxshift(1), ...
                                 -maxshift(2) : maxshift(2));
            this.shifts = [xs(:) ys(:)];
            this.numShifts = size( this.shifts, 1 );
            
            this.noShiftIndex = find( all(this.shifts == 0, 2) );
        end
        
        function [ im_shift, msk_center ] = shiftImageOld( this, im, shiftIdx, sz )
            if( exist( 'sz', 'var') && ~isempty( sz ))
                im = reshape( im, sz );
            else
                sz = size( im );
            end
            msk_center  = padarray( true(sz), this.maxshift, false );
            im_shift = padarray( reshape(im,sz), this.maxshift, this.padmethod );
            im_shift = circshift( im_shift, this.shifts( shiftIdx, : ));
        end
        
        function [ im_shift, msk_center ] = shiftImage( this, im, shiftIdx, sz )
            if( exist( 'sz', 'var') && ~isempty( sz ))
                im = reshape( im, sz );
            else
                sz = size( im );
            end
            msk_center  = padarray( true(sz), this.maxshift, false );
            im_shift = padarray( reshape(im,sz), this.maxshift, this.padmethod );
            im_shift = circshift( im_shift, this.shifts( shiftIdx, : ));
        end
        
        function [ dists, shifts ] = shiftDistances( this, x, Y, sz, isRev )
            % Assumes vector input for x
            % and matrix input for Y [ numObservations x numVariables ]
            % As such, sz is a required input
            
            if( ~exist('isRev', 'var') || isempty( isRev ))
                isRev = false;
            end
            
            numComps = size( Y, 1 );
            dists  = zeros( numComps, 1 );
            shifts = zeros( numComps, 2 );
            for i = 1 : numComps
                if( isRev )
                    [d,s] = this.shiftDistance( x, Y(i,:), sz );
                else
                    [d,s] = this.shiftDistance( Y(i,:), x, sz );
                end
                dists( i ) = d;
                shifts( i , : ) = s;
            end
        end
        
        function [ dist, shift ] = shiftDistance( this, x, y, sz )
            % shift x and hold y constant
            % only compare values in an overlap
            
            issz = exist('sz','var') && ~isempty(sz);
            
            % build a mask
            if( issz )
                msk  = padarray( true(sz), this.maxshift, false );
                xpad = padarray( reshape(x,sz), this.maxshift, this.padmethod );
                ypad = padarray( reshape(y,sz), this.maxshift, nan );
            else
                msk  = padarray( true(size(x)), this.maxshift, false );
                xpad = padarray( x, this.maxshift, this.padmethod );
                ypad = padarray( y, this.maxshift, nan );
            end
            dist = inf;
            
            for i = 1 : this.numShifts
                
                xtmp = circshift( xpad, this.shifts( i, : ));
                
                % shift the mask
                if( this.compOverUnionOnly )
                    msktmp = circshift( msk, this.shifts( i, : ));
                    thismsk = ~isnan( msktmp ) & ~isnan( ypad );
                    if( ~nnz( thismsk ) )
                        fprintf('shiftDistance: warning...empty comparison\n');
                        continue;
                    end
                    
                    % compare values only within the shifted max
                    tmpdist = this.distance( xtmp( thismsk ), ypad( thismsk ));
                else
                    tmpdist = this.distance( xtmp( msk ), ypad( msk ));
                end
%                 % % For debug
%                 if( this.shifts( i, 1 ) == 0 &&  this.shifts( i, 2 ) == 0 )
%                    face = 1;
%                 end
%                 fprintf( 'i=%d : shift( %d %d )   %f v %f\n', i, this.shifts(i,:), tmpdist, dist );

                if( tmpdist < dist ) 
                    dist = tmpdist;
                    shift = this.shifts( i, : );
                elseif( tmpdist == dist && (norm(this.shifts(i,:)) < norm(shift)) )
                    dist = tmpdist;
                    shift = this.shifts( i, : );
                end
            end
        end % shiftDistance
    end % methods
    
end