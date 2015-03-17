classdef PatchCompare < handle
   
    properties( Constant  )
        methodOptions = { 'euc', 'sad', 'ncc' };
        % euc - euclidean distance 
        % sad - sum absolute distances
        % ncc - normalized cross correlation
    end
    
    properties
        method = 'euc';  % euclidean norm
    end
    
    methods
        
        function this = PatchCompare( method )
           if( exist('method','var') && ~isempty(method) )
                if( any( strcmp( this.methodOptions, method)))
                    this.method = method;
                end
           end
        end
        
        function fun = distFunHandle( this )
            switch( this.method )
                case 'euc'
                    fun = @PatchCompare.distEuc;
                case 'ncc' 
                    fun = @PatchCompare.distNCC;
                case 'sad'
                    fun = @PatchCompare.distSAD;
                otherwise
                    error( 'invalid distance method');
            end
        end
        
        function dists = pairDistances( this, X, Y ) 
            if( strcmp( this.method, 'euc'))
                dists = this.distances( X, Y );
            elseif( strcmp( this.method, 'ncc'))
                dists = this.distances( X, Y );
            else
                dists = pdist2( X, Y, this.distFunHandle() );
            end
        end
        
        function dists = distances( this, x, Y )
            
            if( strcmp( this.method, 'euc'))
                %                 fprintf('pdist2\n');
                % use a faster method for euclidean distance
                dists = pdist2( x, Y );
                dists = vecToRowCol( double(dists), 'row');
                
            elseif( strcmp( this.method, 'ncc'))
                dists = pdist2( x, Y, 'correlation' );
                dists = vecToRowCol( double(dists), 'row');
                
            else
                %                 fprintf('distance\n');
                if( numel(x) ~= size( Y, 2))
                    error( 'invalid sized for x and Y' );
                end
                
                N = size( Y, 1 );
                dists = zeros( N, 1 );
                for i = 1:N
                    dists(i) = this.distance( x, Y(i,:) );
                end
                dists = vecToRowCol( double(dists), 'row');
            end
        end
        
        function dist = distance( this, x, y )
            switch( this.method )
                case 'euc'
                    dist = PatchCompare.distEuc( x, y );
                case 'sad'
                    dist = PatchCompare.distSAD( x, y );
                case 'ncc'
                    dist = PatchCompare.distNCC( x, y );
                otherwise
                    error('invalid method');
            end
        end % distance method
    end % methods 
    
    methods( Static )
        
        function dist = distEuc( x, y )
            % Euclidean distance
            % Sqare root of sum of squared differences
            dist = norm( y - x );
        end
        
        function dist = distSAD( x, y )
            % Sum of absolute differences
            dist = sum( abs( x(:) - y(:)));
        end
        
        function dist = distNCC( x, y )
            % Normalized cross correlation
            cc = corrcoef( x, y );
            dist = 1 - cc(1,2);
        end
    end
    
end