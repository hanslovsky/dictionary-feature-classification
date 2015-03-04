classdef PatchCompare < handle
   
    properties
        methodOptions = { 'euc', 'sad', 'ncc' };
        % euc - euclidean distance 
        % sad - sum absolute distances
        % ncc - normalized cross correlation
        
        method = 'euc';  % euclidean norm
    end
    
    methods
        
        function this = PatchCompare( method )
           if( exist('method','var') && ~isempty(method) )
                if( any( strcmp( methodOptions, method)))
                    this.method = method;
                end
           end
        end
        
        function dist = distance( this, x, y )
            switch( this.method )
                case 'euc'
                    dist = distEuc( x, y );
                case 'sad'
                    dist = distSad( x, y );
                case 'ncc'
                    dist = distNcc( x, y );
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
            dist = cc(1,2);
        end
    end
    
end