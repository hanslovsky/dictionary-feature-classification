classdef Drr < handle
    % Drr - "dictionary redundancy remover"
    
    properties ( SetAccess = private )
        tid; % the Tid object
        SimMtx; % the pair-wise similarity matrix
    end
    
    properties
       method = 'degree'; 
    end
    
    methods
        
        function this = Drr( tid )
            this.tid = tid;
        end
        
        function setSimMtx( this, simmtx )
            this.SimMtx = simmtx;
        end
        
        function ri = prune( this, simmtx )
            
            switch ( this.method )
                case 'degree'
                    ri = Drr.pruneDegree( simmtx );
                otherwise
                    error( 'invalid method' );
            end
        end
        
        
    end % methods
    
    methods( Static )
        
        function redundantIndices = pruneDegree( simmtx )
            
            N = size( simmtx, 1 );
            redundantIndices = false( N , 1 );
            
            degrees = sum( simmtx );
            
            % [ ds ] = sort( degrees, 'descend' );
            ds = unique( degrees );
            ds = ds( end: -1 : 1 );
            
            % process nodes in order of decreasing degrees
            for i = 1:length( ds )
                
                %fprintf('%d th degree is %d\n', i, ds(i));
                subIall = (degrees == ds( i ));
                
                if( nnz( subIall ) > 1 )
                    subSimMtx = simmtx( subIall, subIall );
                    M = size( subSimMtx, 1 );
                    
                    rp = find( subIall );
                    rp = rp( randperm( M ));
                else
                    rp = find( subIall );
                end
                
                for jj = 1:length(rp)
                    
                    j = rp(jj);
                    %fprintf('inner iter %d, touching %d\n', jj, j );
                    
                    % make sure we don't bother with
                    % nodes that have already been marked as
                    % redundant 
                    if( redundantIndices(j) ) %|| subIall(j)
                        continue;
                    end
                    
                    simrow = simmtx( j, : );
                    
                    % mark the neighbors of the selected node
                    % as redundant
                    redundantIndices( (simrow > 0) ) = true;
                    
                    
                end % loop over nodes with this degree
            end % loop over degrees
        end % prune
        
    end % static methods
end