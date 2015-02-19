classdef Dict2dTo3dSamplerSub < Dict2dTo3dSampler
    % Dict2dTo3dSamplerSub
    %
    % This class contains logic that builds a high-res 3d dictionary from 
    % 2d patches at a lower resolution. 
    %
    % John Bogovic
    % HHMI
    % February 2015
    
    methods
        
        function this = Dict2dTo3dSamplerSub( D2d, sz, f, overlappingPatches, scaleByOverlap )
            % Constructor
            % D2d -
            % sz  - size of 2d patches
            % f   - downsampling factor
            this = this@Dict2dTo3dSampler( D2d, sz, f, overlappingPatches, scaleByOverlap );
        end
        
        function [ dictIdxs, dictCosts, models ] = bestKdicts( this, x, i, K )
            % returns
            
            % make sure x is a row vector
            x = vecToRowCol( x, 'row');
            
            dists  = nan( this.numDict, 1 );
            models = {};
            
            msk = this.pc.planeMaskLRI( i );
            x = x( msk > 0 );
            
            if( ~isempty(this.intXfmModelType))
                models = cell( this.numDict, 1 );
                for n = 1 : this.numDict
                    bexp      = this.D2d(n,:)';
                    models{n} = fit( bexp, x, this.intXfmModelType );
                    dists(n)  = norm( x - feval( thismodel, bexp ) );
                end
            else
                dists = pdist2( x, this.D2d );
            end
            
            [ dictCosts, is ] = sort( dists );
            dictCosts = dictCosts(1:K);
            
            dictIdxs = is( 1:K );
            
            if( ~isempty(this.intXfmModelType))
                models = models( is(1:K) );
            end
        end
        
        function [ patchParams, modelList, pv, patch ] = solveHR( this, x )
            patchParams = zeros( this.pc.numLocs, 1 );
            modelList   = cell ( this.pc.numLocs, 1 );
            pv = [];
            patch = [];
            
            % make sure x is a column vector
            x = vecToRowCol( x, 'col');
            
            % solve for z-patches only
            zIndices = find(this.pc.dimXyzList(:,1) == 3);
            for k = 1:length( zIndices )
                
                j = zIndices( k );
                fprintf('fitting model for location %d of %d\n', ...
                    j, this.pc.numLocs );
                
                [ idx, ~, model ] = this.fitIdxAndModel( j, x );
                patchParams(j) = idx;
                
                if( ~isempty(this.intXfmModelType))
                    modelList{ j } = model;
                end
            end
        end
        
    end  % methods
    
    methods ( Access = protected )
        
        function iniPatchConstraints( this )
            % initialized the PatchConstraintsSubdim object.
            %   it allows patches to overlap in x and y but not in z
            %   (the assumption being that observed patches are z-projections
            %   of high-res patches)
            this.pc = PatchConstraintsSubdim( this.sz3d(1:2), this.sz3d, this.f, ...
                this.overlappingPatches, this.scaleByOverlap );
        end
    end % protected methods
    
end
