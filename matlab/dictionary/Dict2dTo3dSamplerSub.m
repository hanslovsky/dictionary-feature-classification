classdef Dict2dTo3dSamplerSub < Dict2dTo3dSampler
    % Dict2dTo3dSamplerSub
    %
    % This class contains logic that builds a high-res 3d dictionary from 
    % 2d patches at a lower resolution. 
    %
    % John Bogovic
    % HHMI
    % February 2015
    
    properties( SetAccess = protected )
        D2d_downsampled;
    end
    
    methods
        
        function this = Dict2dTo3dSamplerSub( D2d, sz, f, overlappingPatches, scaleByOverlap )
            % Constructor
            % D2d -
            % sz  - size of 2d patches
            % f   - downsampling factor
            this = this@Dict2dTo3dSampler( D2d, sz, f, overlappingPatches, scaleByOverlap );
            this.D2d_downsampled = this.pc.downsample2dDictByDimSlow( this.D2d );
            this.pc.buildCmtx();
        end
        
        function [ dictIdxs, dictCosts, models, msk ] = bestKdicts( this, xin, i, K  )
            % returns
            
%             if( ~exist( 'isLR', 'var' ) || isempty( isLR ))
%                 isLR = true;
%             end
            
            % make sure x is a row vector
            x = vecToRowCol( xin, 'row');
            
            dists  = nan( this.numDict, 1 );
            models = {};
            
            msk = this.pc.planeMaskLRI( i );
            
%             x = x( msk > 0 );
            x = PatchConstraints.downsampleByMaskDim( x, msk );
            x = squeeze(x);
            if( size(x,1) < size(x,2) )
               x = x'; 
            end
%             size(x)
            x = [x(:)]';
            
            if( numel( x ) == size( this.D2d, 2 ) )
               D = this.D2d; 
            else
               D = this.D2d_downsampled;
            end
            
%             fprintf('length x: %d\n', length(x));
            
            if( ~isempty(this.intXfmModelType))
                models = cell( this.numDict, 1 );
                for n = 1 : this.numDict
                    bexp      = D(n,:)';
                    models{n} = fit( bexp, x, this.intXfmModelType );
                    dists(n)  = norm( x - feval( thismodel, bexp ) );
                end
            else
                dists = pdist2( x, D );
            end
            
            [ dictCosts, is ] = sort( dists );
            dictCosts = dictCosts(1:K);
            
            dictIdxs = is( 1:K );
            
            if( ~isempty(this.intXfmModelType))
                models = models( is(1:K) );
            end
        end
        
        function [ patchParams, modelList, pv, patch ] = solveHR( this, x, K )
            patchParams = zeros( this.pc.numLocs, K );
            modelList   = cell ( this.pc.numLocs, K );
            pv = zeros( prod(this.pc.sz3d), 1);
            
            patch = [];
            
            % make sure x is a column vector
            x = vecToRowCol( x, 'col');
            
            % solve for z-patches only
            zIndices = find(this.pc.dimXyzList(:,1) == 3);
            for k = 1:length( zIndices )
                
                j = zIndices( k );
                fprintf('fitting model for location %d of %d\n', ...
                    j, this.pc.numLocs );
                
                [ dictIdxs, ~, models, msk ] = this.bestKdicts( x, j, K );
                patchParams(k,:) = dictIdxs;
                
                if( ~isempty(this.intXfmModelType))
                    modelList{ j, : } = models;
                end
                
                mskHR = this.pc.planeMaskI( j );
                nnz(msk > 0)
                length( this.D2d( dictIdxs, :))
                 
                % build initial high-res patch
                pv( mskHR > 0 ) = repmat( this.D2d( dictIdxs(1), :), [1 1 this.f]);
                
            end
            
           % now do the rest
           for j = length( zIndices )+1 : d23.pc.numLocs
              
               [ dictIdxs, ~, models ] = this.bestKdicts( j, pv, 0 );
               
               patchParams(k,:) = dictIdxs;
               
               if( ~isempty(this.intXfmModelType))
                   modelList{ j, : } = models;
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
