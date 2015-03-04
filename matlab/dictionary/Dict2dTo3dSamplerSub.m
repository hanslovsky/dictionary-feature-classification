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
            
            % normalize
            this.D2d_downsampled = bsxfun(  @rdivide, ...
                                            this.D2d_downsampled, ...
                                            sqrt(sum( this.D2d_downsampled.^2, 2 ) ));
            
            this.pc.buildCmtx();
        end
        
        function [ dictIdxs, dictCosts, models, x, msk ] = bestKdicts( this, xin, i, K  )
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
            size(x)
            x = [x(:)]';
            
            % normlize x
            x = x ./ norm( x );
            
            if( numel( x ) == size( this.D2d, 2 ) )
               D = this.D2d; 
            else
               D = this.D2d_downsampled;
            end
            
            % check size of K
            K = min( K, size( D, 1 ));
            
            fprintf('length x: %d\n', length(x));
            fprintf('mag x: %d\n',    norm(x));

            % this.intXfmModelType
            
            if( ~isempty(this.intXfmModelType))
                fprintf('fitting model: %s\n', this.intXfmModelType);
                models = cell( this.numDict, 1 );
                for n = 1 : this.numDict
                    bexp      = D(n,:);
                    models{n} = fit( bexp', x', this.intXfmModelType, this.fitParams{:} );
                    dists(n)  = norm( x' - feval( models{n}, bexp' ) );
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
        
        function [ patchParams, modelList, pv, patch ] = solveBestK( this, x, K, cutoffFun )
            if( ~exist( 'cutoffFun', 'var'))
                cutoffFun = [];
            end
            
            % patchParams = zeros( this.pc.numLocs, K );
            patchParams = cell( this.pc.numLocs, 1 );
            
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
                
                [ dictIdxs, dictCosts, models, msk ] = this.bestKdicts( x, j, K );
                if( ~isempty( cutoffFun ))
                    Kfun = cutoffFun( dictCosts, K );
                    patchParams{k} = dictIdxs( 1:Kfun );
                else
                    patchParams{k} = dictIdxs;
                end
                
                if( ~isempty(this.intXfmModelType))
                    modelList( j, : ) = models;
                end
                
                mskHR = this.pc.planeMaskI( j );
                %nnz(msk > 0)
                length( this.D2d( dictIdxs, :))
                
                % build initial high-res patch
                % pv( mskHR > 0 ) = repmat( this.D2d( dictIdxs(1), :), [1 1 this.f]);
                
            end
            
            
        end
        
        function [ patchParams, modelList, pv, patch ] = solveBestK_maybeBad( this, x, K, cutoffFun )
            % Solve for patch parameters from an observation x by:
            %   1)  solving for 2d patches that fit z-plane sub patches
            %   2)  estimate a high-res patch from those 
            %   3)  fit the remaining patches to that HR patch
            %
            %   This seems like a bad idea, actually -JB (2015 Mar 02)
            %   Probably preferable to fit all patches directly to the
            %       observation
            
            if( ~exist( 'cutoffFun', 'var'))
                cutoffFun = [];
            end
            
            % patchParams = zeros( this.pc.numLocs, K );
            patchParams = cell( this.pc.numLocs, 1 );
            
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
                
                [ dictIdxs, dictCosts, models, msk ] = this.bestKdicts( x, j, K );
                if( ~isempty( cutoffFun ))
                    Kfun = cutoffFun( dictCosts, K );
                    patchParams{k} = dictIdxs( 1:Kfun );
                else
                    patchParams{k} = dictIdxs;
                end
                
                if( ~isempty(this.intXfmModelType))
                    modelList( j, : ) = models;
                end
                
                mskHR = this.pc.planeMaskI( j );
                %nnz(msk > 0)
                length( this.D2d( dictIdxs, :))
                
                % build initial high-res patch
                % pv( mskHR > 0 ) = repmat( this.D2d( dictIdxs(1), :), [1 1 this.f]);
                
            end
            
            % now do the rest
            for j = length( zIndices )+1 : this.pc.numLocs
                
                fprintf('fitting model for location %d of %d\n', ...
                    j, this.pc.numLocs );
                
                % DO NOT FIT TO ESTIMATED HIGH-RES PATCH! 
                % [ dictIdxs, distList, models ] = this.fitIdxAndModel( j, pv, 0, 1 );
                
                % FIT TO OBSERVATION
                [ dictIdxs, distList, models ] = this.fitIdxAndModel( j, pv, 0, 1 );
                
                if( ~isempty( cutoffFun ))
                    Kfun = cutoffFun( distList, K );
                    patchParams{ j } = dictIdxs( 1:Kfun );
                else
                    % patchParams( j, : ) = dictIdxs(1:K);
                    patchParams{ j } = dictIdxs(1:K);
                end
                
                if( ~isempty(this.intXfmModelType))
                    modelList( j, : ) = models( 1:K );
                end
            end
            
            if( nargout > 3 )
                patch = reshape( pv, this.sz3d );
            end
            
        end
        
        function [ bestPatchParams, currentModels, pvOut ] = greedySubSearch( this, x, iniPatchParams, iniModels,...
                maxOuterIters, maxNoChangeIters )
            
            if( ~exist('maxOuterIters','var') || isempty(maxOuterIters))
                maxOuterIters = 200;
            end
            if( ~exist('maxNoChangeIters','var') || isempty(maxNoChangeIters))
                maxNoChangeIters = 20;
            end
            
            [ numLocs, K ] = size( iniPatchParams );
            
            if( numLocs ~= this.pc.numLocs )
                error('number of rows of iniPatchParams must equal this.pc.numLocs');
            end
            
%             currentParams = iniPatchParams( :, 1 );
            currentParams = cellfun( @(x)(x(1)), iniPatchParams );
            currentModels = iniModels( :, 1 );
            
            [ pv ] = this.patchFromParams( currentParams, currentModels );
            currentCost = norm( pv - x );
            
            numUnchanged = 0;
            for outeriter = 1:maxOuterIters
                
                fprintf('iteration %d of %d\n', outeriter, maxOuterIters );
                
                % shuffle order to optimize locations
                innerList = randperm( numLocs );
                
                % loop over locations
                changed = false;
                for i = 1:numLocs
                    j = innerList( i );
                    
                    % it is possible there is only one possible dictionary
                    % element for location j.  if so, ski;
                    if( length( iniPatchParams{j}) == 1 )
                        continue;
                    end
                    
                    % loop over the K possible patches at location j
                    if( iscell( currentParams) )
                        bestPatchIdx = currentParams{ j };
                    else
                        bestPatchIdx = currentParams( j );
                    end
                    
                    tmpParams = currentParams;
                    tmpModels = currentModels;
                    
                    changedInner = false;
                    for k = setdiff( 1:K, find( iniPatchParams{j} == bestPatchIdx) )
                        
                        currentIdx = iniPatchParams( j, k );
                        % diff for debug only
%                         diff = zeros( size( tmpParams ));
%                         diff( j ) = 1;
                        
                        tmpParams( j ) = currentIdx;
                        tmpModels{ j } = iniModels{ j, k };
%                         [ tmpParams diff ]
%                         pause;
                        
                        %thismodel = fit( bexp, AxR, this.intXfmModelType );
                        
                        [ pv_tmp ] = this.patchFromParams( tmpParams, tmpModels );
                        cost = norm( pv_tmp - x );
                            
                        if( cost < currentCost )

%                             fprintf( 'old %f vs new %f costs at location %d\n', ...
%                                             currentCost, cost );
                            if( cost < currentCost )
                                fprintf( 'cost improved from %f to %f at location %d\n', ...
                                            currentCost, cost );
                                changedInner = true;
                                currentCost = cost;
                            end
                            
                        end
                    end % parameters
                    
                    if( changedInner )
                        currentParams = tmpParams;
                        currentModels = tmpModels;
                        changed = true;
                    end
                    
                end % inner loop
                
                % increment changed counter
                if( changed )
                    numUnchanged = 0;
                else
                    numUnchanged = numUnchanged + 1;
                end
                
                % nothing has changed for awhile, so break early
                if( numUnchanged > maxNoChangeIters )
                    fprintf( 'no change, breaking early\n');
                    break;
                end
                
            end % outer loop
            
            bestPatchParams = currentParams;
            pvOut = this.patchFromParams( bestPatchParams, currentModels );
            
        end % greedySubSearch
       
        
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
    
    methods( Static )
        
        function [ K ] = chooseCutoffInflection( dictCosts, L ) 
            deriv = diff(dictCosts(1:L));
            deriv2 = diff( deriv );
%             K = find( (deriv(1:end-1) > 0) & (deriv2 < 0) , 1 );
            K = find( (deriv2 < 0) , 1 );
        end 
    end
end
