classdef Dict2dTo3dSampler < Dict2dTo3d
    % Dict2dTo3dSampler
    %
    % This class contains logic that builds a high-res 3d dictionary from 
    % 2d patches at a lower resolution using a sampling technique rather
    % than a Astar-like search.
    %
    % John Bogovic
    % HHMI
    % September 2014
    
    properties
       maxIters  = 100; % max iterations for build
       convIters =  20; % # of iters at a given cost ('flatness')
                        % that indicate convergence
       convEps =  0.001; % convergence epsilon (defines 'flatness')
       
       chooseBestAtIter = 0;
       useSubset = 0; 
       
       intXfmModelType = '';
       paramModels;
       
       scaleDictElems = 0;
       paramScales;
       
       recordParamsOverIters = 0;
       
    end
    
    properties( SetAccess = protected )
       
    end
    
    methods
        
        % Constructor
        % D2d - 
        % sz  - size of 2d patches
        % f   - downsampling factor
        function this = Dict2dTo3dSampler( D2d, sz, f, overlappingPatches )
            this = this@Dict2dTo3d( D2d, sz, f, overlappingPatches );
            if( this.useSubset )
                fprintf('Computing insersection matrix inverses...');
                this.pc.compXsectInverses();
                fprintf('.done\n');
            end
        end
        
        function [ patchParams, iteration, costs ] = build3dPatchIni_old( this, iniPatch, locs )
        % [ patchParams, iteration, costs ] = build3dPatchIni( this, iniPatch, locs )
        % 
        % Inputs
        %   iniPatch - the patch for initialization 
        %   locs     - the locations in the 3d patch
        %              can be specified either as a [Nx2] matrix [ dim, xyz]
        %              or a [Nx1] vector of indices into this.pc.dimXyzList
            
            if( ~exist('locs','var') || isempty( locs ))
                locs = this.iniLocs;
            end
            
            % check format of locs
            if( size( locs, 2) == 2 )
                locidxs = this.pc.locXyzDim2Idx( locs(:,1), locs(:,2));
            elseif( size( locs, 2) == 1 )
                locidxs = locs;
            else
                error( 'invalid format for locs');
            end
            
            numTmpPatches = this.addTemporaryPatchesToDict( iniPatch, 0 );
            
            % build an initialization from the above
            iniParams = zeros( this.pc.numLocs, 1 );
            iniParams( locidxs ) = (this.numDict - numTmpPatches + 1) : this.numDict;
            iniParams( setdiff( 1:this.pc.numLocs, locidxs )) = ...
                    randi( (this.numDict - numTmpPatches), ...
                           (this.pc.numLocs - length(locidxs)), 1);
            
            exclude = false( this.pc.numLocs, 1);
            exclude( locidxs ) = true;
            
            [ patchParams, iteration, costs ] = build3dPatch( this, iniParams, exclude );
            
            this.removeTemporaryPatches( numTmpPatches );
        end
        
        function [ patchParams, iteration, costs ] = build3dPatchIni( this, iniPatch, locs )
        % [ patchParams, iteration, costs ] = build3dPatchIni( this, iniPatch, locs )
        % 
        % Inputs
        %   iniPatch - the patch for initialization 
        %   locs     - the locations in the 3d patch
        %              can be specified either as a [Nx2] matrix [ dim, xyz]
        %              or a [Nx1] vector of indices into this.pc.dimXyzList
            
            if( ~exist('locs','var') || isempty( locs ))
                locs = this.iniLocs;
            end
            
            % check format of locs
            if( size( locs, 2) == 2 )
                locidxs = this.pc.locXyzDim2Idx( locs(:,1), locs(:,2));
            elseif( size( locs, 2) == 1 )
                locidxs = locs;
            else
                error( 'invalid format for locs');
            end
            
            % build an initialization from the above
            iniParams = zeros( this.pc.numLocs, 1 );
            iniParams( locidxs ) = (this.numDict - numTmpPatches + 1) : this.numDict;
            iniParams( setdiff( 1:this.pc.numLocs, locidxs )) = ...
                    randi( (this.numDict - numTmpPatches), ...
                           (this.pc.numLocs - length(locidxs)), 1);
            
            exclude = false( this.pc.numLocs, 1);
            exclude( locidxs ) = true;
            
            [ patchParams, iteration, costs ] = build3dPatch( this, iniParams, exclude );
            
        end 
        
        function [ patchParams, pv, iteration, costs ] = build3dPatch( this, iniPatch, excludeParam )
        % [ patchParams, iteration, costs ] = build3dPatch( this, iniPatch, exclude )
        %
        % Inputs:
        %   this    : this object
        %   iniPatch: initialization parameters
        %   exclude : indices that should not be optimized over
        %             useful if initializing with observations or another
        %             dictionary
        
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            N = this.pc.numLocs;
            isExclusion = 0;
            if( ~exist( 'excludeParam', 'var') || isempty( excludeParam ))
                exclude = false( N, 1 );
            else
                isExclusion = 1;
                if( ~islogical( excludeParam ))
                    exclude = false( this.pc.numLocs, 1);
                    exclude( excludeParam  ) = true;
                else
                    exclude = excludeParam;
                end
                N = nnz( ~exclude ); 
            end
            
            % initialize patch parameters
            % ie which 2d patch goes into which 3d location
            if( exist( 'iniPatch', 'var' ) && ~isempty( iniPatch ))
                patchParams = iniPatch;
            else
                patchParams = randi( this.numDict, N, 1 );
            end
            
             if( this.scaleDictElems )
                 this.paramScales = ones( this.pc.numLocs, 1 );
             end
             
             % initialize models if we're using them
             if( ~isempty( this.intXfmModelType ))
                 this.initializeModels( patchParams );
             end
            
            converged = 0;           
            iteration = 1;
            
%             % randomly generate the patch location 
%             % that will be updated at each iteration
%             randomCoords = randi( N, this.maxIters, 1 );
%             if( isExclusion )
%                 keepInds = 1:this.pc.numLocs;
%                 keepInds = keepInds( ~exclude );
%                 randomCoords = keepInds( randomCoords );
%             end
            randomCoords = this.genRandomCoords( N, exclude );
            
            b = this.pc.constraintValueList( this.D2d, patchParams );
            costs = -1.*ones( this.maxIters, 1 );
            
            lastCost = inf;
            itersAtThisCost = 0;
            
            rcDups = 0;
            while( ~converged )
                fprintf(' iteration: %d\n', iteration );
                if( (iteration - rcDups) > length(randomCoords) )
                    randomCoords = genRandomCoords( this, N, exclude );
                    rcDups = rcDups + this.maxIters;
                end
                i = randomCoords( iteration - rcDups );
                
                if( ~isempty( this.intXfmModelType ) ) 
                    [ bestidx, theseCosts, model ] = ...
                        this.goodPatchConfigModel( b, i, ...
                        this.chooseBestAtIter );
                                          
                    this.paramModels{i} = model;
                    
                elseif( this.scaleDictElems )
                    [ bestidx, scaleCost, scale ] = ...
                        this.goodPatchConfigScale( b, i, ...
                                              this.chooseBestAtIter );
                    
                    this.paramScales(i) = scale;
                else
                    if( this.useSubset )
                        [bestidx, theseCosts] = this.bestPatchConfigSub( b, i );
                    else
                        if( this.chooseBestAtIter )
                            [bestidx, theseCosts] = this.bestPatchConfig( b, i );
                        else
                            [bestidx, theseCosts] = this.goodPatchConfig( b, i );
                        end
                    end
                end
                
                if( isempty( bestidx ))
                    itersAtThisCost = itersAtThisCost + 1;
                    iteration = iteration + 1;
                    
                else
                    % were picking the best patch and updating it
                    patchParams( i ) = bestidx;
                    theseCosts = this.totalPatchCost( patchParams );
                    
                    if( numel( theseCosts ) == 1 )
                        costs(iteration) = theseCosts;
                    else
                        costs(iteration) = theseCosts( bestidx );
                    end
                    
                    % has the cost changed much?
                    if( abs( costs(iteration) - lastCost  ) < this.convEps )
                        itersAtThisCost = itersAtThisCost + 1;
                    else
                        itersAtThisCost = 0;
                    end
                    lastCost = costs(iteration);
                    
                    b = this.pc.updateConstraints( this.D2d, b, i, bestidx );
                end
                
                % converged if we've been at the same cost for awhile
                if( itersAtThisCost == this.convIters )
                    converged = 1;
                end
                
                % force exit after max iteration count
                if(iteration >= this.maxIters)
                    converged = 1;
                else
                    iteration = iteration + 1;
                end
            end % iteration
            
            if( this.verbose )
                fprintf('sampler - buildPatch3d converged after %d iterations\n', (iteration-1) );
            end
            
            costs = costs(1:iteration-1);
            patchParams = vecToRowCol( patchParams, 'col');
            
            patchParams = [ this.pc.dimXyzList, ...
                            patchParams ];
            
            [pv] = this.patchFromParams( patchParams );
            
        end % build3dPatch
        
        function cost = totalPatchCost( this, params )
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            b = this.constraintValue( params );
            x = cmtxInv * b;
            cost = norm( cmtx * x - b );
        end
         
        function randomCoords = genRandomCoords( this, N, exclude )
         % randomly generate the patch location 
         % that will be updated at each iteration
            randomCoords = randi( N, this.maxIters, 1 );
            if( exist('exclude','var') && ~isempty( exclude ))
                keepInds = 1:this.pc.numLocs;
                keepInds = keepInds( ~exclude );
                randomCoords = keepInds( randomCoords );
            end
        end
        
        function iniParams = iniParamsDist( this, N, method )
            switch( method)
                case 'rand'
                    M = this.pc.numLocs;
                    iniParams = mat2cell( ...
                                    randi(  this.numDict, N, M), ...
                                ones( N, 1 ));
                    
                case 'dict'
                    
                    if( isempty( this.Dini ))
                       error('Dini is empty!!'); 
                    end
                    
                    numIni = size( this.Dini, 1 );
                    numPatchesPerIni = length( this.iniLocs );
                     
                    iniIdxHelper = reshape( 1:(numIni*numPatchesPerIni), numPatchesPerIni, [] )';
                    
                    M = this.pc.numLocs;
                    
                    numPatchesPerIni = length( this.iniLocs );
            
                    Dinirs = reshape( this.Dini', ...
                                        [], numIni*numPatchesPerIni)';
                    
                    this.numIniPatches = this.addTemporaryPatchesToDict( Dinirs, 0 );
                    
                    iniParams = zeros( N, M );
                    optLocSet = setdiff( 1:M, this.iniLocs );
                    numOptLoc = length( optLocSet  );
                    iniParams( :, setdiff( 1:M, this.iniLocs )) = ...
                        randi(  this.numDict, N, numOptLoc );
                    
                    iniParams( :, this.iniLocs ) = ...
                        this.numDict + iniIdxHelper( randi( numIni, N, 1), : );
                    
                    iniParams = mat2cell( iniParams, ...
                                          ones( N, 1 ));
                    
                    
                otherwise
                    error('invalid ini method');
                    
            end
        end
        
        function iniParams = iniParamsDictOld()
            % size parameters
            M = size( this.Dini, 1 );
            numPatchesPerIni = length( this.iniLocs );
            
            chosenPatches = this.Dini( randi( M, N, 1), :);
            
            iniParams = mat2cell( ...
                            reshape( chosenPatches', ...
                                     [], N*numPatchesPerIni)', ...
                            length(this.iniLocs)*ones( N, 1) );
        end
        
        function [ idx, dist, cmtx ] = goodPatchConfig( this, b, i )
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            
            % initialize experimental constraint values
            bexp = b;
            
            x = cmtxInv * bexp;
            curdist = norm( cmtx * x - bexp );
            
            % the range of constraint values that will change
            % depending on the patch being tested
            brng = this.pc.constraintVecSubsets(i,:);
            
            idx = [];
            testOrder = randperm( this.numDict );
            for nn = 1:this.numDict
                n = testOrder( nn );
                bexp( brng ) = this.D2d(n,:);
                
                x = cmtxInv * bexp;
                dist = norm( cmtx * x - bexp );
                
                if( dist < curdist )
                   idx = n;
                   return;
                end
            end
            
        end
        
        function initializeModels( this, iniParams, convEps, maxIniIters ) 

            if( ~exist( 'convEps', 'var' ) || isempty( convEps ))
                convEps = 0.001;
            end
            if( ~exist( 'maxIniIters', 'var' ) || isempty( maxIniIters ))
                maxIniIters = 5;
            end
            
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            
            this.paramModels = repmat( ...
                                {cfit(fittype('poly1'), 1, 0 )}, ...
                                this.pc.numLocs, 1 );
            
            i = 1;
            converged = 0;
            lastx = [];
            diff = 1;
            
            while( ~converged )
                b = this.constraintValue( iniParams );
                
                x  = cmtxInv * b;
                x = x ./ norm( x );
                Ax = cmtx * x;
                
                for l = 1:this.pc.numLocs
                    
                    rng = this.pc.constraintVecSubsets(l,:);
                    AxR = Ax( rng );
                    
                    if( var( AxR ) < 0.0001 )
                        AxR = AxR + 0.0001.*randn(size(AxR));
                    end
                    if( var( b(rng) ) < 0.0001 )
                        b(rng) = b(rng) + 0.0001.*randn(size(b(rng)));
                    end
                    
                    thismodel = fit( b(rng), AxR, this.intXfmModelType );
                   
                    this.paramModels{l} = thismodel;
                    
%                     plot( AxR(:),  b(rng), '*' );
%                     thismodel
%                     hold on; plot( AxR(:),  feval( thismodel, b(rng)), '*r' );
%                     pause;
%                     close all;
                end
                
                if( ~isempty(lastx) )
                    diff =  norm( x(:) - lastx(:));
                end
                if( diff < convEps )
                    converged = 1;
                end
                
                lastx = x;
                
                i = i + 1;
                if( i > maxIniIters )
                    break;
                end
            end % while
        end
        
        function [ idx, dist, model ] = goodPatchConfigModel( this, b, i, doBest )
            
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            
            % initialize experimental constraint values
            bexp = b;
            
            % the range of constraint values that will change
            % depending on the patch being tested
            rng = this.pc.constraintVecSubsets(i,:);
            
            x = cmtxInv * bexp;
            
            Ax  = cmtx * x;
            AxR = Ax( rng );
            if( var( AxR ) < 0.0001 )
                AxR = AxR + 0.0001.*randn(size(AxR));
            end
                    
            curdist = norm( AxR - feval( this.paramModels{i}, bexp( rng )) );
            
            idx = [];
            model = this.paramModels{i}; % the current model
            testOrder = randperm( this.numDict );
            for nn = 1:this.numDict
                n = testOrder( nn );
                
                bexp( rng ) = this.D2d(n,:);
                if( var( bexp(rng) ) < 0.0001 )
                        bexp(rng) = bexp(rng) + 0.0001.*randn(size(bexp(rng)));
                end
                    
                    
                thismodel = fit( bexp(rng), AxR, this.intXfmModelType );
                dist = norm( AxR - feval( thismodel, bexp( rng )) );
                
                if( (dist < curdist) )
                    idx   = n;
                    model = thismodel;
                    curdist = dist;
                    if( ~doBest )
                        return;
                    end
                end
                
            end % loop
        end % goodPatchConfigModel
        
        function [ idx, dist, scale ] = goodPatchConfigScale( this, b, i, doBest )
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            
            % initialize experimental constraint values
            bexp = b;
            
            % the range of constraint values that will change
            % depending on the patch being tested
            rng = this.pc.constraintVecSubsets(i,:);
            
            x = cmtxInv * bexp;
            %curdist = norm( cmtx * x - bexp );
            
            Ax  = cmtx * x;
            AxR = Ax( rng );
%             AxRscale = norm( AxR );
            
            scale = this.paramScales( i );
            curdist = norm( AxR - scale.*bexp( rng ) );
            
            idx = [];
            testOrder = randperm( this.numDict );
%             testOrder = 1:this.numDict
            for nn = 1:this.numDict
                n = testOrder( nn );
                
                bexp( rng ) = this.D2d(n,:);
                
%                 cBscale  = norm( bexp( rng ));
                
                scaleTest = Dict2dTo3dSampler.scalingFactor( AxR, bexp(rng ), 'norm');
                dist = norm( AxR - scaleTest.*bexp( rng ) );
                
%                 if( n == 1 )
%                    aaaa = 1; 
%                 end
                if( (dist < curdist) )
                    idx   = n;
                    scale = scaleTest;
                    curdist = dist;
                    if( ~doBest )
                        return;
                    end
                end
                
            end
            
        end
        
        function [ bestidx, sims, cmtx ] = bestPatchConfig( this, b, i )
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            
            % initialize experimental constraint values
            bexp = b;
            
            % the range of constraint values that will change
            % depending on the patch being tested
            brng = this.pc.constraintVecSubsets(i,:);
            
            sims = zeros( this.numDict, 1 );
            for n = 1:this.numDict
                bexp( brng ) = this.D2d(n,:);
                
                x = cmtxInv * bexp;
                sims( n ) = norm( cmtx * x - bexp );
            end
            
            [ ~, bestidx ] = min( sims );
        end
        
        function [ bestidx, sims, cmtx, btot ] = bestPatchConfigSub( this, b, i )
        %  [bestidx, sims] = bestPatchConfig( this, b, i )
        %   b - vector of all constraint values
        
            cmtx  = this.pc.subCmtxAndInvs{i,1};
            cmtxi = this.pc.subCmtxAndInvs{i,2};
            
            brng = this.pc.constraintVecXsectSubsets{i};
            bsub = b( brng );
            
            sims = zeros( this.numDict, 1 );
            for n = 1:this.numDict
               
                % change bsub depending on which patch is being tested
                btot = [ this.D2d(n,:)'; ... 
                         bsub ];
                 
                x = cmtxi * btot;
                sims( n ) = norm( cmtx * x - btot );
                
            end
            [ ~, bestidx ] = min( sims );
        end
    
        function b = constraintValue( this, obj )
            if( isa( obj, 'net.imglib2.algorithms.opt.astar.SortedTreeNode'))
                b = this.pc.constraintValueNode( this.D2d, obj );
            else
                b = this.pc.constraintValueList( this.D2d, obj );
            end
            
            if( ~isempty( this.intXfmModelType ))
                for i = 1:this.pc.numLocs
                    rng = this.pc.constraintVecSubsets(i,:); 
                    b( rng ) = feval(this.paramModels{i}, b( rng ));
                end
            elseif( this.scaleDictElems )
                for i = 1:this.pc.numLocs
                    rng = this.pc.constraintVecSubsets(i,:); 
                    b( rng ) = this.paramScales(i) .* b( rng );
                end
            end
        end
    end
    
    methods( Static )
        
        function s = scalingFactor( truevec, testvec, type )
            switch( type )
                case 'norm'
                    s = norm( truevec ) ./ norm( testvec );
                case 'sum'
                    s = sum( truevec ) ./ sum( testvec );
                case 'mean'
                    s = mean( truevec ) ./ mean( testvec );
                otherwise
                    error('invalid scaling type');
                    
            end
        end
        
        function [err] = patchConfigConsistency( cmtx, x, b, type )
            switch( type )
                case 'ssd'
                    err = norm( cmtx * x - b );
                case 'scaled'
                    
                otherwise
                    error('invalid type parameter');
            end
            
        end
    end
end