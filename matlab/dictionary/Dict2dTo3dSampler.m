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
       alwaysNormD3d    = 0;
       useSubset = 0; 
       
       intXfmModelType = '';
       stopAfterFitParamIni = 1;
       fitParams = { 'Robust', 'LAR' };
       
       scaleDictElems = 0;
       paramScales;
       
       recordParamsOverIters = 0;
       savedParams;

%        updateOrder = 'all';
       updateOrder = 'each';
       verboseEvery=100;
       
    end
    
    properties( SetAccess = protected )
       comparator;
       D2d_downsampled;
    end
    
    methods
        
        
        function this = Dict2dTo3dSampler( D2d, sz, f, overlappingPatches, scaleByOverlap, comparator )
            % Constructor
            % D2d - 
            % sz  - size of 2d patches
            % f   - downsampling factor
            
            this = this@Dict2dTo3d( D2d, sz, f, overlappingPatches, scaleByOverlap );
            if( this.useSubset )
                fprintf('Computing insersection matrix inverses...');
                this.pc.compXsectInverses();
                fprintf('.done\n');
            end
            
            if( this.scaleByOverlap )
                this.paramScales = this.pc.overlapFraction;
            end
            if( ~exist( 'comparator','var'))
                comparator = PatchCompare( 'euc' );
            end
            this.comparator = comparator;
            
            % can't use isa here, because isa returns true 
            % if the argument is a subclass, which is not desired
            if( strcmp( class( this ), 'Dict2dTo3dSampler')) %#ok<STISA>
                this.D2d_downsampled = this.pc.downsample2dDictByDimSlow( this.D2d );
            end
            
            % normalize
            this.D2d_downsampled = bsxfun(  @rdivide, ...
                this.D2d_downsampled, ...
                sqrt(sum( this.D2d_downsampled.^2, 2 ) ));
                
        end
        
        function obj = copy(this)
            obj = Dict2dTo3dSampler( this.D2d, this.sz2d(1), this.f, ...
                                     this.overlappingPatches, this.scaleByOverlap );
            obj.sz3d = this.sz3d;
            
            obj.D3d = this.D3d;
            obj.numDict3d = this.numDict3d;
            
            obj.clone( this );
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
        
        function [ patchParams, pv, iteration, costs, patchesByIter, otherParams ] = ...
                        build3dPatch( this, iniPatch, excludeParam, record )
        % [ patchParams, iteration, costs ] = build3dPatch( this, iniPatch, excludeParam, record )
        %
        % Inputs:
        %   this    : this object
        %   iniPatch: initialization parameters
        %   exclude : indices that should not be optimized over
        %             useful if initializing with observations or another
        %             dictionary
        %   record  : if true, the fifth output will be non-empty and
        %             contain the estimated 3d patch at every iteration 
        
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
%             if( exist( 'objectParams', 'var') && ~isempty( objParams ))
%                 fprintf('setting object parameters\n');
%                 setObjectProperties( this, objParams{:} );
%             end
            
            if( ~exist( 'record', 'var') || isempty( record ))
                record = 0; 
            end
            patchesByIter = {};
            otherParams   = {};
            
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
            modelsSet = 0;
            if( exist( 'iniPatch', 'var' ) && ~isempty( iniPatch ))
                
                if( ~isempty( this.intXfmModelType ))
                    
                    % TODO FIX ME FIX ME!! - WHAT A MESS! :(
                    if( ~isempty( this.Dini ))
                        this.paramModels = cell( this.pc.numLocs, 1);
                        patchParams = zeros( this.pc.numLocs, 1 );
                    
                        x = this.Dini( iniPatch, : );
                        
                        for j = 1:this.pc.numLocs
                            fprintf('fitting model for location %d of %d\n', ...
                                j, this.pc.numLocs );
                            
                            [ idx, ~, model ] = this.fitIdxAndModel( j, x' );
                            patchParams(j) = idx;
                            this.paramModels{ j } = model;
                        end
                        modelsSet = 1;
                        
                        if( this.stopAfterFitParamIni )
                            % stop early and don't iterate
                            [pv] = this.patchFromParams( patchParams );
                            costs = this.totalPatchCost();
                            return;
                        end
                    else
                        this.paramModels = cell( this.pc.numLocs, 1);
                        patchParams = iniPatch;
                    end
                else
                    patchParams = iniPatch;
                end
                
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
            extra = 10;
            randomCoords = this.genRandomCoords( N, exclude, extra );
            
            b = this.pc.constraintValueList( this.D2d, patchParams );
            costs = -1.*ones( this.maxIters, 1 );
            
            lastCost = inf;
            itersAtThisCost = 0;
            
            rcDups = 0;
            while( ~converged )
                
                fprintf(' iteration: %d\n', iteration );
                if( (iteration - rcDups) > this.maxIters )
                    randomCoords = genRandomCoords( this, N, exclude, extra );
                    rcDups = rcDups + this.maxIters;
                end
                
                if( strcmp( this.updateOrder, 'each' ))
                    i = randomCoords( iteration - rcDups );
                elseif( strcmp( this.updateOrder, 'all' ))
                    i = randomCoords( iteration - rcDups, : );
                end
                
                if( ~isempty( this.intXfmModelType ) && ~modelsSet ) 
                    [ bestidx, ~, model ] = ...
                        this.goodPatchConfigModel( b, i, ...
                        this.chooseBestAtIter );

                elseif( this.scaleDictElems )
                    [ bestidx, ~, scale ] = ...
                        this.goodPatchConfigScale( b, i, ...
                                              this.chooseBestAtIter );
                    
                    this.paramScales(i) = scale;
                else
                    if( this.useSubset )
                        [bestidx, ~] = this.bestPatchConfigSub( b, i );
                    else
                        if( this.chooseBestAtIter )
                            [bestidx, ~] = this.bestPatchConfig( b, i );
                        else
                            [bestidx, ~] = this.goodPatchConfig( b, i );
                        end
                    end
                end
                
                if( isempty( bestidx ))
                    itersAtThisCost = itersAtThisCost + 1;
                    iteration = iteration + 1;
                    
                else
                    % were picking the best patch and updating it
                    
                    if( length( i ) > 1 )
                        % bestidx'
                        goodidxs = ( bestidx > 0 );
                        patchParams(goodidxs) = bestidx(goodidxs);
                    else
                        patchParams( i ) = bestidx;
                    end
                    
                    if( ~isempty( this.intXfmModelType ) && ~modelsSet )
                        if( length( model ) > 1 )
                            this.paramModels = model;
                        else
                            this.paramModels{ i } = model;
                        end
                    end
                    
                    [theseCosts,thisPatch] = this.totalPatchCost( patchParams );
                    
                    if( record )
                        patchesByIter{ iteration } = thisPatch;
                    end
                    
                    if( numel( theseCosts ) == 1 )
                        costs(iteration) = theseCosts;
                    else
                        costs(iteration) = theseCosts( bestidx );
                    end
                    
                    % has the cost changed much?
                    if( abs( costs(iteration) - lastCost  ) < this.convEps )
                        itersAtThisCost = itersAtThisCost + 1;
                    else
                        lastCost = costs(iteration);
                        fprintf('updated last cost to: %f\n', lastCost);
                        itersAtThisCost = 0;
                    end
                    
                    b = this.pc.updateConstraints( this.D2d, b, i, bestidx );
                end
                
                % converged if we've been at the same cost for awhile
                if( itersAtThisCost == this.convIters )
                    converged = 1;
                end
                
                % force exit after max iteration count
                if(iteration >= this.maxIters)
                    fprintf('hit max iters!\n');
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
            
            if( ~isempty( this.paramModels ))
                otherParams = this.paramModels;
            end
            
        end % build3dPatch
        
        function saveSupplemental( this, f_out )
            N = size(f_out,1);
            this.savedParams = cell( N , 1 );
            for i = 1:N
                this.savedParams{i} = f_out{i}{1}{6};
            end
        end
        
        function [ gif_fn ] = recordToMovie( this, patchesByIter, destdir )
            
            if( ~exist( destdir, 'dir'))
                mkdir( destdir );
            end
            
            N = length( patchesByIter );
            figure('color','w','WindowStyle','docked');
            for i = 1:N
                
                if( isempty( patchesByIter{i} ))
                    continue;
                end
                im = reshape( patchesByIter{i}, this.sz3d );
                
                nz = this.sz3d(3);
                dispsz = round( sqrt( nz ));
                dispsz = [ dispsz, ceil(nz./dispsz(1)) ];
                imdisp3d( im, 'size', dispsz ) ;
                
                % write to file here
                export_fig(destdir, sprintf('%s/im_%05d.jpg',destdir,i),'-nocrop','-painters');
                clf('reset');
                
            end
            
            % temporarily modify LD_LIBRARY_PATH for ImageMagick
            old_path = getenv('LD_LIBRARY_PATH');
            setenv('LD_LIBRARY_PATH', '/usr/local/cuda/lib64:');
            
            % make the gif
            gif_fn = sprintf('%s.gif',destdir);
            gifdelay = 1
            cmd = sprintf('convert -delay %d -loop 0 %s/*.jpg %s', ...
                gifdelay, destdir, gif_fn);
            fprintf('%s\n', cmd);
            system(cmd);
        end
        
        function [ cost, x ] = totalPatchCost( this, params )
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            b = this.constraintValue( params );
            x = cmtxInv * b;
            cost = norm( cmtx * x - b );
        end
         
        function randomCoords = genRandomCoords( this, N, exclude, extra )
            % randomly generate the patch location
            % that will be updated at each iteration
            %
            if( ~exist('extra','var') || isempty( extra ))
                extra = 1;
            end
            
            switch( this.updateOrder )
                
                case 'each'
                    randomCoords = randi( N, extra.*this.maxIters, 1 );
                    if( exist('exclude','var') && ~isempty( exclude ))
                        keepInds = 1:this.pc.numLocs;
                        keepInds = keepInds( ~exclude );
                        randomCoords = keepInds( randomCoords );
                    end
                case 'all'
                    randomCoords = zeros( extra.*this.maxIters, N );
                    for i = 1:this.maxIters
                        randomCoords( i, : ) = randperm( N );
                    end
                otherwise
                    error('invalid update order');
            end
        end
        
        function [iniParams, iniModels] = iniParamsDist( this, N, method )
            iniModels = [];
            switch( method )
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
                    
                    
                case 'hr'

                    if( isempty( this.Dini ))
                       error('Dini is empty!!'); 
                    end
                    
                    numIni = size( this.Dini, 1 );
                    iniParams = mat2cell( randi( numIni, N, 1 ),...
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
        
        function [ idxList, distList, cmtx ] = goodPatchConfig( this, b, iList )
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            
            % initialize experimental constraint values
            bexp = b;
            
            x = cmtxInv * bexp;
            curdist = norm( cmtx * x - bexp );
            
            num = length( iList );
            if( num > 1 )
                idxList = zeros( length( iList ), 1);
                distList    = zeros( length( iList ), 1);
            else
                idxList  = [];
                distList = [];
            end
            
            for i = iList 

                % the range of constraint values that will change
                % depending on the patch being tested
                brng = this.pc.constraintVecSubsets(i,:);
                
                testOrder = randperm( this.numDict );
                for nn = 1:this.numDict
                    n = testOrder( nn );
                    bexp( brng ) = this.D2d(n,:);
                    
                    x = cmtxInv * bexp;
                    dist = norm( cmtx * x - bexp );
                    
                    if( dist < curdist )
                       
                       if ( num > 1 )
                           idxList( i )  = n;
                           distList( i ) = dist;
                       else
                           idxList  = n;
                           distList = dist;
                       end
                       
                       break;
                    end
                end % over dictionary elements
            end % iList  
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
        
        function [ paramList, modelList, pvList, patchList ] = fitParamsToHR_dist( this, ini, doPatchEst, constrainScalesPos, constrainHrMinMax  )
            global DICTPATH;
            
            if( ~exist( 'doPatchEst', 'var' ) || isempty( doPatchEst ))
                doPatchEst = 1;
            end
            
            % check that we have a valid initialization to go with
            if( isscalar( ini ))
                if( isempty( this.Dini ))
                    error( 'Dini and ini param are empty - cannot run fitParamsToHR_dist');
                end
                iniParam = this.iniParamsDist( ini, 'hr');
            else
                iniParam = ini;
            end
            num_jobs = length( iniParam );
            
            this.save();
            
            fun = @run_obj_method_dist;
            use_gpu = 0;
            num_out_args = 4; % fitParamsToHR has 4 output args
            run_script = fullfile( DICTPATH, 'bin', 'my_run_runObjMethod.sh');
            
            varargin = {  repmat( {this.obj_fn}, num_jobs, 1), ...
                repmat( {'this'}, num_jobs, 1), ...
                repmat( {'fitParamsToHR'}, num_jobs, 1), ...
                repmat( {num_out_args}, num_jobs, 1), ...
                iniParam, ...
                repmat( {doPatchEst}, num_jobs, 1), ...
                repmat( {constrainScalesPos}, num_jobs, 1), ...
                repmat( {constrainHrMinMax}, num_jobs, 1) ...
                };
            
            f_out = qsub_dist(fun, 1, use_gpu, ...
                [], [], run_script, ...
                varargin{:} );
            
            paramList = cell( num_jobs, 1 );
            modelList = cell( num_jobs, 1 );
            pvList    = cell( num_jobs, 1 );
            patchList = cell( num_jobs, 1 );
            for i = 1:size(f_out,1)
                paramList{i} = f_out{i}{1}{1};
                modelList{i} = f_out{i}{1}{2};
                pvList{i}    = f_out{i}{1}{3};
                patchList{i} = f_out{i}{1}{4};
            end
            
        end
        
        function [ patchParams, modelList, pv, patch ] = fitParamsToHR( this, x, doPatchEst, constrainScalesPos, constrainHrMinMax )
            patchParams = zeros( this.pc.numLocs, 1 );
            modelList   = cell ( this.pc.numLocs, 1 );
            pv = [];
            patch = [];
            
            x = vecToRowCol( x, 'col');
            for j = 1:this.pc.numLocs
                
                fprintf('fitting model for location %d of %d\n', ...
                            j, this.pc.numLocs );
                
                [ idx, ~, model ] = this.fitIdxAndModel( j, x );
                patchParams(j) = idx;
                
                if( ~isempty(this.intXfmModelType))
                    modelList{ j } = model;
                end
            end
            
            if( doPatchEst )
                if( nargout >= 4 )
                    [ pv, patch ] = this.patchFromParamsScale( patchParams, constrainScalesPos, constrainHrMinMax  );
                else
                    [ pv ] = this.patchFromParamsScale( patchParams, constrainScalesPos, constrainHrMinMax  );
                end
            end
        end
        
        function [ patchParams, dists, modelList, f_out ] = fitIdxAndModel_dist( this, x, doBest, returnList )
            global DICTPATH;
            
            if( ~exist('returnList','var') || isempty( returnList ))
                returnList = false;
            end
            
            num_jobs = this.pc.numLocs;
            
            if( returnList )
                patchParams = zeros( num_jobs, 1);
                dists       = zeros( num_jobs, 1);
                modelList   = cell( num_jobs, 1);
            else
                patchParams = cell( num_jobs, 1);
                dists       = cell( num_jobs, 1);
                modelList   = cell( num_jobs, 1);
            end
            
            this.save();
            
            fun = @run_obj_method_dist;
            use_gpu = 0;
            num_out_args = 3; % fitIdxAndModel has 2 output args
            run_script = fullfile( DICTPATH, 'bin', 'my_run_runObjMethod.sh');
            
            varargin = {  repmat( {this.obj_fn}, num_jobs, 1), ...
                repmat( {'this'}, num_jobs, 1), ...
                repmat( {'fitIdxAndModel'}, num_jobs, 1), ...
                repmat( {num_out_args}, num_jobs, 1), ...
                num2cell( 1:num_jobs ), ...
                repmat( {x}, num_jobs, 1), ...
                repmat( {doBest}, num_jobs, 1), ...
                repmat( {returnList}, num_jobs, 1 ) ...
                };
            
            f_out = qsub_dist(fun, 1, use_gpu, ...
                [], [], run_script, ...
                varargin{:} );
           
            if( returnList )
                for i = 1:num_jobs
                    patchParams{i} = f_out{i}{1}{1};
                    dists{i}       = f_out{i}{1}{2};
                    modelList{i}   = f_out{i}{1}{3};
                end                
            else
                for i = 1:num_jobs
                    patchParams(i) = f_out{i}{1}{1};
                    dists(i)       = f_out{i}{1}{2};
                    modelList{i}   = f_out{i}{1}{3};
                end
            end
            
        end
            
        function [ idx, curdist, models ] = fitIdxAndModel( this, i, x, doBest, returnList, xMsk )
            if( ~exist( 'doBest', 'var' ) || isempty( doBest ))
                doBest =1;
            end
            if( ~exist( 'returnList', 'var' ) || isempty( returnList ))
                returnList = false;
            end
            if( ~exist( 'xMsk', 'var' ))
                xMsk = [];
            end
%             if( nargout >=4 )
%                allDists = zeros( d23.numDict, 1 );
%             end
            
%             % v1
%             rng  = this.pc.constraintVecSubsets(i,:);
%             cmtx = this.pc.cmtx;
%             Ax   = cmtx * x;
%             AxR  = Ax( rng );

%             % v2
%             rng  = this.pc.constraintVecSubsets(i,:);
%             cmtx = this.pc.cmtx(rng,:);
%             AxR  = cmtx * x;

%           % v3
            if( isempty( xMsk ))
                AxR = this.pc.patchProject( i, x );
            else
                AxR = x;
            end
                
            
            if( var( AxR ) < 0.0001 )
                AxR = AxR + 0.0001.*randn(size(AxR));
            end
            
            if( returnList )
                distList = zeros( this.numDict, 1 );
                if( ~isempty(this.intXfmModelType))
                    models = cell( this.numDict, 1 );
                else
                    models = [];
                end
            else
                idx = [];
                models = [];
            end
            curdist = Inf;
            
            for n = 1:this.numDict
                
                if( mod( n, this.verboseEvery) == 0 )
                    fprintf(' Dict2dTo3dSampler.fitIdxAndModel %d of %d\n', n, this.numDict);
                end
                
                bexp = this.D2d(n,:)';
%                 if( length( bexp ) ~= length( AxR ))
%                     bexp = bexp( rng );
%                 end

                if( ~isempty( xMsk ))
                    if( isscalar( xMsk ))
                        bexp = PatchConstraints.downsampleByMaskDim( bexp, this.pc.planeMaskLRI( i ) );    
                    else
                        bexp = PatchConstraints.downsampleByMaskDim( bexp, xMsk );
                    end
                    bexp = bexp(:);
                end
                
                
                if( var( bexp ) < 0.0001 )
                    bexp = bexp + 0.0001.*randn(size(bexp));
                end
                    
                if( ~isempty(this.intXfmModelType))
                    thismodel = fit( bexp, AxR, this.intXfmModelType, this.fitParams{:} );
                    dist = norm( AxR - feval( thismodel, bexp ) );
                else
                    dist = norm( AxR - bexp );
                end
                
                if( returnList )
                    distList( n ) = dist;
                     if( ~isempty(this.intXfmModelType))
                        models{ n } = thismodel;
                     end
                elseif( (dist < curdist) )
                    idx   = n;
                    curdist = dist;
                    if( ~isempty(this.intXfmModelType) )
                        models = thismodel;
                    end
                    if( ~doBest )
                        return;
                    end
                end
                
            end
            
            if( returnList )
                %fprintf('sorting output to fitIdxAndModel\n');
                [ distList, idx ] = sort( distList );
                curdist = distList;
                if( ~isempty(this.intXfmModelType))
                    models = models( idx );
                end
            end
            
        end

        function [ idx, dist, modelList ] = goodPatchConfigModel( this, b, iList, doBest )
            
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            
            % initialize experimental constraint values
            bexp = b;
            
            x = cmtxInv * bexp;
            Ax  = cmtx * x;
          
            num = length( iList );
            
            
            if( num > 1 )
                idx = zeros( num, 1 ); 
                modelList = cell( num, 1 );
            else
                idx = [];
                modelList = [];
            end
            
            k = 1;
            for i = iList
                % the range of constraint values that will change
                % depending on the patch being tested
                rng = this.pc.constraintVecSubsets(i,:);

                AxR = Ax( rng );
                if( var( AxR ) < 0.0001 )
                    AxR = AxR + 0.0001.*randn(size(AxR));
                end
                        
                curdist = norm( AxR - feval( this.paramModels{i}, bexp( rng )) );
                
                % model = this.paramModels{i}; % the current model
                
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
                        
                        if( num > 1 )
                           idx( i )  = n;
                           modelList{ i } = thismodel;
                        else
                            idx = n;
                            modelList = thismodel;
                        end
                        
                        curdist = dist;
                        if( ~doBest )
                            break;
                        end
                    end
                    
                end % loop
                
                k = k + 1;
                
            end % iList
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
        
        function [ bestidx, sims, cmtx ] = bestPatchConfig( this, b, iList )
            cmtx    = this.pc.cmtx;
            cmtxInv = this.pc.cmtxInv;
            
            % initialize experimental constraint values
            bexp = b;
            
            num = length( iList );
            if( num > 1 )
                bestidx = zeros( length( iList ), 1);
            end
            
            for i = iList
                % the range of constraint values that will change
                % depending on the patch being tested
                brng = this.pc.constraintVecSubsets(i,:);
                
                sims = zeros( this.numDict, 1 );
                for n = 1:this.numDict
                    bexp( brng ) = this.D2d(n,:);
                    
                    x = cmtxInv * bexp;
                    sims( n ) = norm( cmtx * x - bexp );
                end
                
                [ ~, bi ] = min( sims );
                
                if ( num > 1 )
                   bestidx( i ) = bi;
                else
                    bestidx = bi; 
                end

            end
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
    
        function b = constraintValue( this, obj, models )
            if( ~exist('models','var'))
                models = {};
            end
            
            if( isa( obj, 'net.imglib2.algorithms.opt.astar.SortedTreeNode'))
                b = this.pc.constraintValueNode( this.D2d, obj, models );
            else
                b = this.pc.constraintValueList( this.D2d, obj, models );
            end
            
            if( ~isempty( this.intXfmModelType ))
                for i = 1:this.pc.numLocs
                    rng = this.pc.constraintVecSubsets(i,:); 
                    b( rng ) = feval(models{i}, b( rng ));
                end
            end
            if( this.scaleDictElems )
                fprintf('scaling constraints\n');
                for i = 1:this.pc.numLocs
                    rng = this.pc.constraintVecSubsets(i,:); 
                    b( rng ) = this.paramScales(i) .* b( rng );
                end
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
       
        function [ distList, models ] = toFeatures_dist( this, X, num_jobs, isLrObs, doLoad )
            % rows of X are observations
            
            global DICTPATH;
            
            if( ~exist('isLrObs','var') || isempty(isLrObs))
               isLrObs = false; 
            end
            if( ~exist('doLoad','var') || isempty(doLoad))
               doLoad = true; 
            end
            if( isLrObs )
                funName = 'fitToLowResObservation'
            else
                funName = 'toFeatures'
            end
            
            N = size( X, 1 );
            M = size( X, 2 );
            
            if( ~exist('num_jobs','var') || isempty(num_jobs))
                num_jobs = N;
                
            else
                [ ranges, tasksPerWorker ] = divideTasks( N, num_jobs );
            end
            
            % split up the observations of n
            Xarg = mat2cell( X, tasksPerWorker, M);
            
            
            distList = zeros( N, this.pc.numLocs, this.numDict);
            size( distList )
            
            if( ~isempty(this.intXfmModelType))
                models = cell( N, this.pc.numLocs, this.numDict );
            else
                models = {};
            end
            
            this.save();
            
            fun = @run_obj_method_dist;
            use_gpu = 0;
            num_out_args = 2; % toFeatures has 2 output args
            run_script = fullfile( DICTPATH, 'bin', 'my_run_runObjMethod.sh');
            
            varargin = {  repmat( {this.obj_fn}, num_jobs, 1), ...
                repmat( {'this'}, num_jobs, 1), ...
                repmat( {funName}, num_jobs, 1), ...
                repmat( {num_out_args}, num_jobs, 1), ...
                Xarg, ...
                };
            
            if( doLoad )
                f_out = qsub_dist(fun, 1, use_gpu, ...
                    [], [], run_script, ...
                    varargin{:} );
                
                for i = 1:num_jobs
                    
                    range = ranges{i};
                    distList( range(1):range(2), :, : ) = f_out{i}{1}{1};
                    
                    if( ~isempty(this.intXfmModelType))
                        models(ranges(1):ranges(2), :, : ) = f_out{i}{1}{2};
                    end
                end
            else
                 f_out = qsub_dist_noLoad(fun, 1, use_gpu, ...
                    [], [], run_script, ...
                    varargin{:} );
%                  fprintf( '%s\n', f_out{1} );
                 
                 distList = f_out;
                 models = {};
            end
            
        end
        
        function [ distList, models ] = toFeatures( this, X, xMsk, isObs )
            % rows of X are observations
            
            if( ~exist( 'xMsk', 'var' ))
                xMsk = [];
            end
            if( ~exist( 'isObs', 'var' ))
                isObs = false;
            end
            
            numObs = size( X, 1 );
            numLocs = this.pc.numLocs;
            
            if( ~isempty(this.intXfmModelType))
                models = cell( numObs, numLocs, this.numDict );
            else
                models = [];
            end
            
            if( isObs )
                
                %  [ dictIdxs, dictCosts, models, x, ~, ~ ] = this.fitToLowResObservation( x, i  );
                %                 distList( j, i, : ) = dictCosts;
                
                [  distList, models ] = this.fitToLowResObservation( X );
                
            else
                
                distList = zeros( numObs, numLocs, this.numDict);
                
                for j = 1:numObs
                    x = X(j,:)';
                    
                    for i = 1:numLocs
                        
                        % get the sub-patch of the observation
                        if( isempty( xMsk ))
                            AxR = this.pc.patchProject( i, x );
                        else
                            AxR = x;
                        end
                        
                        if( var( AxR ) < 0.0001 )
                            AxR = AxR + 0.0001.*randn(size(AxR));
                        end
                        
                        for n = 1:this.numDict
                            
                            if( mod( n, this.verboseEvery) == 0 )
                                fprintf(' Dict2dTo3dSampler.fitIdxAndModel %d of %d\n', n, this.numDict);
                            end
                            
                            bexp = this.D2d(n,:)';
                            
                            % downsample the dictionary if necessary
                            if( ~isempty( xMsk ))
                                if( isscalar( xMsk ))
                                    bexp = PatchConstraints.downsampleByMaskDim( bexp, this.pc.planeMaskLRI( i ) );
                                else
                                    bexp = PatchConstraints.downsampleByMaskDim( bexp, xMsk );
                                end
                                bexp = bexp(:);
                            end
                            
                            % compute the distance
                            if( ~isempty(this.intXfmModelType))
                                thismodel = fit( bexp, AxR, this.intXfmModelType, this.fitParams{:} );
                                models{ j, i, n } = thismodel;
                                distList( j, i, n ) = this.comparator.distance( AxR, feval( thismodel, bexp ));
                            else
                                distList( j, i, n ) = this.comparator.distance( AxR, bexp );
                            end
                            
                        end
                        
                    end % locations loop
                    
                end % observations loop
            end % isObs?
        end % function
        
        function [ dictCosts, models ] = fitToLowResObservation( this, xin  )
            % observations should be rows of xin
            
            nLocs = this.pc.numLocs;
            N = size( xin, 1 );
            dictCosts = zeros( N, nLocs, this.numDict );
            
            models = [];
            if( ~isempty(this.intXfmModelType))
                models = cell( N, nLocs, this.numDict );
            end
            
            
            for n = 1 : N
                for i = 1 : nLocs
                    [ dictIdxs_n, dictCosts_n, models_n ] = bestKdicts( this, xin(n,:), i, this.numDict );
                    dictCostsTmp = dictCosts_n(1 + length(dictCosts_n) - dictIdxs_n(end:-1:1));
                    dictCosts( n, i, : ) = dictCostsTmp;
                    
                    if( ~isempty(this.intXfmModelType))
                        models( n, i, : ) = models_n; %#ok<AGROW>
                    end
                end
            end
        end
        
        function [ dictIdxs, dictCosts, models, x, msk, didTp ] = bestKdicts( this, xin, i, K  )
            % expects xin to be a row vector
            
            %             if( ~exist( 'isLR', 'var' ) || isempty( isLR ))
            %                 isLR = true;
            %             end
            
            % make sure x is a row vector
            x = vecToRowCol( xin, 'row');
            
            dists  = nan( this.numDict, 1 );
            models = {};
            
%             msk = this.pc.planeMaskLRI( i );
%             [x,~,didTp] = PatchConstraints.downsampleByMaskDim( x(:), msk );
%             x = [x(:)]';

            [x, ~, msk, ~, didTp ] = this.pc.getSubImage( x, i );
            x = x(:)';
            
            if( numel( x ) == size( this.D2d, 2 ) )
                D = this.D2d;
            else
                D = this.D2d_downsampled;
            end
            
            % check size of K
            K = min( K, size( D, 1 ));
            
            if( ~isempty(this.intXfmModelType))
                fprintf('fitting model: %s\n', this.intXfmModelType);
                models = cell( this.numDict, 1 );
                for n = 1 : this.numDict
                    bexp      = D(n,:);
                    models{n} = fit( bexp', x', this.intXfmModelType, this.fitParams{:} );
                    this.comparator.distance( x', feval( models{n}, bexp' ) );
                end
            else
                dists =  this.comparator.distances( x, D );
            end
            
            [ dictCosts, is ] = sort( dists );
            dictCosts = dictCosts(1:K);
            
            dictIdxs = is( 1:K );
            
            if( ~isempty(this.intXfmModelType))
                models = models( is(1:K) );
            end
        end
    
        function [ dictParams, xhat, dictCost, xhatup ] = dictParamsLasso( this, xin, loci, params )
            % uses the spams toolbox's implementation of lasso to estimate
            % a sparse coefficient vector that explains the observation 'x'
            % at location 'loci' using the 2d dictionary
            
            x = vecToRowCol( xin, 'row');
%             msk = this.pc.planeMaskLRI( loci );
%             [x,~,didTp] = PatchConstraints.downsampleByMaskDim( x(:), msk );
            [x, ~, ~, ~, ~ ] = this.pc.getSubImage( x, loci );
            x = x(:);
            
            isLowRes = false;
            if( numel( x ) == size( this.D2d, 2 ) )
                D = this.D2d;
            else
                D = this.D2d_downsampled;
                isLowRes = true;
            end
            
            dictParams = mexLasso( x, D', params );
            if( nargout > 1 )
                xhat = D' * dictParams;
            end
            if( nargout > 2 )
                dictCost =  norm( x - xhat );
            end
            if( nargout > 3 )
                if( isLowRes )
                    xhatup = this.D2d' * dictParams;
                else
                    % if the observation is high res, 
                    % then the current xhat is already high res
                    xhatup = xhat; 
                end
            end
        end
           
    end
    
    methods( Static )
        
        function DiniHR = upsampleInitialObservations( DiniLR, iniSz, outSz, dsFactor, interpOpts )
            Nini = size( DiniLR, 1 );    
            DiniHR = zeros( Nini, prod(outSz));
            for i = 1:Nini
                imtmp = reshape( DiniLR( i, :), iniSz );
                im_re = upsampleObservationForIni( imtmp, iniSz, dsFactor, interpOpts{:} );
                DiniHR(i,:) = im_re(:); 
            end
        end
        
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
