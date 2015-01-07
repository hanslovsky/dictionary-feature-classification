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
       paramModels;
       stopAfterFitParamIni = 1;
       
       scaleDictElems = 0;
       paramScales;
       
       recordParamsOverIters = 0;
       savedParams;

%        updateOrder = 'all';
       updateOrder = 'each';
       
    end
    
    properties( SetAccess = protected )
       
    end
    
    methods
        
        % Constructor
        % D2d - 
        % sz  - size of 2d patches
        % f   - downsampling factor
        function this = Dict2dTo3dSampler( D2d, sz, f, overlappingPatches, scaleByOverlap )
            this = this@Dict2dTo3d( D2d, sz, f, overlappingPatches, scaleByOverlap );
            if( this.useSubset )
                fprintf('Computing insersection matrix inverses...');
                this.pc.compXsectInverses();
                fprintf('.done\n');
            end
            
            if( this.scaleByOverlap )
                this.paramScales = this.pc.overlapFraction;
            end
                
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
        
        function [ paramList, modelList, pvList ] = fitParamsToHR_dist( this, ini, doPatchEst, constrainScalesPos, constrainHrMinMax  )
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
            num_out_args = 3; % fitParamsToHR has 2 output args
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
            pvList = cell( num_jobs, 1 );
            for i = 1:size(f_out,1)
                paramList{i} = f_out{i}{1}{1};
                modelList{i} = f_out{i}{1}{2};
                pvList{i}    = f_out{i}{1}{3};
            end
            
        end
        
        function [ patchParams, modelList, pv ] = fitParamsToHR( this, x, doPatchEst, constrainScalesPos, constrainHrMinMax )
            patchParams = zeros( this.pc.numLocs, 1 );
            modelList   = cell ( this.pc.numLocs, 1 );
            pv = [];
            
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
               [ pv ] = this.patchFromParamsScale( patchParams, constrainScalesPos, constrainHrMinMax  );
            end
        end
        
        function [ idx, curdist, model ] = fitIdxAndModel( this, i, x, doBest )
            if( ~exist( 'doBest', 'var' ) || isempty( doBest ))
                doBest =1;
            end
            
            rng  = this.pc.constraintVecSubsets(i,:);
            cmtx = this.pc.cmtx;
            Ax   = cmtx * x;
            AxR  = Ax( rng );
            
            if( var( AxR ) < 0.0001 )
                AxR = AxR + 0.0001.*randn(size(AxR));
            end
            
            idx = [];
            model = [];
            curdist = Inf;
            for n = 1:this.numDict

                bexp = this.D2d(n,:)';
%                 if( length( bexp ) ~= length( AxR ))
%                     bexp = bexp( rng );
%                 end
                
                if( var( bexp ) < 0.0001 )
                    bexp = bexp + 0.0001.*randn(size(bexp));
                end
                    
                if( ~isempty(this.intXfmModelType))
                    thismodel = fit( bexp, AxR, this.intXfmModelType );
                    dist = norm( AxR - feval( thismodel, bexp ) );
                else
                    dist = norm( AxR - bexp );
                end
                
                if( (dist < curdist) )
                    idx   = n;
                    curdist = dist;
                    if( ~isempty(this.intXfmModelType))
                        model = this.paramModels;
                    end
                    if( ~doBest )
                        return;
                    end
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
            end
            if( this.scaleDictElems )
                fprintf('scaling constraints\n');
                for i = 1:this.pc.numLocs
                    rng = this.pc.constraintVecSubsets(i,:); 
                    b( rng ) = this.paramScales(i) .* b( rng );
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
