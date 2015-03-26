classdef Tid < handle
    % Tid - transformationally invariant dictionary
    %
    % John Bogovic
    % HHMI
    % August 2014

    properties ( SetAccess = private )

        alpha;  % The dictionary coefficients
        
        numSamples;
        numVars;
        
        iMergeDxfms;    % indicates which original dictionary element (D)
                        % the elements of Dxfm came from
        mergeFun = @max;
        
        params; % dictionary learning parameters
        
        % distance parameters
        distanceFun;
        
        % transformation helpers
        tSubPatches;

        patchSize;
        dictPatchSize;

        rotRng;
        rotXfmIdx;
        
        dictFun;

        % persistence
        obj_fn;
    end
    
    properties ( Transient, SetAccess = private )
        X;      % The data
        D;      % The dictionary elements
        Dxfm;   % The full dictionary
    end
    
    properties
        distTol  = 0.02;
        
        % optimization parameters
        bigIters = 1;
        nJobs    = 20;
    end

    methods

        % Constructor
        function this = Tid( X, patchSize, dictPatchSize, params, distTol )

            if( ~exist( 'distTol', 'var') || isempty(distTol))
                this.distTol = 0.1;    
            else
                this.distTol = distTol;
            end
            this.obj_fn = [];
                
            this.X = X;  % the data
            [ this.numSamples, this.numVars ] = size( X );

            % check that patchSize is compatible with dictionary size
            if ( prod( patchSize ) ~= this.numVars )
                error( 'patch size not compatible with dictionary size');
            end
            
            this.patchSize = patchSize;
            if( isempty( dictPatchSize ))
                this.dictPatchSize = this.patchSize;
            else
                this.dictPatchSize = dictPatchSize;
            end
            
            if( ~exist('params','var') || isempty( params ))
                this.params = Tid.defaultParams();
            else 
                this.params = params;
            end

            genVectorTransformations( this );
            
            this.dictFun = Tid.getDictionaryBuildingFunction();
            this.setComparator();
            
            % TODO - this assumes that patchSizes are odd.
            % should we check this? if, so where?
            this.tSubPatches = Tid.translationSubPatches( this.patchSize, this.dictPatchSize );
            
        end % constructor

        function obj_fn = save( this, fn )
            global DFEVAL_DIR;
            if( exist( 'fn', 'var') && ~isempty( fn ) )
                this.obj_fn = fn;
            elseif( isempty( this.obj_fn ))
                this.obj_fn = fullfile( DFEVAL_DIR, sprintf('Dict2dTo3d_%s.mat', ...
                    datestr( now, 30 ) ));
            end
            
            obj_fn = this.obj_fn;
            if( ~exist( obj_fn, 'file'))
                save( this.obj_fn, 'this' ); 
            end
        end

        function thecopy = copy( this )
            thecopy = Tid( this.X, this.patchSize, this.params, this.distTol );
            thecopy.X = [];
        end

        function setComparator( this, method )
            if( ~exist( 'method','var'))
                method = 'euc';
            end
            this.distanceFun = PatchCompare( method );
            
        end
        
        % Very experimental
        function D = buildDictionary(this)
            
            this.D = this.dictFun( this.X, this.params );

            for iter = 1:this.bigIters
               
                % initialize with previous dictionary
                this.params.D = this.D;

                D = this.D;
%                 this.makeDictRotInv();
            end

            this.Dxfm = this.xfmDict( );

        end % the magic

        function alpha = getFeatures( this, X, f, params )
            
            if( ~exist('f','var'))
                f = [];
            end
            if( ~exist('params','var') || isempty(params))
                prm = this.params;
            else
                prm = params;
            end
            
            if( ~isempty( f ))
               ds = Tid.getDownsampler3dzRs();
               dict = ds( this.Dxfm, this.patchSize, f ); 
               dict = dict ./ repmat( sqrt(sum(dict.*dict)), size(dict,1), 1);
            else
               dict = this.Dxfm; 
            end
            
            try
                alpha = mexLasso( X, dict, prm );
            catch e
                disp('no mexLasso - random features');
                %e % print error 
                alpha = randn( size(X,2), size(dict,2) );
            end
            
        end
        
        function alpha = trainingFeatures( this )
            alpha = getFeatures( this, this.X );
            this.alpha = alpha;
        end

        function alpha = getFeaturesOrig( this, X, f, params )
            
            if( ~exist('f','var'))
                f = [];
            end
            if( ~exist('params','var') || isempty(params))
                prm = this.params;
            else
                prm = params;
            end
            
            if( ~isempty( f ))
               ds = Tid.getDownsampler3dzRs();
               % szds = this.patchSize ./ [1 1 f];
               dict = ds( this.D, this.patchSize, f );
               dict = dict ./ repmat( sqrt(sum(dict.*dict)), size(dict,1), 1);
            else
               dict = this.D; 
            end
            
            try
                alpha = mexLasso( X, dict, prm );
            catch e
                disp('no mexLasso - random features ');
                %e % print error 
                alpha = randn( size(X,2), size(dict,2) );
            end
            
        end
        
        function alpha = trainingFeaturesOrig( this )
            alpha = getFeaturesOrig( this, this.X );
            this.alpha = alpha;
        end

        % Make the dictionary rotationally invariant
        function makeDictRotInv( this )

            dxsz  = size( this.D, 2); % size of transformed dictionary
            rxi= cell2mat(this.rotXfmIdx);
            
            allIdxs = 1:dxsz;
            invariantIdxs = true(1,dxsz);
            
            for i = allIdxs
                % check if this index is still valid
                % it may have been removed at a previous iteration
                if( ~invariantIdxs(i) )
                    fprintf('skipping %d\n',i);
                    continue;
                end
                
                j = setdiff( allIdxs, i );
                
                thisXfmDict = this.D( :, i );
                thisXfmDict = thisXfmDict( rxi );
                Draw = this.D( :, j );
                
                pd = min( pdist2( thisXfmDict', Draw' ), [], 1);
                similar = ( pd < this.distTol );
                
                fprintf('excluding %d elems at %d\n', nnz(similar), i );
                
                % update indices of dict elems that are
                % transformations of this dict element
                invariantIdxs( allIdxs(j(similar)) ) = false;
                
            end

            this.D = this.D( :, invariantIdxs );
            
        end % make dict rot inv

        function Dxfm = xfmDict( this )
            [r,N] = size( this.D );
            M = length( this.rotXfmIdx );
            Dxfm = zeros( r, N*M );
            j = reshape( 1:(N*M), M, [] );
            k =  cell2mat( this.rotXfmIdx );
            
            for i = 1:N
                tmp = this.D(:,i);
                Dxfm(:,j(:,i)) = tmp(k);
            end
            
            this.iMergeDxfms = reshape(repmat( 1:N, M, 1 ), [], 1);
        end

        function alphaM = mergeFeatures( this, alpha )
            if( isempty( this.iMergeDxfms))
                alphaM = [];
                return;
            end
            
            N =  size(this.D, 2);
            S = size(alpha,2); % number of samples
            alphaM = zeros( N, S ); 
            
            for i = 1:N
                alphaM(i,:) = this.mergeFun( alpha( (this.iMergeDxfms == i), :));
            end
            
        end

        function genVectorTransformations( this )
            % possible transformations include

            % comupte the index look-ups for this transformation 
            % dont reshape them.
            
            % check if we're working in 2d or 3d
            patchSzTmp = this.patchSize( this.patchSize > 1);
            ndim = length( patchSzTmp );
            
            if( ndim == 3 )
                this.rotXfmIdx = xfmToIdx( cubeSymmetry(), patchSzTmp, 0 );
            elseif (ndim == 2)
                this.rotXfmIdx = xfmToIdx( squareSymmetry(), patchSzTmp, 0 );
            end
            
            if( ~isempty( this.rotRng ) )
                this.rotXfmIdx = this.rotXfmIdx( this.rotRng );
            end

        end
        
        function [D, dictiterList, simRiList ] = buildTranslationInvariantDictionary( ...
                this, iniD, outputDictAtEveryIter )
            
            wasInitialization = true;
            if( ~exist( 'iniD','var') || isempty( iniD ))
                iniD = this.X( randperm( this.numSamples, this.params.K ), : );
                wasInitialization = false;    
            end
            this.D = iniD;
            
            if( ~exist( 'outputDictAtEveryIter','var') || isempty( outputDictAtEveryIter ))
                outputDictAtEveryIter = false;
            end
            
            drr = Drr( this );
            
            % initial optimization
            if( ~wasInitialization )
                % set initial dictionary
                this.params.D = this.D';
                fprintf('start of initial opt\n');
                this.D = mexTrainDL( this.X', this.params );
                this.D = this.D'; % make sure dictionary elements are in the rows
                fprintf('end of initial opt\n');
            end
            
            %             try
            %                   this.makeDictRotInv();
            %             catch e
            %                 disp('no mexTrainDL - choosing dictionary by random sampling');
            %                 %e % print error
            %                 this.D = this.X( :, randperm( this.numSamples, this.params.K));
            %                 this.makeDictRotInv();
            %             end
            
            if( outputDictAtEveryIter )
                k = 1;
                m = 1;
                if( ~wasInitialization )
                    dictiterList = cell( this.bigIters + 1, 1 );
                    simRiList = cell( this.bigIters, 1 );
                    
                    dictiterList{ k } = this.D;
                    k = k + 1;
                else
                    dictiterList = cell( this.bigIters, 1 );
                    simRiList = cell( this.bigIters - 1, 1 );
                end
            end
            
            for iter = 1:this.bigIters
                
                % initialize with previous dictionary
                % ( may be the input initial dictionary )
                
                % prune redundant entries
                %  THE GENERAL, BUT SLOW VERSION
                %                 Ddists = this.distanceFun.pairDistances( this.D, this.D );
                
                [Ddists] = this.translationInvariantHelper_dist( this.D, this.nJobs );
                ri = drr.prune( Ddists );
                
                this.D( ri, : ) = this.X( randperm( this.numSamples, nnz(ri) ), :);
                
                % optimize
                this.D = mexTrainDL( this.X', this.params );
                this.D = this.D'; % make sure dictionary elements are in the rows
                
                if( outputDictAtEveryIter )
                    dictiterList{ k } = this.D;
                    
                    simRiList{ m, 1 } = Ddists;
                    simRiList{ m, 2 } =  ri;
                    m = m + 1;
                
                    k = k + 1;
                end
                
                fprintf('end of opt %d\n', iter);
            end
            
            D = this.D;
        end % buildTranslationInvariant
        
        function [ newD, sim ] = translationInvariantHelperDeep( this, D )
            N = size( D, 1 );
            
%             N2 = N.*N;
%             NN =  max( N2, min(10, 0.05.*N2) );
%             sim = spalloc( N, N, NN );

            numSp = size(this.tSubPatches,1);
            mi = ( numSp - 1 )./2;
            
            bignum = this.distTol .* 1000;
            sim = bignum.*ones( N, N, numSp, numSp );
            
            % why not start m at (n+1)
            %   to compare different subpatches of the two
            for n = 1:N
                for j = 1:numSp
                    
                    midPatch = D( n, this.tSubPatches(j,:) );
                    
                    for m = 1:N
                        
                        if( m == n )
                            continue;
                        end
                        
                        for k = 1:numSp
                            dist = this.distanceFun.distance( midPatch, D(m,this.tSubPatches(k,:)) );
%                             fprintf('%d %d %d - %f\n', n, m , k, dist );
                            if( dist < this.distTol )
                                sim(n,m,j,k) = dist;
                            end
                        end
                    end
                end
            end
            
            fprintf('num nan: %d\n',nnz( isnan( sim )));
            
%             sim(sim == bignum) = 0;
%             aa = sum( sim, 3 )
%             sum(aa,1)

            newD = [];
%             simOut = sparse( sim );
        end
        
        function [ sim ] = translationInvariantHelper_dist( this, D, nJobs )
            global DICTPATH;
            
            N = size( D, 1 );
            
            % build ranges
            dimrangeList = cell( nJobs, 1 );
            bnds = ceil(linspace( 1, N+1, nJobs+1 ));
            for i = 1 : length(bnds) - 1
                dimrangeList{i} = bnds(i) : bnds(i+1)-1;
            end
            
            this.save();
            
            fun = @run_obj_method_dist;
            use_gpu = 0;
            num_out_args = 1; % translationInvariantHelperSparse has 1 output args
            run_script = fullfile( DICTPATH, 'bin', 'my_run_runObjMethod.sh');
           
            varargin = {  repmat( {this.obj_fn}, nJobs, 1), ...
                repmat( {'this'}, nJobs, 1), ...
                repmat( {'translationInvariantHelperSparse'}, nJobs, 1), ...
                repmat( {num_out_args}, nJobs, 1), ...
                repmat( {D}, nJobs, 1), ...
                dimrangeList
            };
                %num2cell( 1:nJobs ), ...
            
            f_out = qsub_dist(fun, 1, use_gpu, ...
                [], [], run_script, ...
                varargin{:} );

            % assume that cat is efficient with sparse matrices
            % and try to call it just once
            % * first grab the interesting elements and drop them into a 
            %   cell array

            tmp = cell( nJobs, 1 );
            sumsz = 0;
            for i = 1:nJobs
                tmp{i} = f_out{i}{1}{1};
                size(tmp{i})
                sumsz = sumsz + size( tmp{i}, 1);
            end
            
            sim = cat( 1, tmp{:} );
            
        end
        
        function [ sim ] = translationInvariantHelperSparse( this, D, dimrange )

            N = size( D, 1 );
            if( ~exist('dimrange','var') || isempty( dimrange ))
               dimrange = 1:N; 
            end
            M = length( dimrange );
            
            numSp = size(this.tSubPatches,1);
            
%             if( numSp == 1 )
%                 mi = 1;
%             else
%                 mi = ( numSp - 1 )./2;
%             end
            
            M2 = M*N;
            MM =  max( M2, min(10, 0.01.*M2) );
            sim = spalloc( M, N, MM );
            
            totOverlapRange =  2.*(this.patchSize - 1) + 1;
            szDiff = this.patchSize - this.dictPatchSize;
            tmpLo = ((totOverlapRange - 1)/2 - szDiff);
            tmpHi = ((totOverlapRange - 1)/2 + szDiff);
            
            [vorX, vorY] = ndgrid( tmpLo(1):tmpHi(1), tmpLo(2):tmpHi(2));
            validOverlapInds = sub2ind( totOverlapRange, vorX(:), vorY(:));
            
            % get range around center
            
            % TODO think about whether explicitly storing
            % the symmetrical elements of the matrix is worth it...
            %   storing them requires extra logic elsewhere
            %   not storing them is cheaper
            for n = 1:M
                
                if( mod(n,10) == 0 )
                    fprintf('%d\n', n);
                end

                % basePatch = reshape( D( dimrange(n), this.tSubPatches(mi,:) ), this.dictPatchSize );
                basePatch = reshape( D( dimrange(n), : ), this.patchSize );
                
                %for m = 1:N
                for m = (dimrange(n)+1) : N

%                     fprintf('n: %d, m: %d\n', dimrange(n), m);
                    
                    %if( m == n )
                    %    continue;
                    %end
                  
                    % first compute max cross-correlationa between dictionary elements
%                     xc = xcorr2( basePatch, reshape( D( m, : ), this.patchSize ));
                    xc = normxcorr2_general( basePatch, reshape( D( m, : ), this.patchSize ));
                    
                    % distance is  -max(cross correlation) 
                    dist = 1 - max( xc( validOverlapInds ));

                    if( dist < this.distTol )
                        sim(n,m) = dist;  %#ok<SPRIX>
                    end
                end
            end

        end

        function [ sim ] = translationInvariantHelperSparse_old( this, D, dimrange )
            
            N = size( D, 1 );
            
            if( ~exist('dimrange','var') || isempty( dimrange ))
               dimrange = 1:N; 
            end
            
            numSp = size(this.tSubPatches,1);
            mi = ( numSp - 1 )./2;
            
            M = length( dimrange );
            M2 = M*N;
            MM =  max( M2, min(10, 0.01.*M2) );
            sim = spalloc( M, N, MM );
            %MM;
            %nzmax(sim)

            % why not start m at (n+1)
            %   to compare different subpatches of the two
            for n = dimrange
                fprintf('%d\n', n);
                midPatch = D( n, this.tSubPatches(mi,:) );
                
                % size( midPatch )
                for m = 1:N
                    
                    if( m == n )
                        continue;
                    end
                    
                    for k = 1:numSp
                        dist = this.distanceFun.distance( midPatch, D(m,this.tSubPatches(k,:)) );
                        % fprintf('%d %d %d - %f\n', n, m , k, dist );
                        if( dist < this.distTol )
                            sim(n,m,k) = dist;  %#ok<SPRIX>
                        end
                    end
                end
            end
            
            fprintf('num nan: %d\n',nnz( isnan( sim )));
            
%             sim(sim == bignum) = 0;
%             aa = sum( sim, 3 )
%             sum(aa,1)
%             simOut = sparse( sim );
        end
        
        function [ sim ] = translationInvariantHelper( this, D, dimrange )
            N = size( D, 1 );
            
            numSp = size(this.tSubPatches,1);
            mi = ( numSp - 1 )./2;
            
            bignum = this.distTol .* 1000;
            sim = bignum.*ones( N, N, numSp );
            
            % why not start m at (n+1)
            %   to compare different subpatches of the two
            for n = 1:N
                fprintf('%d\n', n);
                midPatch = D( n, this.tSubPatches(mi,:) );
                
                for m = 1:N
                    
                    if( m == n )
                        continue;
                    end
                    
                    for k = 1:numSp
                        dist = this.distanceFun.distance( midPatch, D(m,this.tSubPatches(k,:)) );
%                         fprintf('%d %d %d - %f\n', n, m , k, dist );
                        if( dist < this.distTol )
                            sim(n,m,k) = dist;
                        end
                    end
                end
            end
            
            fprintf('num nan: %d\n',nnz( isnan( sim )));
            
        end
        
    end % methods

    methods( Static )
        
        function [tx,ty] = inferTranslations( patchSize, dictPatchSize )            
            diff = patchSize - dictPatchSize;
            [ tx, ty ] = ndgrid( 1:diff(1)+1, 1:diff(2)+1);
            tx = tx(:);
            ty = ty(:);
        end
        
        function [ sp, midIdx ] = translationSubPatches( patchSize, dictPatchSize )
            [tx, ty ] = Tid.inferTranslations( patchSize, dictPatchSize );
            N = length( tx );
            sp = false( N, prod( patchSize ));
            
            for n = 1:N
               tmp = false( patchSize );
               tmp( tx(n) : tx(n) + dictPatchSize(1) - 1, ...
                    ty(n) : ty(n) + dictPatchSize(2) - 1 ) = 1;
                sp( n, : ) =  reshape( tmp, 1, [] );
            end
            midIdx = (N-1)./2;
        end
        
        function dictFun = getDictionaryBuildingFunction()
            dat = rand( 10, 10 );
            
            param.K = 5;  % dictionary size
            param.lambda = 0.1;
            param.numThreads = 1; % number of threads
            param.verbose = 0;
            param.iter = 400;  % num iterations.

            try
                D = mexTrainDL( dat, param );
            catch e
                % e
                dictFun = @( X, params )( ... 
                    X( :, randperm( size(X,2), params.K)) );
                return;
            end
            dictFun = @( X, params )(mexTrainDL( X, params ));
        end
            
        function param = defaultParams()

            param.K = 250;  % dictionary size
            param.lambda=0.15;
            param.numThreads = 4; % number of threads
            param.batchsize = 100;
            param.iter = 100;
            param.verbose = 0;

        end
        
        function [xi,yi,zi] = repz( sz, f, isvec )
            if( ~exist('isvec', 'var') || isempty( isvec ))
                isvec = 0;
            end
            zr = repmat(1:sz(3),f,1);
            [xi,yi,zi] = ndgrid( 1:sz(1), 1:sz(2), zr(:));
            if( isvec )
                xi = sub2ind( sz, xi, yi, zi );
                xi = xi(:);
                clear yi zi;
            end
            
        end
            
        function s2rFun = patchSlicesToRows3d()
            s2rFun = @(X)( reshape(X, [], size(X,3) )');
        end
        
        function r2sFun = patchRowsToSlices3d()
            r2sFun = @(X, sz2d)( reshape(X', [sz2d size(X,1)]));
        end
        
        function downsampler = getDownsampler2dx()
            downsampler = @(X,f)(reshape(mean(reshape( X', f, [] )), size(X')./[f 1]))';
        end
        
        function downsampler = getDownsampler2dy()
            downsampler = @(X,f)(reshape(mean(reshape( X, f, [] )), size(X)./[f 1]));
        end
        
        function downsampler = getDownsampler2dxy()
            downsamplerX = @(X,f)(reshape(mean(reshape( X', f, [] )), size(X')./[f 1]))';
            downsamplerY = @(X,f)(reshape(mean(reshape( X, f, [] )), size(X)./[f 1]));
            downsampler = @(X,f)(downsamplerX( downsamplerY(X,f) , f));
        end
        
        function summer = sum2dx()
            summer = @(X,f)(reshape(sum(reshape( X', f, [] )), size(X')./[f 1]))';
        end
        
        function summer = sum2dy()
            summer = @(X,f)(reshape(sum(reshape( X, f, [] )), size(X)./[f 1]));
        end
        
        function summer = sum2dxy()
            sumX = Tid.sum2dx();
            sumY = Tid.sum2dy();
            summer = @(X,f)(sumX( sumY(X,f) , f));
        end
        
        function downsampler = getDownsampler3dz()
            downsampler = @(X,sz,f)(permute(reshape(mean(reshape( permute(X, [3 2 1]), f, [] )), sz([3 2 1])./[f 1 1]), [3 2 1]));
        end
        
        function downsamplerRs = getDownsampler3dzRs()
            downsamplerRs = @(X,sz,f) (reshape(permute(reshape(mean(reshape(permute( reshape(X, [sz size(X,2)] ), [ 3 2 1 4 ]), f, [])), [sz(3)./f sz(1) sz(2) size(X,2)]), [ 3 2 1 4 ]), [], size(X,2) )); 
        end
        
        function X_ds = downsampleZAvgRs( X, sz, dsFactor )
            N = size( X, 1 );
            M = prod( sz );
            Mds = M ./ dsFactor;
            X_ds = zeros( N, Mds );
            for n = 1:N
               im = reshape( X(n,:), sz );
               imds = Tid.downsampleZAvg( im, dsFactor );
               X_ds( n, : ) = imds( : );
            end
        end
        
        function im_ds = downsampleZAvg( im, dsFactor )
            sz = size( im );
            szds = sz ./ [1 1 dsFactor];
            szdsp = szds( [3 1 2 ]);
            imv = reshape(permute( im, [3 1 2] ), dsFactor, [] );
            im_ds = permute( reshape(sum( imv, 1 ), szdsp ), [ 2 3 1]);
        end
        
    end % static methods
end % Tid class
