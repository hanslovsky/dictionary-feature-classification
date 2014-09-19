classdef Tid_from2D < handle
    % Tid_from2D - transformationally invariant dictionary
    %
    % John Bogovic
    % HHMI
    % August 2014

    properties ( SetAccess = private )


        
        iMergeDxfms;    % indicates which original dictionary element (D)
                        % the elements of Dxfm came from
        mergeFun = @max;
        
        params; % dictionary learning parameters
        bigIters = 1;
        dictDist = 'euclidean';
        distTol  = 0.1;

        alpha;  % The dictionary coefficients

        numSamples;
        numVars;

        patchSize2d;
        patchSize3d;

        rotRng2d;
        rotXfmIdx2d;

        rotRng3d;
        rotXfmIdx3d;

    end

    properties
        X;      % The data
        
        D;      % The dictionary elements
        Dxfm;   % The full dictionary
        
    end

    methods

        function this = Tid_from2D( D, patchSize, params, distTol )
            if( ~exist( 'distTol', 'var') || isempty(distTol))
                this.distTol = 0.1;
            else
                this.distTol = distTol;
            end
            
            if( ~exist('params','var') || isempty( params ))
                this.params = Tid.defaultParams();
            else 
                this.params = params;
            end
            
            this.patchSize2d = [patchSize patchSize];
            this.patchSize3d = [patchSize patchSize patchSize];
            
%             % check that patchSize is compatible with dictionary size
%             if ( prod( patchSize ) ~= this.numVars )
%                 error( 'patch size not compatible with dictionary size');
%             end

            this.D = D;

        end % constructor
        
%         % Constructor
%         function this = Tid( X, patchSize, params, distTol )
%                 
%             this.X = X;  % the data
%             [this.numVars, this.numSamples] = size( X );
% 
%             % check that patchSize is compatible with dictionary size
%             if ( prod( patchSize ) ~= this.numVars )
%                 error( 'patch size not compatible with dictionary size');
%             end
% 
%             this.patchSize = patchSize;
% 
%             if( ~exist('params','var') || isempty( params ))
%                 this.params = Tid.defaultParams();
%             else 
%                 this.params = params;
%             end
% 
%             genVectorTransformations( this );
% 
%         end % constructor

%         function thecopy = copy( this )
%             thecopy = Tid( this.X, this.patchSize, this.params, this.distTol );
%             thecopy.X = [];
%         end

        % Very experimental
        function D = buildDictionary(this)
            
            try
                this.D = mexTrainDL( this.X, this.params );
%                 this.makeDictRotInv();
            catch e
                disp('no mexTrainDL - choosing dictionary by random sampling');
                %e % print error 
                this.D = this.X( :, randperm( this.numSamples, this.params.K));
                this.makeDictRotInv();
            end

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
        function makeDictRotInv( this, dim )

            dxsz  = size( this.D, 2); % size of transformed dictionary
            
            if( dim == 2 )
                rxi= cell2mat(this.rotXfmIdx2d);
            elseif( dim == 3)
                rxi= cell2mat(this.rotXfmIdx3d);
            end
            
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
                
                pd = pdist2( thisXfmDict', Draw', this.dictDist );
                size( pd )
                pd = min( pd, [], 1 );
                
                
                similar = ( pd < this.distTol );
                
                fprintf('excluding %d elems at %d\n', nnz(similar), i );
                
                % update indices of dict elems that are
                % transformations of this dict element
                invariantIdxs( allIdxs(j(similar)) ) = false;
                
            end

            this.D = this.D( :, invariantIdxs );
            
        end % make dict rot inv

        function Dxfm = xfmDict( this, dim )
            
            if( dim == 2 )
                rxi= this.rotXfmIdx2d;
            elseif( dim == 3)
                rxi= this.rotXfmIdx3d;
            end
            
            [r,N] = size( this.D );
            M = length( rxix );
            Dxfm = zeros( r, N*M );
            j = reshape( 1:(N*M), M, [] );
            k =  cell2mat( rxi );
            
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
            patchSzTmp2d = this.patchSize2d( this.patchSize2d > 1);
            patchSzTmp3d = this.patchSize3d( this.patchSize3d > 1);
            
            this.rotXfmIdx2d = xfmToIdx( squareSymmetry(), patchSzTmp2d, 0 );
            this.rotXfmIdx3d = xfmToIdx( cubeSymmetry(), patchSzTmp3d, 0 );
            
            if( ~isempty( this.rotRng2d ) )
                this.rotXfmIdx2d = this.rotXfmIdx2d( this.rotRng2d );
            end
            if( ~isempty( this.rotRng3d ) )
                this.rotXfmIdx3d = this.rotXfmIdx3d( this.rotRng3d );
            end
        end

    end % methods

    methods( Static )

        function param = defaultParams()

            param.K = 250;  % dictionary size
            param.lambda=0.15;
            param.numThreads = 4; % number of threads
            param.batchsize = 100;
            param.iter = 100;
            param.verbose = 0;

        end
        
        function downsampler = getDownsampler3dz()
            downsampler = @(X,f)(permute(reshape(mean(reshape( permute(x, [3 2 1]), f, [] )), sz([3 2 1])./[f 1 1]), [3 2 1]));
        end
        
        function downsamplerRs = getDownsampler3dzRs()
            downsamplerRs = @(X,sz,f) (reshape(permute(reshape(mean(reshape(permute( reshape(X, [sz size(X,2)] ), [ 3 2 1 4]), f, [])), [sz(3)./f sz(1) sz(2) size(X,2)]), [ 3 2 1 4]), [], size(X,2) )); 
        end
        
    end % static methods
end % Tid class
