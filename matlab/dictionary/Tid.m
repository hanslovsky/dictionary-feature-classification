classdef Tid < handle
    % Tid - transformationally invariant dictionary
    %
    % John Bogovic
    % HHMI
    % August 2014

    properties ( SetAccess = private )

        D;      % The dictionary elements
        Dxfm;   % The full dictionary

        params; % dictionary learning parameters
        bigIters = 1;
        dictDist = 'euclidean';
        distTol  = 0.1;

        X;      % The data
        alpha;  % The dictionary coefficients

        numSamples;
        numVars;

        patchSize;

        rotRng;
        rotXfmIdx;

    end

    methods

        % Constructor
        function this = Tid( X, patchSize, params, distTol )

            if( ~exist( 'distTol', 'var') || isempty(distTol))
                this.distTol = 0.1;    
            else
                this.distTol = distTol;
            end
                
            this.X = X;  % the data
            [this.numVars, this.numSamples] = size( X );

            % check that patchSize is compatible with dictionary size
            if ( prod( patchSize ) ~= this.numVars )
                error( 'patch size not compatible with dictionary size');
            end

            this.patchSize = patchSize;

            if( ~exist('params','var') || isempty( params ))
                this.params = Tid.defaultParams();
            else 
                this.params = params;
            end

            genVectorTransformations( this );

        end % constructor


        % Very experimental
        function D = buildDictionary(this)
            
            try
                this.D = mexTrainDL( this.X, this.params );
            catch e
                disp('no mexTrainDL - choosing dictionary by random sampling');
                %e % print error 
                this.D = this.X( :, randperm( this.numSamples, this.params.K));
            end

            for iter = 1:this.bigIters
               
                % initialize with previous dictionary
                this.params.D = this.D;

                D = this.D;
            end

        end % the magic

        function alpha = getFeatures( this )
            try
                this.D = mexTrainDL( this.X, this.params );
                alpha = mexLasso( this.X, this.D, this.params );
            catch e
                disp('no mexLasso - determining features ');
                %e % print error 
                alpha = randn( size(this.X,2), size(this.D,2) );
            end
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
        end


        function genVectorTransformations( this )
            % possible transformations include

            % comupte the index look-ups for this transformation 
            % dont reshape them.
            
            % check if we're working in 2d or 3d
            patchSzTmp = this.patchSize( this.patchSize > 1);
            ndim = length( patchSzTmp );
            
            if( ndim ==3 )
                this.rotXfmIdx = xfmToIdx( cubeSymmetry(), patchSzTmp, 0 );
            elseif (ndim == 2)
                this.rotXfmIdx = xfmToIdx( squareSymmetry(), patchSzTmp, 0 );
            end
            
            if( ~isempty( this.rotRng ) )
                this.rotXfmIdx = this.rotXfmIdx( this.rotRng );
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
            param.verbose=1;

        end

    end % static methods
end % Tid class
