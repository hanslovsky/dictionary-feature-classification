classdef Tid < handle
    % Tid - transformationally invariant dictionary
    %
    % John Bogovic
    % HHMI
    % August 2014

    properties ( SetAccess = private )

        D;      % The dictionary elements
        params; % dictionary learning parameters
        bigIters = 1;
        dictDist = 'euclidean';
        distTol  = 0.1;

        X; % The data

        numSamples;
        numVars;

        patchSize;

        rotRng;
        rotXfmIdx;

    end

    methods

        % Constructor
        function this = Tid( X, patchSize, params )

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

        end % constructor


        % Very experimental
        function D = workSomeMagic(this)

            this.D = mexTrainDL( this.X, this.params );

            for iter = 1:bigIters
               
                % initialize with previous dictionary
                this.params.D = this.D;

                D = this.D;
            end

        end % the magic

        % Make the dictionary rotationally invariant
        function makeDictRotInv( )

            K = this.params.K;

            numInvalid = 0;
            invalid =  zeros(1,K);

            for i = 1:K

                % check if this index is still valid
                % it may have been removed at a previous iteration
                if ( nnz(invalid == i) > 0 )
                    continue;
                end

                Dxfm = D(:, i ); 
                Draw = this.D(:, setdiff( valid, i) );

                pd = pdist2( Dxfm', Draw' );

                % if this is invalid
                if ( 0 )
                    numInvalid = numInvalid + 1; 
                    invalid( numInvalid ) = i;
                end

            end


        end % make dict rot inv

        function Dxfm = xfmDict( d )

        end

        function fillInDict( ) 

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
