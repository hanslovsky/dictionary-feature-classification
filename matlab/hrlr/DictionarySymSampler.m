classdef DictionarySymSampler < handle
    % DictionarySymSampler Class samples 
    %
    % John Bogovic
    % HHMI
    % August 2014

    properties (SetAccess = private )

        D;  % The dictionary elements:
            % rows are variables, 
            % columns are observations
            
        numDictElems;   % num dictionary elements
        numVars;        % num variables 

        patchSize;  % patch size
        reshapeFunction = @reshape;
        %reshapeFunction = @reshapePy;

        dictProbs;  % probability of observing

        rotRng;     % range of rotations
        shiftRng;   % range of shifts 
        rotProbs;   % probabilities of possible rotations
        shiftProbs; %  probabilities of possible shifts

        rotXfm;
        flipXXfm;
        flipZXfm;

    end

%   % NONE SO FAR
%   properties ( SetAccess = public )
% 
%   end

    methods

        % Constructor
        function this = DictionarySymSampler( D, dictProb, patchSize, ...
                rotRng, robProb, ...
                shiftRng, shiftProb  )

            this.D = D;
            this.patchSize = vecToRowCol( patchSize, 'row' );

            this.numVars      = size(D,1); 
            this.numDictElems = size(D,2); 

            % check that patchSize is compatible with dictionary size
            if ( prod( patchSize ) ~= this.numVars )
                error( 'patch size not compatible with dictionary size');
            end

            % default uniform probability over dictionary elements
            if( ~exist('prob','var') || isempty( prob ))
                this.dictProbs = ones( this.numDictElems, 1) ./ this.numDictElems;
            end

            % Fill in default values
            if( ~exist('rotRng','var') || isempty( rotRng ))
                rotRng = [ 90 180 270 ];
            else
                this.rotRng = rotRng;
            end
            
            if( ~exist('rotProb','var') || isempty( rotProb ))
                this.rotProbs = ones( length(rotRng), 1) ./ length(rotRng);
            else
                this.rotProbs = rotProb;
            end
            
            if( ~exist('shiftRng','var') || isempty( shiftRng ))
                this.shiftRng = [];
            else
                this.shiftRng = shiftRng;
            end
            
            if( ~exist('shiftProb','var') || isempty( shiftProb ))
                this.shiftProbs = ones( length(shiftRng), 1) ./ length(shiftRng);
            else
                this.shiftProbs = shiftProb;
            end

        end % Constructor
        
        % Returns N observations under the distribution described by
        % this object
        function X = sample(this,N)

            elemSamps  = logical(mnrnd( 1, this.dictProbs, N ));
            rngSamps   = logical(mnrnd( 1,   this.rngProb, N ));
            shiftSamps = logical(mnrnd( 1, this.shiftProb, N ));

            m = [1:N]';
            elemSamps  = elemSamps  * m;

            X = xfmPatches( this, D(:, elemSamps), rotSamps, flipSamps );
                            
        end

        % Transform the patch
        function X_xfm = xfmPatches( this, X, i, xfmj )
            X_xfm = X;
            
        end
    
        function genVectorTransformations( this )
            % possible transformations include
            this.rotXfm = xfmToIdx( cubeSymmetry(), this.patchSize );

        end

        function genAtomicVectorTransformations( this )

            imr = 1:this.numVars;
            im  = this.reshapeFunction( imr, this.patchSize );

            this.rotXfm{1} = permute( im( :, end:-1:1, : ), [2 1 3]); 
            this.rotXfm{1} = this.rotXfm{1}(:);

            this.rotXfm{2} = im( end:-1:1, end:-1:1, : ); 
            this.rotXfm{2} = this.rotXfm{2}(:);

            this.rotXfm{3} = permute( im( end:-1:1, :, : ), [2 1 3]); 
            this.rotXfm{3} = this.rotXfm{3}(:);

            this.flipXXfm = im( end:-1:1, :, : );
            this.flipXXfm = this.flipXXfm(:);

            this.flipZXfm = im( :, :, end:-1:1 );
            this.flipZXfm = this.flipZXfm(:);

        end

    end % Methods

    methods( Static )

        function D = standardDictionary1()

            corner = [ 0 0 0; 0 1 0; 0 1 1 ];
            edgeM  = [ 0 0 0; 1 1 1; 0 0 0 ];
            edgeE  = [ 1 1 1; 0 0 0; 0 0 0 ];
            edgeD  = [ 1 0 0; 0 1 0; 0 0 1 ];
            blob   = [ 1 1 0; 1 1 0; 0 0 0 ];

        end

    end

end
