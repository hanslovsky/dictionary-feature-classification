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
        maxshift;

        rotXfmIdx;
        shiftXfmIdx;

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

            [ this.numVars, this.numDictElems ] = size(D); 

            % check that patchSize is compatible with dictionary size
            if ( prod( patchSize ) ~= this.numVars )
                error( 'patch size not compatible with dictionary size');
            end

            % TODO check that patch size is compatible with shift
            % magnitude

            % default uniform probability over dictionary elements
            if( ~exist('dictProb','var') || isempty( dictProb ))
                this.dictProbs = ones( this.numDictElems, 1) ./ this.numDictElems;
            else
                this.dictProbs = dictProb;
            end

            % Fill in default values
            if( ~exist('rotRng','var') || isempty( rotRng ))
                %rotRng = [ 90 180 270 ];
                [mi,mj,mk] = meshgrid( 0:1, 0:1, 0:11);
                v = 32*mi + 16*mj + mk;
                rotRng = sort(v(:));
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
                this.maxshift = zeros( size( patchSize ));
            else
            %    this.shiftRng = vecToRowCol(shiftRng, 'col');
                 this.shiftRng = shiftRng;
                 this.maxshift = max( abs(shiftRng) );
            end
            
            if( ~exist('shiftProb','var') || isempty( shiftProb ))
                this.shiftProbs = ones( length(shiftRng), 1) ./ length(shiftRng);
            else
                this.shiftProbs = shiftProb;
            end

            genVectorTransformations( this );

        end % Constructor
        
        % Returns N observations under the distribution described by
        % this object
        function X = sample(this,N)

            elemSamps  = logical(mnrnd( 1, this.dictProbs,  N ));
            rotSamps   = logical(mnrnd( 1,  this.rotProbs,  N ));

            if ( isempty ( this.shiftProbs ))
                shiftSamps = -1;
            else
                shiftSamps = logical(mnrnd( 1, this.shiftProbs, N ));
            end

            m = [ 1:this.numDictElems ]';
            elemSamps  = elemSamps  * m;

            X = xfmPatches( this, this.D(:, elemSamps), rotSamps, shiftSamps );
                            
        end

        % Transform the patch
        function X_xfm = xfmPatches( this, X, rotIdx, shiftIdx )
            X_xfm = X;
           
            % shift first, then rotate
            % TODO - make this work for arbitrary dictionary / patch
            N = size( X, 2 );
            for i = 1:N 

                patchTmp = padarray( reshape (  X(:,i), this.patchSize ), ...
                                this.maxshift );

                if( shiftIdx < 1 )
                    X_xfm(:,i) = reshape( patchTmp(this.rotXfmIdx{ rotIdx(i,:)}), [], 1);
                else
                    X_xfm(:,i) = reshape( patchTmp(this.shiftXfmIdx{shiftIdx(i,:)}(this.rotXfmIdx{ rotIdx(i,:)})), [], 1);
                end

            end
            
            % X_xfm(:,i) = patchTmp( this.shiftXfmIdx{ shiftIdx(1,:) }(this.rotXfmIdx) );

        end
    
        function genVectorTransformations( this )
            % possible transformations include

            this.rotXfmIdx = xfmToIdx( cubeSymmetry(), this.patchSize );

            if( ~isempty( this.rotRng ) )
                this.rotXfmIdx = this.rotXfmIdx( this.rotRng );
            end

            if( ~isempty( this.shiftRng ))
                this.shiftXfmIdx = shiftToIdx( this.shiftRng, 2.*this.maxshift + this.patchSize );
            end
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

        function D = standardDictionary2d()
            D = zeros( 9, 4 );

            D(:,1) = reshape([ 0 0 0; 0 1 0; 0 1 1 ], [], 1); % corner
            D(:,2) = reshape([ 0 0 0; 1 1 1; 0 0 0 ], [], 1); % edgeM
            D(:,3) = reshape([ 1 1 1; 0 0 0; 0 0 0 ], [], 1); % edgeE
            D(:,4) = reshape([ 1 1 0; 1 1 0; 0 0 0 ], [], 1); % blob

            edgeD  = [ 1 0 0; 0 1 0; 0 0 1 ]; % edgeD
            
        end
        function D = simpleDictionary3d()
            corner = cat(3, [ 0 0 0; 0 1 1; 0 1 0 ], ...
                            [ 0 0 0; 0 1 1; 0 1 0 ], ...
                            [ 0 0 0; 0 1 1; 0 1 0 ] );

            blob = cat(3, [ 0 0 0; 0 0 0; 0 0 0 ], ...
                          [ 0 0 0; 0 1 1; 0 1 1 ], ...
                          [ 0 0 0; 0 1 1; 0 1 1 ] );

            edgeM = cat(3, [ 0 0 0; 0 0 0; 0 0 0 ], ...
                           [ 1 1 1; 1 1 1; 1 1 1 ], ...
                           [ 0 0 0; 0 0 0; 0 0 0 ] );

            edgeE = cat(3, [ 1 1 1; 1 1 1; 1 1 1 ], ...
                           [ 0 0 0; 0 0 0; 0 0 0 ], ...
                           [ 0 0 0; 0 0 0; 0 0 0 ] );

            D = zeros( 27, 4 );

            D(:,1) = corner(:); 
            D(:,2) = edgeM(:); 
            D(:,3) = edgeE(:); 
            D(:,4) = blob(:);
            
        end

    end

end
