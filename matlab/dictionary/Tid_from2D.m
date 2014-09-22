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
        
        D2d;      % The dictionary elements
        D3d;      % The dictionary elements
        Dxfm;     % The full dictionary
        
        dictDist = 'euclidean';
    end

    methods

        function this = Tid_from2D( X, patchSize, params, distTol )
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

            this.X = X;

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
        function buildDictionary(this)
            
            this.genVectorTransformations();
            
%             try
%                 this.D2d = mexTrainDL( this.X, this.params );
% %                 this.makeDictRotInv();
%             catch e
%                 disp('no mexTrainDL - choosing dictionary by random sampling');
                %e % print error 
                this.D2d = this.X;
                this.D2d
                this.makeDictRotInv( 2 );
                
                this.dict2dTo3d();
                this.makeDictRotInv_approx( 3 );
%             end

%             for iter = 1:this.bigIters
%                
%                 % initialize with previous dictionary
%                 this.params.D = this.D2d;
% 
%                 D = this.D2d;
% %                 this.makeDictRotInv();
%             end

            this.Dxfm = this.xfmDict( 3 );

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

        function dict2dTo3d( this, dorep )
            
            if( ~exist( 'dorep', 'var') || isempty( dorep ))
               dorep = 1; 
            end
            
            N = size( this.D2d, 1 );
            
            i = reshape( 1:prod(this.patchSize2d), this.patchSize2d );
            irepz = repmat( i, [ 1 1 this.patchSize3d(3) ]);
            irot = repmat( permute( i, [ 3 1 2]), [this.patchSize3d(3) 1 1 ]);
            
            jrep = irepz(:);
            jrot = irot(:);
            
            if( dorep )
                [ii,jj] = ndgrid( 1:N, 1:N );
                ii = ii(:);
                jj = jj(:);
            else
                ii = zeros( N.*(N-1)./2, 1);
                jj = zeros( N.*(N-1)./2, 1);
                n = 1;
                for a = 1:N
                    for b = a+1:N
                        ii(n) = a;
                        jj(n) = b;
                        n = n + 1;
                    end
                end
            end
            
%             this.D3d =  this.D2d(ii,jrep) .*  this.D2d(jj,jrot);
            this.D3d =  this.D2d(jrep,ii) .*  this.D2d(jrot,jj);
            
        end
        
        function makeDictRotInv( this, dim )

            if( dim == 2 )
                rxi = cell2mat(this.rotXfmIdx2d);
                D = this.D2d;
            elseif( dim == 3)
                rxi = cell2mat(this.rotXfmIdx3d);
                D = this.D3d;
            end
            rxi = rxi';

            nD   = size( D, 2 );
            nRxi = size( rxi, 1 );

            Dbig = zeros( size(D,1), nD.*nRxi);

            k = 1;
            for i = 1:nD
                dd = D(:,i);
                dd = dd(rxi);

                %pd = pdist( dd );
                %[n,m] = ind2subTri( size(dd,1), find( pd < 0.1 ));
                
                [~, ki ] = dissimilarVecs( dd, 0.1 );
                num2add = length( ki );

                Dbig( :, k:k+num2add-1 ) = dd(ki,:)'; 

                k = k + num2add;
            end

            Dbig = Dbig( :, 1:(k-1) );
            Dbig
            if( dim == 2 )
                this.D2d = Dbig;
            elseif( dim == 3)
                this.D3d = Dbig;
            end
            
        end

        % Make the dictionary rotationally invariant
        function makeDictRotInv_approx( this, dim )

            if( dim == 2 )
                rxi= cell2mat(this.rotXfmIdx2d);
                D = this.D2d;
            elseif( dim == 3 )
                rxi= cell2mat(this.rotXfmIdx3d);
                D = this.D3d;
            end
            
            dxsz  = size( D, 2 ); % size of transformed dictionary
            
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
                Draw = D( :, j );
                
                thisXfmDict = D( :, i );
                thisXfmDict = thisXfmDict( rxi );
                
                pd = pdist2( thisXfmDict', Draw', this.dictDist );
                size( pd )
                pd = min( pd, [], 1 );
                
                similar = ( pd < this.distTol );
                fprintf('excluding %d elems at %d\n', nnz(similar), i );
                
                % update indices of dict elems that are
                % transformations of this dict element
                invariantIdxs( allIdxs(j(similar)) ) = false;
                
            end

            if( dim == 2 )
                this.D2d = this.D2d( :, invariantIdxs );
            elseif( dim == 3)
                this.D3d = this.D3d( :, invariantIdxs );
            end
            
        end % make dict rot inv

        function Dxfm = xfmDict( this, dim )
            
            if( dim == 2 )
                rxi= this.rotXfmIdx2d;
                D = this.D2d;
            elseif( dim == 3)
                rxi= this.rotXfmIdx3d;
                D = this.D3d;
            end
            
            [r,N] = size( D );
            M = length( rxi );
            Dxfm = zeros( r, N*M );
            j = reshape( 1:(N*M), M, [] );
            k =  cell2mat( rxi );
            
            for i = 1:N
                tmp = D(:,i);
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
