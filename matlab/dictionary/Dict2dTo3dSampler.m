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
       maxIters  = 500; % max iterations for build
       convIters =  20; % # of iters at a given cost ('flatness')
                        % that indicate convergence
       convEps =  0.001; % convergence epsilon (defines 'flatness')
       
       useSubset = 0; 
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
        
        function [ patchParams, iteration, costs ] = build3dPatch( this, iniPatch, excludeParam )
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
            
            % initialize if given something
            if( exist( 'iniPatch', 'var' ) && ~isempty( iniPatch ))
                patchParams = iniPatch;
            else
                patchParams = randi( this.numDict, N, 1 );
            end
            
            converged = 0;           
            iteration = 1;
            
            % randomly generate the patch location 
            % that will be updated at each iteration
            randomCoords = randi( N, this.maxIters, 1 );
            if( isExclusion )
                keepInds = 1:this.pc.numLocs;
                keepInds = keepInds( ~exclude );
                randomCoords = keepInds( randomCoords );
            end
            
            b = this.pc.constraintValueList( this.D2d, patchParams );
            costs = -1.*ones( this.maxIters, 1 );
            
            lastCost = inf;
            itersAtThisCost = 0;
            
            while( ~converged )
                
                i = randomCoords( iteration );
                
                if( this.useSubset )
                    [bestidx, theseCosts] = this.bestPatchConfigSub( b, i );    
                else
                    [bestidx, theseCosts] = this.bestPatchConfig( b, i );
                end
                patchParams( i ) = bestidx;
                
                costs(iteration) = theseCosts( bestidx );
                
                % has the cost changed much?
                if( abs( costs(iteration) - lastCost  ) < this.convEps )
                   itersAtThisCost = itersAtThisCost + 1;
                else
                    itersAtThisCost = 0;
                end
                lastCost = costs(iteration);
                
                % converged if we've been at the same cost for awhile
                if( itersAtThisCost == this.convIters )
                    converged = 1;
                end
                
                % force exit after max iteration count
                if(iteration == this.maxIters)
                    converged = 1;
                else
                    iteration = iteration + 1;
                end
                
                b = this.pc.updateConstraints( this.D2d, b, i, bestidx );
                
            end % iteration
            
            if( this.verbose )
                fprintf('sampler - buildPatch3d converged after %d iterations\n', (iteration-1) );
            end
            
            costs = costs(1:iteration-1);
            patchParams = vecToRowCol( patchParams, 'col');
            
            patchParams = [ this.pc.dimXyzList, ...
                            patchParams ];
            
        end % build3dPatch
        
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
    end
end