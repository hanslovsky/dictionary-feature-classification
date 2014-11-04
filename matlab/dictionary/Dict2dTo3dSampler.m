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
       maxIters  = 1000; % max iterations for build
       convIters =   20; % # of iters at a given cost
                         % that
       convEps =  0.001; % convergence epsilon
       
       useSubset = 1; 
    end
    
    properties ( SetAccess = private )
       pc;
    end
    
    methods
        
        % Constructor
        % D2d - 
        % sz  - size of 2d patches
        % f   - downsampling factor
        function this = Dict2dTo3dSampler( D2d, sz, f )
            this = this@Dict2dTo3d( D2d, sz, f );
            
            % use a PatchConstraints object to compute
            % the constraint matrix once 
            this.pc = PatchConstraints( sz, f );
            this.pc.buildCmtx();
            this.pc.compXsectInverses();
        end
        
        function [ patchParams, iteration, costs ] = build3dPatch( this, iniPatch, num2exclude )
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            N = size( this.dimXyzList, 1 );
            
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
        
        function iniParams = iniParamsDist( this, N )
            M = size(this.dimXyzList,1);
            iniParams = mat2cell( ...
                            randi(  this.numDict, N, M), ...
                        ones( N, 1 ));
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