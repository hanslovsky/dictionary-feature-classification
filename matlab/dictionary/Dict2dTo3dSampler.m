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
    

    properties ( SetAccess = private )
       maxIters = 50; % quadratic program options
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
        end
        
        function [ patchParams, iteration ] = build3dPatch( this, iniPatch )
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            N = size( this.dimXyzList, 1 );
            
            % initialize if given something
            if( exist( 'iniPatch', 'var' ) && ~isempty( iniPatch ))
                patchParams = iniPatch;    
            else
                patchParams = randi( this.numDict, 1, 1 );
            end
            
            converged = 0;           
            iteration = 1;
            
            % generate
            randomCoords = randi( N, this.maxIters, 1 );
            
            while( ~converged )
                
                i = randomCoords( iteration );
                
                j = pc.
                cmtx
                
                iteration = iteration + 1;
                
                % force exit after max iteration count
                if(iteration > this.maxIters)
                    converged = 1;
                end
            end % iteration
        end % build3dPatch
        
    end
end