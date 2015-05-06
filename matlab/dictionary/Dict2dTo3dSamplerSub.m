classdef Dict2dTo3dSamplerSub < Dict2dTo3dSampler
    % Dict2dTo3dSamplerSub
    %
    % This class contains logic that builds a high-res 3d dictionary from 
    % 2d patches at a lower resolution.  Uses PatchConstraintsSubdim.
    %
    % John Bogovic
    % HHMI
    % February 2015
    
    methods
        
        function this = Dict2dTo3dSamplerSub( D2d, sz, f, overlappingPatches, scaleByOverlap, comparator )
            % Constructor
            % D2d -
            % sz  - size of 2d patches
            % f   - downsampling factor
            this = this@Dict2dTo3dSampler( D2d, sz, f, overlappingPatches, scaleByOverlap, comparator );
            
            this.pc.buildCmtx();
        end
    end  % methods
    
    methods ( Access = protected )
        
        function iniPatchConstraints( this )
            % initialized the PatchConstraintsSubdim object.
            %   it allows patches to overlap in x and y but not in z
            %   (the assumption being that observed patches are z-projections
            %   of high-res patches)
            this.pc = PatchConstraintsSubdim( this.sz3d(1:2), this.sz3d, this.f, ...
                this.overlappingPatches, this.scaleByOverlap );
        end
    end % protected methods
    
    methods( Static )
        
        function [ K ] = chooseCutoffInflection( dictCosts, L ) 
            deriv = diff(dictCosts(1:L));
            deriv2 = diff( deriv );
%             K = find( (deriv(1:end-1) > 0) & (deriv2 < 0) , 1 );
            K = find( (deriv2 < 0) , 1 );
        end 
    end
end
