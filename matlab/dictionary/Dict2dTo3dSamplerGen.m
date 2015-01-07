classdef Dict2dTo3dSamplerGen < Dict2dTo3dSampler
    % Dict2dTo3dSamplerGen
    %
    % This class contains logic that builds a high-res 3d dictionary from 
    % 2d patches at a lower resolution using a sampling technique rather
    % than a Astar-like search.
    %
    % John Bogovic
    % HHMI
    % January 2015
    
    properties
        tmp = [];
    end
    
    methods
        
        function this = Dict2dTo3dSamplerGen( D2d, sz2d, sz3d, f, overlappingPatches, scaleByOverlap, doInvCmtx )
            
            if( ~exist('overlappingPatches','var') || isempty(overlappingPatches))
                overlappingPatches = 1;
            end
            if( ~exist('scaleByOverlap','var') || isempty(scaleByOverlap))
                scaleByOverlap = 0;
            end
            if( ~exist('doInvCmtx','var') || isempty(doInvCmtx))
                doInvCmtx = 1;
            end
            
            this = this@Dict2dTo3dSampler( D2d, sz2d, f, overlappingPatches, scaleByOverlap );
            this.sz3d = sz3d;
            this.iniPatchConstraints( doInvCmtx );
            
            if( this.scaleByOverlap )
                this.paramScales = this.pc.overlapFraction;
            end
                
        end
        
        function obj = copy(this)
            obj = Dict2dTo3dSampler( this.D2d, this.sz2d(1), this.f, ...
                                     this.overlappingPatches, this.scaleByOverlap );
            obj.sz3d = this.sz3d;
            
            obj.D3d = this.D3d;
            obj.numDict3d = this.numDict3d;
            
            obj.clone( this );
        end

    end
    
    methods ( Access = protected )
        function iniPatchConstraints( this, doInvCmtx )
            
%             this = PatchConstraintsGeneralZ( sz2d, sz3d, f, overlappingPatches, scaleByOverlap  )
            
            this.sz2d
            this.sz3d
            
%             if( isequal( this.sz2d, this.sz3d(1:2)))
%                 this.pc = PatchConstraints( this.sz3d(1), this.f, this.overlappingPatches, this.scaleByOverlap );
%             else
                this.pc = PatchConstraintsGeneralZ( this.sz2d, this.sz3d, this.f, this.overlappingPatches, this.scaleByOverlap );
%             end
            if( doInvCmtx )
                fprintf('Computing constraint matrix inverse...');
                this.pc.buildCmtx();
                fprintf('.done\n');
            end
        end
    end

end
