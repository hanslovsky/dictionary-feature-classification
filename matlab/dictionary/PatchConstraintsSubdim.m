classdef PatchConstraintsSubdim < PatchConstraints
    % PatchConstraintsSubdim contains lo
    %
    % John Bogovic
    % HHMI 
    % February 2015

    properties
        pairLocRngXY;
    end
    
    methods
        
        function this = PatchConstraintsSubdim( sz2d, sz3d, f, overlappingPatches, scaleByOverlap )
            % Constructor
            this = this@PatchConstraints( sz3d(1), f, 0, scaleByOverlap, []  );
            
            this.sz2d = sz2d;
            this.sz3d = sz3d;
            
            this.iniLocs();
        end
        
        function iniLocs( this )
            
            this.pairLocRngXY = (1) : (this.sz3d(1)- this.f + 1);
            this.pairLocRng   = (1) : this.f : this.sz3d(1);
            
            [dListXY, xyList] = ndgrid( 1:2, this.pairLocRngXY );
            [dListZ ,  zList] = ndgrid(   3, this.pairLocRng   );
            
            this.dimXyzList = [ dListZ(:) ,  zList(:); ...
                                dListXY(:), xyList(:); ];
                            
            this.numLocs = size( this.dimXyzList, 1 );
            
            this.numConstraints = prod( this.sz2d ) *  this.numLocs;
        end
        
    end
    
end
