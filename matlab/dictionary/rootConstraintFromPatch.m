function [rootConstraint, leafConstraint] = rootConstraintFromPatch( idxList, num, dim, downsampleFactor )
% Usage:
%   rootConstraint = rootConstraintFromPatch( idxList, num, dim, downsampleFactor )

import net.imglib2.algorithms.patch.*;
import net.imglib2.algorithms.opt.astar.*;

spl = SubPatch2dLocation( dim, 1, idxList(1), 0 );
rootConstraint = SortedTreeNode( spl );
           
parent = rootConstraint;

for i = 2:num
    
    xyz = 1 + (i-1) * downsampleFactor;
    newSpl  = SubPatch2dLocation( dim, xyz, idxList(i), 0 );
    newNode = SortedTreeNode( newSpl, parent );
    
    parent = newNode;
end

leafConstraint = newNode; 

end