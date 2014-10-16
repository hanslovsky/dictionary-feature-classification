function [rootConstraint, leafConstraint] = rootConstraintFromArray( paramList )
% Usage:
%   rootConstraint = rootConstraintFromPatch( idxList, num, dim, downsampleFactor )

num = size( paramList, 1);

dimList = paramList(:,1);
xyzList = paramList(:,2);
idxList = paramList(:,3);

import net.imglib2.algorithms.patch.*;
import net.imglib2.algorithms.opt.astar.*;

spl = SubPatch2dLocation( dimList(1), xyzList(1), idxList(1), 0 );
rootConstraint = SortedTreeNode( spl );
           
parent = rootConstraint;

for i = 2:num
    
    newSpl  = SubPatch2dLocation( dimList(i), xyzList(i), idxList(i), 0 );
    newNode = SortedTreeNode( newSpl, parent );
    
    parent = newNode;
end

leafConstraint = newNode; 

end