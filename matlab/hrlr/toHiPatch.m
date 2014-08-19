function out = toHiPatch( elem, row, isrot )
% out = toHiPatch( elem, row, isrot )

mxsz = max(size(elem));
out = zeros([ mxsz mxsz ] );
out(row,:) = elem;

if( isrot )
   out = out'; 
end