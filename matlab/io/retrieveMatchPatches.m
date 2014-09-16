function [patches, coords, isPos] = retrieveMatchPatches( f, im, patchSize )
% function p = retrieveMatchPatches( f, im, patchSize )

csvdat = csvread( f, 1, 0 );
pts = reshape( csvdat', 3, [] )';

isPos = ones( size(pts,1),1);
isPos(2:2:end) = 0;

coords = { pts(:,1), pts(:,2), pts(:,3) };
[patches, coords] = grabPatchesSimple( im, patchSize, [], coords, [] );