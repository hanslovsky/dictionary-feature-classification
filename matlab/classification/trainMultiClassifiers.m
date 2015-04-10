function [ classifiers ] = trainMultiClassifiers( X, Y, classifierFun, classifierArgs )
% Expects classifierFun to take arguments in a format like
%   classifierFun( X, Y, <args> );
%
% X - the features (rows are observations)
% Y - the outputs  (rows are observations)

nOut = size( Y, 2 );
classifiers = cell( nOut, 1 );

for i = 1:nOut 
    
    classifiers{i} = classifierFun( X, Y(:,i), classifierArgs{:});
    
end
