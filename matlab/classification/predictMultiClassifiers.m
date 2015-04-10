function [ Y,S ] = predictMultiClassifiers( X, classifierList )
% Usage:
%   [ classifiers ] = predictMultiClassifiers( X, classifierList )

nSamples = size(X,1);
nOutputs = length( classifierList );

Y = zeros( nSamples, nOutputs );
S = zeros( nSamples, nOutputs );

for i = 1:nOutputs 
    [Yi,Si]= classifierList{i}.predict( X );
    
    if( i == 1 )
        S = zeros( nSamples, nOutputs, size(Si,2) );
    end
    Y(:,i) = Yi;
    S(:,i,:) = Si;
end
