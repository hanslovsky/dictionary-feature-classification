function [ svmClassifier, trainPred, trainErr, ...
            testPred, testErr ] = classifyEvaluateSvm( alphaTrain, alphaTest, ...
                                    labels_train, labels_test, opts )
% Usage: 
% [ svmClassifier, trainPred, trainAcc, ...
%            testPred, testAcc ] = classifyEvaluateSvm( alphaTrain, alphaTest, ...
%                                    labels_train, labels_test, opts )
% 
% Inputs:
%   alphaTrain - training features
%   alphaTest  - testing features

svmClassifier = svmtrain(  alphaTrain, labels_train, opts{:} );

trainPred = svmclassify( svmClassifier, alphaTrain );
trainErr = nnz( trainPred - labels_train ) ./ length( labels_train );

testPred = svmclassify( svmClassifier, alphaTest );
testErr = nnz( testPred - labels_test ) ./ length( labels_test );
