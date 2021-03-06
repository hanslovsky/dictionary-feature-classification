%% multiHypoConcept
% run_script('multiHypoConcept', '20k patches, context-full');
% dbstop if error; run_script('multiHypoConcept', '20k patches, context-prob');
% dbstop if error; run_script('multiHypoConcept', 'save data');

%%
global SAVEPATH
global SAVEPREFIX

%% helper functions

downsamplerRs = Tid.getDownsampler3dzRs();
normRows = @(x)( x./repmat( sqrt(sum( x.*x, 2 )), 1, size(x,2)) );

classifySVM = @( X, Y, opts )( svmtrain(  X, Y, opts{:} ) );
evaluateSVM = @( svm, X )( svmclassify( svm, X ) );

classifyRF = @( X, Y, numTrees, opts )( TreeBagger(numTrees, X, Y, opts{:}) );
evaluateRF = @( rf, X )( cellfun( @str2double, rf.predict(X)) );

%%

dataLoc = '';
dataTestLoc = '';

destdir='/data-ssd1/john/projects/dictionary/multiHypo';

%% parameters

% for reproducibility
rng( 12 );

N = 20000;
patchSize  = [9 9 9];

dsfactor = 3;
patchSizeD = patchSize ./ [1 1 dsfactor];
numNeighbors = 27;

excludeNeighbors = 1;
poolingType = '';

donorm = 0;
distMeasure = 'correlation';

params = Tid.defaultParams();
params.K = 500;
params.iter = 200;
params.batchsize = 100;
params.verbose = 1;

classify = @( X, Y, opts )( classifyRF( X, Y, 200, opts) );
evaluate = evaluateRF;
classopts = {};

% classify = @( X, Y, opts )( classifySVM( X, Y, opts ) );
% evaluate = evaluateSVM;
% classopts = {'kernel_function', 'rbf', 'options', struct('MaxIter', 15000)};

%% load data

training = 18;
test = 3;
ds = ds_medulla(training, test);

im = read_image_stack( ds.data_fn{1} );
lb = read_image_stack( ds.labels_fn{1} );
lb = (mean( lb , 4 ) > 0.5);
lb = uint8 ( lb );

msk = read_image_stack( ds.mask_fn{1} );
msk = max( msk, [], 4 );

%% grab patches
nPos = N./2;
nNeg = N - nPos; 

% training
[ppos,~]  = grabPatchesNeighborhood( im, patchSize, nPos, [], (msk &  lb) );
[pneg,~]  = grabPatchesNeighborhood( im, patchSize, nNeg, [], (msk & ~lb) );

% testing
[pposTest,~]  = grabPatchesNeighborhood( im, patchSize, nPos, [], (msk &  lb) );
[pnegTest,~]  = grabPatchesNeighborhood( im, patchSize, nNeg, [], (msk & ~lb) );

clear im msk lb;

%% build dictionary

if( excludeNeighbors )
    p = [ ppos(1:numNeighbors:end,:) ; ...
          pneg(1:numNeighbors:end,:) ];
      
    Y = [ ones( size(ppos,1)./numNeighbors, 1); ...
          zeros( size(pneg,1)./numNeighbors, 1) ];
          
    pTest = [ pposTest(1:numNeighbors:end,:) ; ...
              pnegTest(1:numNeighbors:end,:) ];          
    YTest = [ ones( size(pposTest,1)./numNeighbors, 1); ...
              zeros( size(pnegTest,1)./numNeighbors, 1) ];
else
    p = [ ppos; pneg ];
    Y = [ ones( size(ppos,1), 1); ...
          zeros( size(pneg,1), 1) ];
      
    pTest = [ pposTest; pnegTest ];
	YTest = [ ones( size(pposTest, 1), 1); ...
              zeros( size(pnegTest,1), 1) ];
end

% clear ppos pneg;

D = mexTrainDL( p', params );

%% 

X = mexLasso( p', D, params );
X = X';

XTest = mexLasso( pTest', D, params );
XTest = XTest';

if( isempty(dataLoc))
    save( sprintf('%s/%s_feats_trn.m', SAVEPATH, SAVEPREFIX ), X);
    save( sprintf('%s/%s_dict.m', SAVEPATH, SAVEPREFIX ), D);
end
if( isempty(dataTestLoc))
    save( sprintf('%s/%s_feats_trn.m', SAVEPATH, SAVEPREFIX ), XTest);
end

% return;

%% quick analysis
% 
% dictPairDist = pdist(  D', distMeasure );
% sf = squareform( dictPairDist );
% 
% [~,k]=sort( sf(:) );
% [i,j] = ind2sub( size(sf), k);
% 
% 
% for n = 1:length(k)
%     
%     m = k(n);
%     i(m)
%     j(m)
%     sf(i(m), j(m))
%     
% 	figure;
%     imdisp( permute(reshape( D(:,i(m)), patchSize ), [1 2 4 3]), 'border', 0.1 );
%     
%     figure;
%     imdisp( permute(reshape( D(:,j(m)), patchSize ), [1 2 4 3]), 'border', 0.1 );
%     
%     pause;
%     close all;
%     
% end

%% classification level 1
% X = full( X );
% XTest = full( XTest );
% 
% classifier = classify( X, Y, classopts );
% 
% pred       = evaluate( classifier, X );
% trainingAcc = (length(Y) - nnz(Y - pred))./(length(Y));
% 
% predTest       = evaluate( classifier, XTest );
% testAcc = (length(YTest) - nnz(YTest - predTest))./(length(YTest));
% 
% fprintf( '(1) training accuracy: %f\n', trainingAcc );
% fprintf( '(1) test     accuracy: %f\n', testAcc );

% clear X XTest p pTest

%% setup for classification level 2

p2 = [ ppos; pneg ];
Y2 = [ ones( size(ppos,1), 1); ...
       zeros( size(pneg,1), 1) ];

pTest2 = [ pposTest; pnegTest ];
YTest2 = [ ones( size(pposTest, 1), 1); ...
           zeros( size(pnegTest,1), 1) ];

% clear ppos pneg
% clear pposTest pnegTest

X2 = mexLasso( p2', D, params );
X2 = X2';

XTest2 = mexLasso( pTest2', D, params );
XTest2 = XTest2';

if( isempty(dataLoc))
    save( sprintf('%s/%s_featsNbr_trn.m', SAVEPATH, SAVEPREFIX ), X2);
end
if( isempty(dataTestLoc))
    save( sprintf('%s/%s_featsNbr_trn.m', SAVEPATH, SAVEPREFIX ), XTest2);
end

return;

X2 = full( X2 );
XTest2 = full( XTest2 );

%% make predictions on all neighborhood patches

[pred_1c2, prob_1c2 ]        = classifier.predict( X2 );
[predTest_1c2, probTest_1c2] = classifier.predict( XTest2 );


%% reshape and build final feature vector

% X2     = reshape( X2', size(X2,2).*numNeighbors, []  )';
% XTest2 = reshape( XTest2', size(XTest2,2).*numNeighbors, [] )';

% X2 = [ X2(1:27:end,:), reshape(prob_1c2(:,2)', numNeighbors, [])'];
% XTest2 = [XTest2(1:27:end,:), reshape(probTest_1c2(:,2)', numNeighbors, [])'];

X2 = reshape(prob_1c2(:,2)', numNeighbors, [])';
XTest2 = reshape(probTest_1c2(:,2)', numNeighbors, [])';

size( X2 )
size( XTest2 )

Y2     = Y2(1:numNeighbors:end);
YTest2 = YTest2(1:numNeighbors:end);

% pool if requested
if( ~isempty( poolingType ))
   X2 = poolFeatures( poolingType, X2, numNeighbors );
   XTest2 = poolFeatures( poolingType, XTest2, numNeighbors );
end


% classification level 2

classifier2 = classify( X2, Y2, classopts );

pred2       = evaluate( classifier2, X2 );
trainingAcc2 = (length(Y2) - nnz(Y2 - pred2))./(length(Y2));

predTest2       = evaluate( classifier2, XTest2 );
testAcc2 = (length(YTest2) - nnz(YTest2 - predTest2))./(length(YTest2));

fprintf( '(2) training accuracy: %f\n', trainingAcc2 );
fprintf( '(2) test     accuracy: %f\n', testAcc2 );
