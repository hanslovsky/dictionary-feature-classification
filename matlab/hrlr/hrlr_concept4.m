% hrlr_concept4

%% helper functions

downsampler = Tid.getDownsampler3dz();
downsamplerRs = Tid.getDownsampler3dzRs();

%% 3d proof-of-concept ( segmentation )

f = 3; % downsampling factor
Ntrain = 2000;
Ntest  = 2000;

classifierOpts = {  'kernel_function', 'rbf' };

patchSize = [9 9 9];
segProbs = [ 0.5 0.5 ];
numDict   = 8;
distTol   = 0.3;
patchSizeDs = patchSize ./ [1 1 f];
[Dtrue, ks] = DictionarySymSampler.dctDictionary3d( patchSize, numDict );

params = Tid.defaultParams();
params.K = 24;
params.lambda = 0.2;


%%
% deal with a binary segmentation for now

rng( 6 );

seg = mnrnd( 1, segProbs, numDict ); 
seg = seg(:,1);

rotRng = [ 1 2 3 ];
dss = DictionarySymSampler( Dtrue, [], patchSize, rotRng, [], [], [] );
[X,di] = dss.sample( Ntrain );
[X_test,di_test] = dss.sample( Ntest );

Xds = downsamplerRs( X, patchSize, f ); 
Xds_test = downsamplerRs( X_test, patchSize, f ); 

labels_train = seg(di);
labels_test  = seg(di_test);



%% the high res dictionary

tid = Tid( X, patchSize, params, distTol );
tid.buildDictionary();
size( tid.D )

tid.makeDictRotInv();
size( tid.D )

%% vis

% D = tid.D;
% N = size( D, 2);
% for i = 1:N
%     figure;
%    imdisp( permute( reshape( D(:,i), patchSize ), [ 1 2 4 3] ), ...
%        'border', 0.1 );
% %    pause;
% end

%% do a low res dictionary for comparison

tidLR = Tid( Xds, patchSizeDs, params );
tidLR.buildDictionary();
tidLR.makeDictRotInv();


%%
% hypothesis - building a dictionary as described above should
% be better for segmentationt that a dictionary build directly from 
% the low-res data

%%  classifiers on original HR data for baseline 

trn_alpha = tid.getFeaturesOrig( X );
tst_alpha = tid.getFeaturesOrig( X_test );

[ svm, trn_pred, trn_err, ...
            tst_pred, tst_err ] = classifyEvaluateSvm( trn_alpha', tst_alpha', ...
                                                labels_train, labels_test, classifierOpts );

fprintf('high-res err rate (train) : %f\n',   trn_err );
fprintf('high-res error rate (test)  : %f\n', tst_err );

%%  classifiers on original HR data xfm dict

trn_alpha_xfm = tid.getFeatures( X );
tst_alpha_xfm = tid.getFeatures( X_test );

[ svm, trn_pred_xfm, trn_err_xfm, ...
       tst_pred_xfm, tst_err_xfm ] = classifyEvaluateSvm( trn_alpha_xfm', tst_alpha_xfm', ...
                                                labels_train, labels_test, classifierOpts );

fprintf('high-res xfm error rate (train) : %f\n', trn_err_xfm );
fprintf('high-res xfm error rate (test)  : %f\n', tst_err_xfm );

%%  classifiers on original HR data xfm pooled dict

trn_alpha_xfmM = tid.mergeFeatures( tid.getFeatures( X ));
tst_alpha_xfmM = tid.mergeFeatures( tid.getFeatures( X_test ));

[ svm, trn_pred_xfmM, trn_err_xfmM, ...
       tst_pred_xfmM, tst_err_xfmM ] = classifyEvaluateSvm( trn_alpha_xfmM', tst_alpha_xfmM', ...
                                                labels_train, labels_test, classifierOpts );

fprintf('high-res xfm pooled error rate (train) : %f\n', trn_err_xfmM );
fprintf('high-res xfm pooled error rate (test)  : %f\n', tst_err_xfmM );

%%
fprintf('\nEND HR BASELINE\nSTART LR MUST BEAT\n');

%% classification on LR data ( we have to beat this )

trn_alpha_LR = tidLR.getFeaturesOrig( Xds );
tst_alpha_LR = tidLR.getFeaturesOrig( Xds_test );

[ svm_LR, trn_pred_LR, trn_err_LR, ...
            tst_pred_LR, tst_err_LR ] = classifyEvaluateSvm( trn_alpha_LR', tst_alpha_LR', ...
                                                labels_train, labels_test, classifierOpts );

fprintf('low-res err rate (train) : %f\n',   trn_err_LR );
fprintf('low-res error rate (test)  : %f\n', tst_err_LR );

%% classification on LR data xfm ( we have to beat this )

trn_alpha_LR_xfm = tidLR.getFeatures( Xds );
tst_alpha_LR_xfm = tidLR.getFeatures( Xds_test );

[ svm_LR_xfm, trn_pred_LR_xfm, trn_err_LR_xfm, ...
              tst_pred_LR_xfm, tst_err_LR_xfm ] = classifyEvaluateSvm( trn_alpha_LR_xfm', tst_alpha_LR_xfm', ...
                                                labels_train, labels_test, classifierOpts );

fprintf('low-res xfm error rate (train): %f\n', trn_err_LR_xfm );
fprintf('low-res xfm error rate (test) : %f\n', tst_err_LR_xfm );

%% classification on LR data xfm ( we have to beat this )

trn_alpha_LR_xfm = tidLR.getFeatures( Xds );
tst_alpha_LR_xfm = tidLR.getFeatures( Xds_test );

[ svm_LR_xfm, trn_pred_LR_xfm, trn_err_LR_xfm, ...
              tst_pred_LR_xfm, tst_err_LR_xfm ] = classifyEvaluateSvm( trn_alpha_LR_xfm', tst_alpha_LR_xfm', ...
                                                labels_train, labels_test, classifierOpts );

fprintf('low-res xfm error rate (train): %f\n', trn_err_LR_xfm );
fprintf('low-res xfm error rate (test) : %f\n', tst_err_LR_xfm );

%%
fprintf('\nEND HR BASELINE\nSTART LR MUST BEAT\n');

%% classification on LR data HR dictionary

trn_alpha_LRHR = tid.getFeaturesOrig( Xds, f );
tst_alpha_LRHR = tid.getFeaturesOrig( Xds_test, f );

[ svm_LRHR, trn_pred_LRHR, trn_err_LRHR, ...
            tst_pred_LRHR, tst_err_LRHR ] = classifyEvaluateSvm( trn_alpha_LRHR', tst_alpha_LRHR', ...
                                                labels_train, labels_test, classifierOpts );

fprintf('low-res HR dict error rate (train): %f\n', trn_err_LRHR );
fprintf('low-res HR dict error rate (test) : %f\n', tst_err_LRHR );

%% classification on LR data HR dictionary

trn_alpha_LRHR_xfm = tid.getFeatures( Xds, f );
tst_alpha_LRHR_xfm = tid.getFeatures( Xds_test, f );

[ svm_LRHR_xfm, trn_pred_LRHR_xfm, trn_err_LRHR_xfm, ...
                tst_pred_LRHR_xfm, tst_err_LRHR_xfm ] = classifyEvaluateSvm( trn_alpha_LRHR_xfm', tst_alpha_LRHR_xfm', ...
                                                labels_train, labels_test, classifierOpts );

fprintf('low-res HR dict xfm error rate (train): %f\n', trn_err_LRHR_xfm );
fprintf('low-res HR dict xfm error rate (test) : %f\n', tst_err_LRHR_xfm );

%% classification on LR data HR dictionary

trn_alpha_LRHR_xfmM = tid.mergeFeatures( tid.getFeatures( Xds, f ));
tst_alpha_LRHR_xfmM = tid.mergeFeatures( tid.getFeatures( Xds_test, f ));

[ svm_LRHR_xfmM, trn_pred_LRHR_xfmM, trn_err_LRHR_xfmM, ...
                 tst_pred_LRHR_xfmM, tst_err_LRHR_xfmM ] = classifyEvaluateSvm( trn_alpha_LRHR_xfmM', tst_alpha_LRHR_xfmM', ...
                                                labels_train, labels_test, classifierOpts );

fprintf('low-res HR dict pooled-xfm error rate (train): %f\n', trn_err_LRHR_xfmM );
fprintf('low-res HR dict pooled-xfm error rate (test) : %f\n', tst_err_LRHR_xfmM );
