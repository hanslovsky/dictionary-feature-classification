% hrlr_concept3

%% helper functions
% downsampler = @(x,f) reshape( mean(reshape(permute(x,[3 1 2]),3,[])), size(x)./[1 1 f]);
% downsamplerRs = @(X,sz,f) reshape(mean(reshape( X,[ sz size(X,2)]),3),[],size(X,2));
            
% downsampler = @(X,f)(permute(reshape(mean(reshape( permute(x, [3 2 1]), f, [] )), sz([3 2 1])./[f 1 1]), [3 2 1]));
% downsamplerRs = @(X,sz,f) (reshape(permute(reshape(mean(reshape(permute( reshape(X, [sz size(X,2)] ), [ 3 2 1 4]), f, [])), [sz(3)./f sz(1) sz(2) size(X,2)]), [ 3 2 1 4]), [], size(X,2) ));

downsampler = Tid.getDownsampler3dz();
downsamplerRs = Tid.getDownsampler3dzRs();

%% 3d proof-of-concept ( segmentation )

f = 3; % downsampling factor
Ntrain = 1000;
Ntest  = 1000;

classifierOpts = {  'kernel_function', 'rbf' };

patchSize = [9 9 9];
segProbs = [ 0.5 0.5 ];
numDict   = 12;
distTol   = 0.1;
patchSizeDs = patchSize ./ [1 1 f];
[Dtrue, ks] = DictionarySymSampler.dctDictionary3d( patchSize, numDict );

params = Tid.defaultParams();
params.K = 100;
params.lambda = 0.1;



%% visualize true dictionary 

% for i = 1:numDict
%    imdisp( permute( reshape( Dtrue(:,i), patchSize ), [ 1 2 4 3] ), ...
%        'border', 0.1 );
%    pause;
% end

%%
% deal with a binary segmentation for now
seg = mnrnd( 1, segProbs, numDict ); 
seg = seg(:,1);

rotRng = [ ];
dss = DictionarySymSampler( Dtrue, [], patchSize, rotRng, [], [], [] );
[X,di] = dss.sample( Ntrain );
[X_test,di_test] = dss.sample( Ntest );

Xds = downsamplerRs( X, patchSize, f ); 
Xds_test = downsamplerRs( X_test, patchSize, f ); 

labels_true = seg(di);
labels_test = seg(di_test);

%% the high res dictionar

tid = Tid( X, patchSize, params, distTol );
tid.buildDictionary();

%%

% tid.trainingFeatures( );
% 
% alpha_ds = tid.getFeaturesOrig( Xds, f ); 
% alpha_ds_test = tid.getFeaturesOrig( Xds_test, f );
% 
% alpha_ds_xfm = tid.getFeatures( Xds, f );
% alpha_ds_test_xfm = tid.getFeatures( Xds_test, f );


%% visualize the dictionary (for fun)

% D = tid.D;
% N = size( D, 2);
% for i = 1:N
%    imdisp( permute( reshape( D(:,i), patchSize ), [ 1 2 4 3] ), ...
%        'border', 0.1 );
%    pause;
% end

%% do a low res dictionary for comparison

tidLR = Tid( Xds, patchSizeDs, params );
tidLR.buildDictionary();

%%

% alphaLR = tidLR.getFeaturesOrig( Xds );
% alphaLR_test = tidLR.getFeaturesOrig( Xds_test );
% 
% alphaLR_xfm = tidLR.getFeatures( Xds );
% alphaLR_test_xfm = tidLR.getFeatures( Xds_test );


%% visualize the dictionary (for fun)

% Dlr = tidLR.D;
% N = size( Dlr, 2);
% for i = 1:N
%    imdisp( permute( reshape( Dlr(:,i), tidLR.patchSize ), [ 1 2 4 3] ), ...
%        'border', 0.1 );
%    pause;
% end

%%
% hypothesis - building a dictionary as described above should
% be better for segmentationt that a dictionary build directly from 
% the low-res data

%%  classifiers on original HR data for baseline 

[ svm, trn_pred, trn_err, ...
            tst_pred, tst_err ] = classifyEvaluateSvm( alpha', alpha_test', ...
                                                labels_true, labels_test, classifierOpts );

fprintf('high-res xfm pooled err rate (train) : %f\n',   trn_err );
fprintf('high-res xfm pooled error rate (test)  : %f\n', tst_err );

%% classifiers on original HR data for baseline 

% alpha = tid.getFeaturesOrig( X );
% alpha_test = tid.getFeaturesOrig( X_test );
% 
% svm   = svmtrain(  alpha', labels_true, 'kernel_function', 'rbf' );
% 
% lab_pred_trn = svmclassify( svm, alpha' );
% trn_err_hr = ( lab_pred_trn - labels_true );
% trn_errRate_hr = nnz( trn_err_hr ) ./ length( trn_err_hr );
% 
% lab_pred_test = svmclassify( svm, alpha_test' );
% err_hr = ( lab_pred_test - labels_test );
% errRate_hr = nnz( err_hr ) ./ length( err_hr );
% 
% fprintf('high-res error rate (train) : %f\n', trn_errRate_hr );
% fprintf('high-res error rate (test)  : %f\n', errRate_hr );

%% classifiers on original HR data,  for baseline 

% xfm_lasso_prm = params;
% xfm_lasso_prm.lambda = 0.25;
% 
% alpha_xfm = tid.getFeatures( X, [], xfm_lasso_prm );
% alpha_test_xfm = tid.getFeatures( X_test, [], xfm_lasso_prm );
% 
% svm_xfm = svmtrain(  alpha_xfm', labels_true, 'kernel_function', 'rbf' );
% 
% lab_pred_trn_xfm = svmclassify( svm_xfm, alpha_xfm' );
% trn_err_hr_xfm = ( lab_pred_trn_xfm - labels_true );
% trn_errRate_hr_xfm = nnz( trn_err_hr_xfm ) ./ length( trn_err_hr_xfm );
% 
% lab_pred_test_xfm = svmclassify( svm_xfm, alpha_test_xfm' );
% err_hr_xfm = ( lab_pred_test_xfm - labels_test );
% errRate_hr_xfm = nnz( err_hr_xfm ) ./ length( err_hr_xfm );
% 
% fprintf('high-res xfm error rate (train) : %f\n', trn_errRate_hr_xfm );
% fprintf('high-res xfm error rate (test)  : %f\n', errRate_hr_xfm );
% 
% %% classifiers on original HR data xfm Pool
% 
% alpha_xfmM      = tid.mergeFeatures( alpha_xfm );
% alpha_test_xfmM = tid.mergeFeatures( alpha_test_xfm );
% 
% svm_xfmM = svmtrain(  alpha_xfmM', labels_true, 'kernel_function', 'rbf' );
% 
% lab_pred_trn_xfmM = svmclassify( svm_xfmM, alpha_xfmM' );
% trn_err_hr_xfmM = ( lab_pred_trn_xfmM - labels_true );
% trn_errRate_hr_xfmM = nnz( trn_err_hr_xfmM ) ./ length( trn_err_hr_xfmM );
% 
% lab_pred_test_xfmM = svmclassify( svm_xfmM, alpha_test_xfmM' );
% err_hr_xfmM = ( lab_pred_test_xfmM - labels_test );
% errRate_hr_xfmM = nnz( err_hr_xfmM ) ./ length( err_hr_xfmM );
% 
% fprintf('high-res xfm pooled error rate (train) : %f\n', trn_errRate_hr_xfmM );
% fprintf('high-res xfm pooled error rate (test)  : %f\n', errRate_hr_xfmM );
% 
% %% classifiers on original HR data xfm Pool
% 
% [ svm_xfmM, trn_pred_xfmM, trn_err_xfmM, ...
%             tst_pred_xfmM, tst_err_xfmM ] = classifyEvaluateSvm( alpha_xfmM', alpha_test_xfmM', ...
%                                                 labels_true, labels_test, classifierOpts );
% 
% fprintf('high-res xfm pooled err rate (train) : %f\n',   trn_err_xfmM );
% fprintf('high-res xfm pooled error rate (test)  : %f\n', tst_err_xfmM );

%% classifiers on LR data HR dict,  for baseline 

% svm_ds = svmtrain(  alpha_ds', labels_true, 'kernel_function', 'rbf' );
% 
% lab_pred_test = svmclassify( svm_ds, alpha_ds_test' );
% err_hr_ds = ( lab_pred_test - labels_test );
% errRate_hr_ds = nnz( err_hr_ds ) ./ length( err_hr_ds );
% 
% fprintf('high-res ds error rate: %f\n', errRate_hr_ds );
% 
% %% classifiers on LR data HR dict xfm,  for baseline 
% 
% svm_ds_xfm = svmtrain(  alpha_ds_xfm', labels_true, 'kernel_function', 'rbf' );
% 
% lab_pred_test = svmclassify( svm_ds_xfm, alpha_ds_test_xfm' );
% err_hr_xfm = ( lab_pred_test - labels_test );
% errRate_hr_xfm = nnz( err_hr_xfm ) ./ length( err_hr_xfm );
% 
% fprintf('high-res ds xfm error rate: %f\n', errRate_hr_xfm );
% 
% %% classifier on LR data
% 
% svmLR   = svmtrain(   alphaLR', labels_true, 'kernel_function', 'rbf' );
% lab_LR_pred_test = svmclassify( svmLR, alphaLR_test' );
% err_lr = ( lab_LR_pred_test - labels_test );
% errRate_lr = nnz( err_lr ) ./ length( err_lr );
% 
% fprintf('low-res  error rate: %f\n', errRate_lr );
% 
% %% classifier on LR data xfm
% 
% svmLR_xfm   = svmtrain(   alphaLR_xfm', labels_true, 'kernel_function', 'rbf' );
% lab_LR_pred_test_xfm = svmclassify( svmLR_xfm, alphaLR_test_xfm' );
% err_lr_xfm = ( lab_LR_pred_test_xfm - labels_test );
% errRate_lr_xfm = nnz( err_lr_xfm ) ./ length( err_lr_xfm );
% 
% fprintf('low-res  xfm error rate: %f\n', errRate_lr_xfm );

