%% segWith2dTo3dDict

%% helper functions 

classifySVM = @( X, Y, opts )( svmtrain(  X, Y, opts{:} ) );
evaluateSVM = @( svm, X )( svmclassify( svm, X ) );

classifyRF = @( X, Y, numTrees, opts )( TreeBagger(numTrees, X, Y, opts{:}) );
evaluateRF = @( rf, X )( cellfun( @str2double, rf.predict(X)) );

%% params and data

ds_training = 18;
ds_test     =  3;
ds = ds_medulla(ds_training, ds_test);

% dictionary building parameters
param.K = 100;  % dictionary size
param.lambda=0.1;
param.numThreads=4; % number of threads
param.verbose=0;
param.iter = 250;  % num iterations.

num = 50;
patchSize = [9 9];
dsFactor  = 3;

N = 10000;
numTrees = 100;

basedir = '/groups/saalfeld/home/bogovicj/projects/dictionary/dict2dTo3d';

%% load training images

im_trn = ds.data_fn{1};
lb_trn = ds.labels_fn{1};
mk_trn = ds.mask_fn{1};

im_test = ds.data_fn{3};
lb_test = ds.labels_fn{3};
mk_test = ds.mask_fn{3};


%% grab data for building dictionary

% [X_test, coords] = grabPatchesSimple( im, patchSize, N, [], mk );

[~, X_trn_cell, xyz_trn ] = getTrainTestData( im_trn, lb_trn, mk_trn, ...
                                              [N N], [], [patchSize dsFactor], param, dsFactor );

X_trn = [ X_trn_cell{1}; X_trn_cell{2} ];
Y_trn = [ zeros(N,1); ones(N,1) ];

%%

D = mexTrainDL( X_trn', param );

destdir = sprintf('%s%sdict_%s', basedir, filesep, datestr(now, 30));
mkdir( destdir );

save( fullfile( destdir, 'dict'), 'D', 'param');

%% train classifier

alpha_trn = encoding( X_trn', D, 'sc', [], param );
alpha_trn = full(alpha_trn');

fprintf('training classifier\n');
classifier = classifyRF( alpha_trn, Y_trn, numTrees, {});

pred_trn  = evaluateRF( classifier, alpha_trn );
trn_errors = ( pred_trn - Y_trn );

trn_fpc    = nnz( trn_errors ==  1 );
trn_fnc    = nnz( trn_errors == -1 );
trn_tpc    = nnz( trn_errors == 0 & Y_trn == 1 );
trn_tnc    = nnz( trn_errors == 0 & Y_trn == 0 );

trn_acc = nnz( trn_errors == 0 )./length(Y_trn);
trn_acc

%% evaluate classifier 

[alpha_tst_cell, X_tst_cell, xyz_tst] = getTrainTestData( ds.data_fn{3}, ds.labels_fn{3}, ds.mask_fn{3}, ...
                                                          [N N], D, [patchSize dsFactor], param, dsFactor );
Y_tst = [ zeros(N,1); ones(N,1) ];

alpha_tst = full([ alpha_tst_cell{1}'; alpha_tst_cell{2}' ]);
pred_tst  = evaluateRF( classifier, alpha_tst );

% tst_errors is: 0 for correct examples
%                1 for false positives
%               -1 for false negatives
tst_errors = ( pred_tst - Y_tst );

tst_fpc    = nnz( tst_errors ==  1 );
tst_fnc    = nnz( tst_errors == -1 );
tst_tpc    = nnz( tst_errors == 0 & Y_tst == 1 );
tst_tnc    = nnz( tst_errors == 0 & Y_tst == 0 );

tst_acc = nnz( tst_errors == 0 )./length(Y_tst);
tst_acc

%% try building a 3d dictionary from the 2d one

d23 = Dict2dTo3d( D', patchSize(1), dsFactor );
