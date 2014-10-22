%% segWith2dTo3dDict
%
% dbstop if error; run_script('segWith2dTo3dDict', 'k=200, N=100000, iters=500');
% dbstop if error; run_script('segWith2dTo3dDict', 'k=200, N=100000, iters=500, build 3d dict( 200 )');
% dbstop if error; run_script('segWith2dTo3dDict', 'k=200, N=100000, iters=500, build 3d dict( 200 ), constrained min');
% dbstop if error; run_script('segWith2dTo3dDict', 'k=1000, N=1000000, iters=500, build 3d dict dist( 1000 ), unconstrained min');
%
% dbstop if error; run_script('segWith2dTo3dDict', 'test consistency UC');
% dbstop if error; run_script('segWith2dTo3dDict', 'test consistency C');
% dbstop if error; run_script('segWith2dTo3dDict', 'test distributed');

global SAVEPATH
global SAVEPREFIX

%% helper functions 

classifySVM = @( X, Y, opts )( svmtrain(  X, Y, opts{:} ) );
evaluateSVM = @( svm, X )( svmclassify( svm, X ) );

classifyRF = @( X, Y, numTrees, opts )( TreeBagger(numTrees, X, Y, opts{:}) );
evaluateRF = @( rf, X )( cellfun( @str2double, rf.predict(X)) );

%% params and data

% set random seed
rng(42);

ds_training = 18;
ds_test     =  3;
ds = ds_medulla(ds_training, ds_test);

% dictionary building parameters
param.K = 1000;  % dictionary size
param.lambda=0.1;
param.numThreads=3; % number of threads
param.verbose=0;
param.iter = 500;  % num iterations.

num = 50;
patchSize = [9 9];
dsFactor  = 3;

N = 1000000;

numTrees = 100;

dictSize3d = 1000;
pe_batch = 1;

doClassification = 0;
constrainedMin   = 0;

Dfile = [];
Dfile = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0156_segWith2dTo3dDict/exp0156_dict.mat';

% basedir = '/groups/saalfeld/home/bogovicj/projects/dictionary/dict2dTo3d';

%% load training images

im_trn = ds.data_fn{1};
lb_trn = ds.labels_fn{1};
mk_trn = ds.mask_fn{1};

im_test = ds.data_fn{3};
lb_test = ds.labels_fn{3};
mk_test = ds.mask_fn{3};

%% grab data for building dictionary

if( isempty( Dfile ))
    % [X_test, coords] = grabPatchesSimple( im, patchSize, N, [], mk );
    
    [~, X_trn_cell, xyz_trn ] = getTrainTestData( im_trn, lb_trn, mk_trn, ...
        [N N], [], [patchSize dsFactor], param, dsFactor );
    
    X_trn = [ X_trn_cell{1}; X_trn_cell{2} ];
    Y_trn = [ zeros(N,1); ones(N,1) ];
    
    %%
    
    
    D = mexTrainDL( X_trn', param );
    
    % destdir = sprintf('%s%sdict_%s', basedir, filesep, datestr(now, 30));
    % mkdir( destdir );
    
    save( fullfile( SAVEPATH, [SAVEPREFIX,'_dict']), 'D', 'param');
else
    fprintf('loading dictionary from file: %s\n', Dfile);
    load( Dfile );
end

%% train classifier
if ( doClassification )
    
alpha_trn = encoding( X_trn', D, 'sc', [], param );
alpha_trn = full(alpha_trn');

fprintf('training classifier\n');
classifier = classifyRF( alpha_trn, Y_trn, numTrees, {});

fprintf('saving classifier\n');
save( fullfile( SAVEPATH, [SAVEPREFIX,'_classifierRF']), 'classifier' );

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

end

%% try building a 3d dictionary from the 2d one

clear d23; 
if( constrainedMin )
    fprintf('Constrained min\n');
    d23 = Dict2dTo3dConstr( D', patchSize(1), dsFactor );
else
    fprintf('Un-constrained min\n');
    d23 = Dict2dTo3d( D', patchSize(1), dsFactor );    
end

% d23.build3dDictionary( dictSize3d );
[patches, out2 ] = d23.build3dDictionaryDist( dictSize3d, pe_batch );

% clear the saved object
system( sprintf('rm -v %s', d23.obj_fn ));

save( fullfile( SAVEPATH, [SAVEPREFIX,'_d23']), 'd23', 'param');

%%

% printme = patchParams
% params = [  patchParams.getData().dim, ...
%             patchParams.getData().xyz, ...
%             patchParams.getData().idx ];
%         
% i = 1;
% while( ~printme.isRoot())
%     printme = printme.getParent()
%     params = [ params; ...
%                 printme.getData().dim, ...
%                 printme.getData().xyz, ...
%                 printme.getData().idx ];
%                 
% end
% params

