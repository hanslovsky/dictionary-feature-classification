%% recon2dTo3dDict
%
% dbstop if error; run_script('recon2dTo3dDict', 'k=100, N=1000 );
% dbstop if error; run_script('recon2dTo3dDict', 'k=200, N=100000 );

global SAVEPATH
global SAVEPREFIX

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

N = 1000;
numTrees = 100;
dictSize3d = 200;

doClassification = 0;

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

% generate downsampled patches
[~, X_trn_cell, xyz_trn ] = getTrainTestData( im_trn, lb_trn, mk_trn, ...
                                              [N N], [], [patchSize dsFactor], param, dsFactor );

X_trn = [ X_trn_cell{1}; X_trn_cell{2} ];

D = mexTrainDL( X_trn', param );

% save( fullfile( SAVEPATH, [SAVEPREFIX,'_dict']), 'D', 'param');

%% try building a 3d dictionary from the 2d one


clear d23; clc;
d23 = Dict2dTo3d( D', patchSize(1), dsFactor );

% d23.build3dDictionary( dictSize3d );
% save( fullfile( SAVEPATH, [SAVEPREFIX,'_d23']), 'd23', 'param');


[rc, lc] = rootConstraintFromPatch( [10 9 8], 3, 3, 3 );
iniParams = Dict2dTo3d.patchParamNodeToArray( lc )

patchParams = d23.build3dPatch( lc );

params = Dict2dTo3d.patchParamNodeToArray( patchParams )

isConsistent = 1;
colset = [ 1:3; 1:3:7; ]
for col = 1:2
    for sv = 1:3
        isConsistent = isConsistent & ( nnz(params(:,col) == colset(col, sv)) == 3 );
    end
end
isConsistent

%% test supplementing old dictionary with new patches

D = rand( 81, 10 );
patchSize = [ 9 9 9 ];
dsFactor = 3;

clear d23; clc;
d23 = Dict2dTo3d( D', patchSize(1), dsFactor );

patch = rand( 3, 81 );

numAdded = d23.addTemporaryPatchesToDict( patch );
iniIdx = d23.numDict-numAdded+1:d23.numDict;

paramList = [ 3 3 3; 1 4 7; iniIdx]'

[rc, lc] = rootConstraintFromArray( paramList );

rc.printAsTree()

%%

newSims = d23.addToAllSims( numAdded );

patchParams = d23.build3dPatch( lc, numAdded );

[pv, patch] = d23.patchFromParams( patchParams );
imdisp3d( patch )

% d23.removeTemporaryPatches( numAdded );

