% genAndWriteDictionaryFeatures

% depends on 

%% param
ds_training = 18;
ds_test     =  3;

dict_fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0060_emSparseCoding_varyClusters/exp0060_results.mat';

N = [2000 2000];

%% data set & dictionary

ds = ds_medulla(ds_training, ds_test);

load( dict_fn, 'D_list', 'param', 'patch_size');
D = D_list{5};
clear D_list;

% make sure num clusters is correct
param.K = size(D,2);

%%

feats = getTrainTestData( ds.data_fn{1}, ds.labels_fn{1}, ds.mask_fn{1}, ...
                            N, D, patch_size, param );

success = writeFeaturesH5( feats, [], '/nobackup/saalfeld/john/forPhilipp/ds_18_dict1k.h5' );

%%

sample_dir = '/groups/saalfeld/home/bogovicj/dev/dawmr/dawmr_lib_public/projects/sample_medulla/';
data_fn_samp = fullfile(sample_dir, 'medulla_sub1_data.h5');

feats = getTrainTestData( data_fn_samp, [], [], ...
                            [], D, patch_size, param );
                        
success = writeFeaturesH5( feats, [], '/nobackup/saalfeld/john/forPhilipp/ds_tst3_dict1k.h5' );