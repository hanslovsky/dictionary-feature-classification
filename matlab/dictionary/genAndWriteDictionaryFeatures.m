% genAndWriteDictionaryFeatures
%
% dbstop if error; genAndWriteDictionaryFeatures;

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
param.verbose = true;

%%

% feats = getTrainTestData( ds.data_fn{1}, ds.labels_fn{1}, ds.mask_fn{1}, ...
%                            N, D, patch_size, param );
% 
% success = writeFeaturesH5( feats, [], '/nobackup/saalfeld/john/forPhilipp/ds_18_dict1k_faster.h5' );
% success

%%

sample_dir = '/groups/saalfeld/home/bogovicj/dev/dawmr/dawmr_lib_public/projects/sample_medulla/';
data_fn_samp_trn = fullfile(sample_dir, 'medulla_sub1_data.h5');
labels_fn_samp_trn = fullfile(sample_dir, 'medulla_sub1_labels_1.h5');
mask_fn_samp_trn = fullfile(sample_dir, 'medulla_sub1_mask_1.h5');

feats_trn = getTrainTestData( data_fn_samp_trn, labels_fn_samp_trn, mask_fn_samp_trn, ...
                            [], D, patch_size, param );
                        
success_trn = writeFeaturesH5( feats_trn, [], '/nobackup/saalfeld/john/forPhilipp/med_sample_train_wLabels.h5' );
success_trn 

%
sample_dir = '/groups/saalfeld/home/bogovicj/dev/dawmr/dawmr_lib_public/projects/sample_medulla/';
data_fn_samp_tst = fullfile(sample_dir, 'medulla_sub2_data.h5')
labels_fn_samp_tst = fullfile(sample_dir, 'medulla_sub2_labels_1.h5')
mask_fn_samp_tst = fullfile(sample_dir, 'medulla_sub2_mask_1.h5')

feats_tst = getTrainTestData( data_fn_samp_tst, labels_fn_samp_tst, mask_fn_samp_tst, ...
                            [], D, patch_size, param );
                        
success_tst = writeFeaturesH5( feats_tst, [], '/nobackup/saalfeld/john/forPhilipp/med_sample_test_wLabels.h5' );
success_tst