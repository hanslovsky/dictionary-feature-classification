% quick start sample script
% medulla data set, small subvolumes
%
% cd ~/dev/main/dictionary-feature-classification/matlab/exps/
%
% dbstop if error; run_script('qs_full_medulla_feats',' LR image, 1k clusters from HR');
% dbstop if error; run_script('qs_full_medulla_feats',' LR image, 1k clusters from LR');

global SAVEPATH
global SAVEPREFIX

%% exp params

%saved_dict_fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0080_qs_full_medulla_mlp/exp0080_qs_full_medulla_mlp_learnedFeatures.mat';
saved_dict_fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0084_qs_full_medulla_mlp/exp0084_qs_full_medulla_mlp_learnedFeatures.mat';
% saved_dict_fn = '';

use_downsampled = 1;
downsample_dict = 0;
factor = [1 1 3];

num_workers = 50;

num_training_patches = 25000;
num_testing_patches  = 5000;

feature_dir_base = '/nobackup/saalfeld/john/medulla_features';

%% define the data set

AgToBp = true;
if( AgToBp )
    bpSuffix = '_toBnd';
else
    bpSuffix = '';
end

training = 18;
test = 3;
ds = ds_medulla(training, test);

data_fn_train = ds.data_fn{1};
data_fn_test  = ds.data_fn{3};

mask_fn_train = ds.mask_fn{1};
mask_fn_test  = ds.mask_fn{3};

labels_fn_train = ds.labels_fn{1};
labels_fn_test  = ds.labels_fn{3};

datdir = '/groups/saalfeld/home/bogovicj/projects/aniso/downsamp/medulla/dsdat';
[~,trnVolName] = fileparts( data_fn_train );
[~,tstVolName] = fileparts( data_fn_test );
trnVolDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d%s.h5',trnVolName,factor,bpSuffix));
tstVolDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d%s.h5',tstVolName,factor,bpSuffix));

[~,trnMskName] = fileparts( mask_fn_train );
[~,tstMskName] = fileparts( mask_fn_test );
trnMskDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d%s.h5',trnMskName,factor,bpSuffix));
tstMskDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d%s.h5',tstMskName,factor,bpSuffix));

[~,trnLabName] = fileparts( labels_fn_train );
[~,tstLabName] = fileparts( labels_fn_test );
trnLabDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d%s.h5',trnLabName,factor,bpSuffix));
tstLabDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d%s.h5',tstLabName,factor,bpSuffix));

if ( use_downsampled )
    ds = dawmr_set( { trnVolDsFn, '', tstVolDsFn   }, ...
                    { trnLabDsFn, '', tstLabDsFn }, ...
                    { trnMskDsFn, '', tstMskDsFn   }, ...
                    trnVolDsFn );
end

ds.data_fn{1}
ds.data_fn{3}
ds.affinity_edges = [];

%% specify unsupervised architecture
c_thresh_pol_kmeans = 3;
c_ave_pooling       = 0;
c_max_pooling       = 1;
patch_dim           = [9 9 3]

pooling_radius      = 2;
pooling_type        = c_max_pooling;

num_clusters        = 1000;
num_patches_kmeans  = 10000;
num_train_instances = Inf;
num_test_instances  = Inf;

feature_normalization = 3;
% feature_normalization = 0;

dc = dawmr_clustering(patch_dim, num_clusters);
dp_cen  = dawmr_pooling(pooling_type, 0, [0 0 0]);
dc.add_dp(dp_cen);

downsampling_type = 1;

%% specify mlp parameters

mlp_init = mlp( 100 );
mlp_init.num_updates_default = 5e5;
mlp_init.minibatch_size_default = 40;
mlp_init.use_gpu = 1;
mlp_init.pos_ratio = 0.5;
mlp_init.margin = 0.4;
mlp_init.eta_w_start = [0.02 0.02];
mlp_init.eta_b_start = [0.02 0.02];
mlp_init.eta_b_start = [0.02 0.02];
mlp_init.loss_func = 1;

%%

if( isempty(saved_dict_fn) )
    
    dm = dawmr(ds, 3, ec_mlp_h5(mlp_init), script_name);
    dm = dawmr_add_multilayer(dm, dc, [1]);
    dm.learn_features(num_patches_kmeans, 2, [],[],[], 50);
    
    if(~isempty(SAVEPATH))
        save(sprintf('%s/%s_%s_learnedFeatures.mat', ...
            SAVEPATH, SAVEPREFIX, script_name), ...
            'dm','-v7.3');
    end
else
    % load the saved dm w/ dictionary params
    load( saved_dict_fn );
    dm.ds = ds
end
dm.end_classifier =  ec_mlp_h5(dm.end_classifier.mlp_init);


if( downsample_dict )
    dm = resampleDawmrClusters( dm, factor );
    dm.dds.dcs
end

%%

[accs_train, labels_gt_train, ~, labels_pd_train] = ...
    dm.comp_features(1, num_train_instances, ...
                  [],[],[],[], num_workers, [], ...
                  fullfile(feature_dir_base,'train'));
              
if(~isempty(SAVEPATH))
  save(sprintf('%s/%s_%s_dawmrObj.mat', ...
               SAVEPATH, SAVEPREFIX, script_name), ...
       'dm','-v7.3');
end

% swap image data
tmp = dm.ds.data_fn{1};
dm.ds.data_fn{1} = dm.ds.data_fn{3};
dm.ds.data_fn{3} = tmp;

% swap label data
tmp = dm.ds.labels_fn{1};
dm.ds.labels_fn{1} = dm.ds.labels_fn{3};
dm.ds.labels_fn{3} = tmp;

% swap mask data
tmp = dm.ds.mask_fn{1};
dm.ds.mask_fn{1} = dm.ds.mask_fn{3};
dm.ds.mask_fn{3} = tmp;

[accs_tst, labels_gt_test, ~, labels_pd_test ] = ...
    dm.comp_features(1, num_test_instances, ...
                  [],[],[],[], num_workers, 0, ...
                  fullfile(feature_dir_base,'test'));
te = toc(ts);

% temporary save
if(~isempty(SAVEPATH))
  save(sprintf('%s/%s_%s_results.mat', ...
               SAVEPATH, SAVEPREFIX, script_name), ...
       '-v7.3');
end

