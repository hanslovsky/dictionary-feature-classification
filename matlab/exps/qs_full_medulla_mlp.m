% quick start sample script
% medulla data set, small subvolumes
%
% cd ~/dev/main/dictionary-feature-classification/matlab/exps/
%
% dbstop if error; run_script('qs_full_medulla_mlp','dictionary only, 2k clusters');
% dbstop if error; run_script('qs_full_medulla_mlp','dictionary only, downsampZ 3, patch[9 9 3], 1k clusters');
% dbstop if error; run_script('qs_full_medulla_mlp','dictionary only, downsampZ 3, patch[9 9 3], 2k clusters');
%
% dbstop if error; run_script('qs_full_medulla_mlp','mlp100, 1k clust');
% dbstop if error; run_script('qs_full_medulla_mlp','ds3, mlp100, 1k clust 9-9-3');
%
% dbstop if error; run_script('qs_full_medulla_mlp','ds3, mlp100, 1k clust 9-9-3 dsTo 9-9-3');


global SAVEPATH
global SAVEPREFIX

%% exp params

saved_dict_fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0080_qs_full_medulla_mlp/exp0080_qs_full_medulla_mlp_learnedFeatures.mat';
% saved_dict_fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0084_qs_full_medulla_mlp/exp0084_qs_full_medulla_mlp_learnedFeatures.mat';
% saved_dict_fn = '';

do_classifier   = 1;
use_downsampled = 1;
downsample_dict = 1;
factor = [1 1 3];

num_workers = 50;

num_training_patches = 25000;
num_testing_patches  = 5000;

%% define the data set

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
trnVolDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d.h5',trnVolName,factor));
tstVolDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d.h5',tstVolName,factor));

[~,trnMskName] = fileparts( mask_fn_train );
[~,tstMskName] = fileparts( mask_fn_test );
trnMskDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d.h5',trnMskName,factor));
tstMskDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d.h5',tstMskName,factor));

[~,trnLabName] = fileparts( labels_fn_train );
[~,tstLabName] = fileparts( labels_fn_test );
trnLabDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d.h5',trnLabName,factor));
tstLabDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d.h5',tstLabName,factor));

if ( use_downsampled )
    ds = dawmr_set( { trnVolDsFn,   '', tstVolDsFn   }, ...
                    { trnLabDsFn, '', tstLabDsFn }, ...
                    { trnMskDsFn,   '', tstMskDsFn   }, ...
                    trnVolDsFn );
end

ds.data_fn{1}
ds.data_fn{3}


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

% dm.end_classifier = ec_mlp_fs(mlp_init);

%%

if( isempty(saved_dict_fn) )
    
    dm = dawmr(ds, 3, ec_mlp_fs(mlp_init), script_name);
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
end

if( downsample_dict )
    dm = resampleDawmrClusters( dm, factor );
    dm.dds.dcs
end

if( do_classifier )
    ts = tic;
    % train 
    [accs_train, labels_gt_train, ~, labels_pd_train] = ...
        dm.classifier(1, num_train_instances, ...
            [],[],[],[], num_workers);
      
    % test 
    [accs, labels_gt, ~, labels_pd, aucs] = ...
        dm.classifier(3, num_test_instances, ...
            [],[],[],[], num_workers);
    te = toc(ts);
    
    fprintf('first model\n');
    dawmr_set_print_stats(accs_train, accs, aucs, te, [], ...
        labels_gt_train, labels_pd_train, ...
        labels_gt, labels_pd, ...
        SAVEPREFIX, script_name, dm, num_clusters);
    
    % results save
    if(~isempty(SAVEPATH))
        save(sprintf('%s/%s_%s_finalResults.mat', ...
            SAVEPATH, SAVEPREFIX, script_name), ...
            '-v7.3');
    end
    
end
