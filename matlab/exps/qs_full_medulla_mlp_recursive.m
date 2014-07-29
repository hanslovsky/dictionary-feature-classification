% quick start sample script
% medulla data set, small subvolumes
%
% cd ~/dev/main/dictionary-feature-classification/matlab/exps/
%
% dbstop if error; run_script('qs_full_medulla_mlp_recursive','HR dawmr recursive iter 2, 2k clusters, 5e5 iters');
% dbstop if error; run_script('qs_full_medulla_mlp_recursive','LR dawmr recursive iter 2, 2k clusters, 5e5 iters');
% 
% dbstop if error; run_script('qs_full_medulla_mlp_recursive','dictionary only, 2k clusters');


global SAVEPATH
global SAVEPREFIX

%% exp params
 
% HR MLP result
saved_mlp = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0102_qs_full_medulla_mlp/exp0102_qs_full_medulla_mlp_finalResults.mat';

% LR MLP result
% saved_mlp = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0101_qs_full_medulla_mlp/exp0101_qs_full_medulla_mlp_finalResults.mat'

% HIGH RES DICTIONARY 1k clusters
% saved_dict_fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0080_qs_full_medulla_mlp/exp0080_qs_full_medulla_mlp_learnedFeatures.mat';

% LOW RES DICTIONARY 1k clusters
% saved_dict_fn = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0084_qs_full_medulla_mlp/exp0084_qs_full_medulla_mlp_learnedFeatures.mat';

saved_dict_fn = '';

do_classifier   = 1;
do_inference    = 1;

use_downsampled = 1;
downsample_dict = 0;

factor = [1 1 3];

do_foveated = 1;
do_cn       = 0;

num_workers = 50;

num_clusters        = 2000;
num_patches_kmeans  = 10000;
feat_iters 			= 300;


%% define the data set

training = 18;
test = 3;
ds = ds_medulla(training, test);
suffix = '_toBnd';

data_fn_train = ds.data_fn{1};
data_fn_test  = ds.data_fn{3};

mask_fn_train = ds.mask_fn{1};
mask_fn_test  = ds.mask_fn{3};

labels_fn_train = ds.labels_fn{1};
labels_fn_test  = ds.labels_fn{3};

datdir = '/groups/saalfeld/home/bogovicj/projects/aniso/downsamp/medulla/dsdat';
[~,trnVolName] = fileparts( data_fn_train );
[~,tstVolName] = fileparts( data_fn_test );
trnVolDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d%s.h5',trnVolName,factor,suffix));
tstVolDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d%s.h5',tstVolName,factor,suffix));

[~,trnMskName] = fileparts( mask_fn_train );
[~,tstMskName] = fileparts( mask_fn_test );
trnMskDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d%s.h5',trnMskName,factor,suffix));
tstMskDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d%s.h5',tstMskName,factor,suffix));

[~,trnLabName] = fileparts( labels_fn_train );
[~,tstLabName] = fileparts( labels_fn_test );
trnLabDsFn = fullfile(datdir, sprintf('trn_%s_ds%d-%d-%d%s.h5',trnLabName,factor,suffix));
tstLabDsFn = fullfile(datdir, sprintf('tst_%s_ds%d-%d-%d%s.h5',tstLabName,factor,suffix));

if ( use_downsampled )
    ds = dawmr_set( { trnVolDsFn,   '', tstVolDsFn   }, ...
                    { trnLabDsFn, '', tstLabDsFn }, ...
                    { trnMskDsFn,   '', tstMskDsFn   }, ...
                    trnVolDsFn );
end

ds.data_fn{1}
ds.data_fn{3}


%% specify unsupervised architecture

alphas    = 0;

c_thresh_pol_kmeans = 3;
c_ave_pooling       = 0;
c_max_pooling       = 1;
patch_dim           = [9 9 3];

encoding_type      = c_thresh_pol_kmeans;
pooling_radius     = [2 2 1];
pooling_type       = c_max_pooling;

num_train_instances = 0.1;
num_test_instances  = Inf;

feature_normalization = 3;
% feature_normalization = 0;

if (do_foveated)
    dc = dawmr_clustering(patch_dim, num_clusters, do_cn, ...
        [], [],[],[], [],[], 0);
    dc = dc_foveated(dc, pooling_radius, pooling_type, encoding_type );
else
    dc = dawmr_clustering(patch_dim, num_clusters);
    dp_cen  = dawmr_pooling(pooling_type, 0, [0 0 0]);
    dc.add_dp(dp_cen);
end

downsampling_type = 1;


%% specify mlp parameters

mlp_init = mlp( 200 );
mlp_init.num_updates_default = 5e5;
mlp_init.minibatch_size_default = 40;
mlp_init.use_gpu = 1;
mlp_init.pos_ratio = 0.5;
mlp_init.margin = 0.4;
mlp_init.eta_w_start = [0.02 0.02];
mlp_init.eta_b_start = [0.02 0.02];
mlp_init.loss_func = 1;

% dm.end_classifier = ec_mlp_fs(mlp_init);

%%

if( ~isempty(saved_mlp) )
    disp('loading mlp') 
	oldMlp = load( saved_mlp, 'dm', 'SAVEPATH', 'SAVEPREFIX','script_name'); 
	SAVEPATH_old = oldMlp.SAVEPATH;
	SAVEPREFIX_old = oldMlp.SAVEPREFIX;
	script_name_old = oldMlp.script_name;
	dm = oldMlp.dm;


	ag_fn_training = sprintf('%s/%s_%s_train_infer.h5', ...
	   SAVEPATH_old, SAVEPREFIX_old, ...
	   script_name_old);
	ag_fn_testing  = sprintf('%s/%s_%s_test_infer.h5', ...
	   SAVEPATH_old, SAVEPREFIX_old, ...
	   script_name_old);
end

%disp( 'paused');
%pause; 

if( isempty( saved_dict_fn))

    dm2 = dawmr(ds, 3, ec_mlp_fs(mlp_init), script_name);
    dm2 = dawmr_add_multilayer(dm2, dc, [1]);
	dm2.ds.data_extra_channels_fn = ag_fn_training;
    dm2.learn_features(num_patches_kmeans, 2, [],[],[], feat_iters);

    if(~isempty(SAVEPATH))
        save(sprintf('%s/%s_%s_learnedFeatures.mat', ...
            SAVEPATH, SAVEPREFIX, script_name), ...
            'dm2','-v7.3');
    end
else
    % load the saved dm w/ dictionary params
    load( saved_dict_fn );
end

if( downsample_dict )
    dm2 = resampleDawmrClusters( dm2, factor );
    dm2.dds.dcs
end

if( do_classifier )
    ts = tic;
    % train 
    [accs_train, labels_gt_train, ~, labels_pd_train] = ...
        dm2.classifier(1, num_train_instances, ...
            [],[],[],[], num_workers);
      
        
	dm2.ds.data_extra_channels_fn = ag_fn_testing;
    % test 
    [accs, labels_gt, ~, labels_pd, aucs] = ...
        dm2.classifier(3, num_test_instances, ...
            [],[],[],[], num_workers);
    te = toc(ts);
    
    fprintf('second model\n');
    dawmr_set_print_stats(accs_train, accs, aucs, te, [], ...
        labels_gt_train, labels_pd_train, ...
        labels_gt, labels_pd, ...
        SAVEPREFIX, script_name, dm2, num_clusters);
    
    % results save
    if(~isempty(SAVEPATH))
        save(sprintf('%s/%s_%s_finalResults.mat', ...
            SAVEPATH, SAVEPREFIX, script_name), ...
            '-v7.3');
    end
    
end

if( do_inference )
    train_infer_fn = sprintf('%s/%s_%s_train_infer.h5', ...
            SAVEPATH, SAVEPREFIX, script_name)

    test_infer_fn = sprintf('%s/%s_%s_test_infer.h5', ...
            SAVEPATH, SAVEPREFIX, script_name)

    % training
    dm.infer([],[],train_infer_fn,...
            [],[],[],[20 20 20],[], 1);
        
    % testing
    dm.infer([],[],test_infer_fn,...
            [],[],[],[20 20 20],[], 3);

end
