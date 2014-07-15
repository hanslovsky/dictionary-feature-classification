% quick start sample script
% medulla data set, small subvolumes
%
% run_script('qs_sample_medulla_sub_feats','test');

global SAVEPATH
global SAVEPREFIX
global DAWMRLIBPATH

script_name = mfilename;

%% distributed computing parameters
num_workers = 50;

%% define the data set
base_dir = sprintf('%s/projects/sample_medulla', DAWMRLIBPATH);

% mask_fn_train   = sprintf('%s/medulla_sub1_mask.h5',   base_dir);
% mask_fn_test    = sprintf('%s/medulla_sub2_mask.h5',   base_dir);
% labels_fn_train = sprintf('%s/medulla_sub1_labels.h5', base_dir);
% labels_fn_test  = sprintf('%s/medulla_sub2_labels.h5', base_dir);

data_fn_train   = sprintf('%s/medulla_sub1_data.h5',   base_dir);
mask_fn_train   = sprintf('%s/medulla_sub1_mask_1.h5',   base_dir);
labels_fn_train = sprintf('%s/medulla_sub1_labels_1.h5', base_dir);

data_fn_test    = sprintf('%s/medulla_sub2_data.h5',   base_dir);
mask_fn_test    = sprintf('%s/medulla_sub2_mask_1.h5',   base_dir);
labels_fn_test  = sprintf('%s/medulla_sub2_labels_1.h5', base_dir);

ds = dawmr_set( { data_fn_train,   '', data_fn_test   }, ...
                { labels_fn_train, '', labels_fn_test }, ...
                { mask_fn_train,   '', mask_fn_test   }, ...
                data_fn_train );
% ds = dawmr_set( { data_fn_test,   '', data_fn_train   }, ...
%                 { labels_fn_test, '', labels_fn_train }, ...
%                 { mask_fn_test,   '', mask_fn_train   }, ...
%                 data_fn_test );
ds.affinity_edges = []

%% specify unsupervised architecture
c_thresh_pol_kmeans = 3;
c_ave_pooling       = 0;
c_max_pooling       = 1;

patch_dim           = 5;

pooling_radius      = 2
pooling_type        = c_max_pooling;

num_clusters        = 1000
num_patches_kmeans  = 10000;
num_train_instances = Inf;
num_test_instances  = Inf;

dc = dawmr_clustering(patch_dim, num_clusters);
dc = dc_foveated(dc, pooling_radius, pooling_type, ...
                 c_thresh_pol_kmeans);

%% specify mlp parameters
mlp_init = mlp(100);
mlp_init.num_updates_default = 5e4;
mlp_init.minibatch_size_default = 40;
mlp_init.use_gpu = 1;
mlp_init.pos_ratio = 0.5;
mlp_init.margin = 0.4;
mlp_init.eta_w_start = [0.02 0.02];
mlp_init.eta_b_start = [0.02 0.02];
mlp_init.eta_b_start = [0.02 0.02];
mlp_init.loss_func = 1;

%% segmentation parameters
thresh           = [0.5 0.7 0.8 0.9 0.95 0.99];
watershed_thresh = [0];


%% set up model, do learning
ts = tic;

dm = dawmr(ds, 3, ec_mlp_h5(mlp_init), script_name);
dm = dawmr_add_multilayer(dm, dc, [1]);
dm.learn_features(num_patches_kmeans, 2, [],[],[],50);

if(~isempty(SAVEPATH))
  save(sprintf('%s/%s_%s_learnedFeatures.mat', ...
               SAVEPATH, SAVEPREFIX, script_name), ...
       '-v7.3');
end

[accs_train, labels_gt_train, ~, labels_pd_train] = ...
    dm.comp_features(1, num_train_instances, ...
                  [],[],[],[], num_workers);
              
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
                  [],[],[],[], num_workers, 0);
te = toc(ts);

% temporary save
if(~isempty(SAVEPATH))
  save(sprintf('%s/%s_%s_results.mat', ...
               SAVEPATH, SAVEPREFIX, script_name), ...
       '-v7.3');
end


