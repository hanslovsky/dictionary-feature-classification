%% buildDict2d
% dbstop if error; run_script('buildDict2d', 'test bubble data, f=1');
%
% dbstop if error; run_script('buildDict2d', 'test microtubules, build D2d, ps[5], N=1m, K=500, f=5');
%
% dbstop if error; run_script('buildDict2d', 'rick-0512-6, TEST, build D2d, ps[9], N=1m, K=1k');

global SAVEPATH
global SAVEPREFIX

%% 

% srcdir = '/groups/saalfeld/home/bogovicj/projects/dictionary/dict2dTo3d/sim/bubbles';
% train_fn = fullfile( srcdir, 'Bubbles_r3_b1_snr10_1.tif');
% test_fn  = fullfile( srcdir, 'Bubbles_r3_b1_snr10_2.tif');
% mskFun = @(im)( im > 0 );
% dsFor2dDict = 1;

train_fn = fullfile( '/groups/saalfeld/saalfeldlab/data/virtual/microtubules/1_5', '*.png');
test_fn  = fullfile( '/groups/saalfeld/saalfeldlab/data/virtual/microtubules/1_1', '*.png');
mskFun = @(im)(im < 255 );
dsFor2dDict = 0;

% train_fn = '/groups/saalfeld/saalfeldlab/data/rick/0512-6_Stack4ECS_Segmentation.corrected.0.02-8-5-0.001-0.001-0.001.tif';
% im = readMultiTiff( train_fn );
% train_im = downSampleGaussianImglib( im, [2 2 1], [0.5 0.5 0], [0.5 0.5 0]);
% clear im;
% mskFun = @(im)( true(size(im)) );
% dsFor2dDict = 0;

%% params

N = 1e6;
% N = 1e3;

% dsFactor = 3;
% patchSize = [5 5 dsFactor];

dsFactor = 5;
patchSize = [9 9 dsFactor];

% dsFactor = 5;
% patchSize = [15 15 dsFactor];

% dictionary building parameters
param.K = 1000;  % dictionary size
param.lambda = 0.1;
param.numThreads = 1; % number of threads
param.verbose = 1;
param.iter = 500;  % num iterations.

normData = 1;

useNMF = false;

%% build 2d dictionary

if( ~exist( 'Dfile', 'var') || isempty( Dfile ))
    
    if( ~exist( 'train_im', 'var') )
        if( strcmp( train_fn(end-2:end), 'png' ))
            train_im = readImageStack( train_fn );
        else
            train_im = readMultiTiff( train_fn );
        end
    end
    
    if( dsFor2dDict && (dsFactor > 1 ))
        X_trn = grabPatchesSimple( train_im, patchSize, N, [], mskFun( train_im ) );
        fprintf('downsampling for 2d dictionary building\n');
        downsampler = Tid.getDownsampler3dzRs();
        X_trn = [downsampler( X_trn', patchSize, dsFactor )]';
    else
        X_trn = grabPatchesSimple( train_im, [patchSize(1:2) 1], N, [], mskFun( train_im ) );
    end
    
    if( normData )
        X_trn = bsxfun( @rdivide, X_trn, sqrt(sum(X_trn.^2,2)));
    end
    
    varlist = { 'D2d', 'param', 'useNMF' };
    if( useNMF )
        [D2d, V] = nmf( X_trn', param );
        varlist{4} = 'V';
    else
        D2d = mexTrainDL( X_trn', param );
    end
    
    save( fullfile( SAVEPATH, [SAVEPREFIX,'_dict2d']), varlist{:});
    
end

%%
% figure;
% for i = 1:20
% %     imdisp( reshape( X_trn(i,:), patchSize(1:2)));
%     imdisp( reshape( D2d(:,i), patchSize(1:2)));
%     pause;
%     clf;
% end

