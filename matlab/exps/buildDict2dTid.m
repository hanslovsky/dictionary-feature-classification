%% buildDict2dTid
% dbstop if error; run_script('buildDict2dTid', 'test');
% dbstop if error; run_script('buildDict2dTid', 'ds_medulla, build D2d, norm obs, ps[9x9], N=4m, K=4k, iters=4k, lam=0.4');
% dbstop if error; run_script('buildDict2dTid', 'ds_medulla, build D2d, norm obs, ps[9x9], N=4m, K=4k, iters=4k, lam=0.4, bigIters=3');

global SAVEPATH
global SAVEPREFIX

%% 

% srcdir = '/groups/saalfeld/home/bogovicj/projects/dictionary/dict2dTo3d/sim/bubbles';
% train_fn = fullfile( srcdir, 'Bubbles_r3_b1_snr10_1.tif');
% test_fn  = fullfile( srcdir, 'Bubbles_r3_b1_snr10_2.tif');
% mskFun = @(im)( im > 0 );
% dsFor2dDict = 1;

% train_fn = fullfile( '/groups/saalfeld/saalfeldlab/data/virtual/microtubules/1_5', '*.png');
% test_fn  = fullfile( '/groups/saalfeld/saalfeldlab/data/virtual/microtubules/1_1', '*.png');
% mskFun = @(im)(im < 255 );
% dsFor2dDict = 0;

% train_fn = '/groups/saalfeld/saalfeldlab/data/rick/0512-6_Stack4ECS_Segmentation.corrected.0.02-8-5-0.001-0.001-0.001.tif';
% im = readMultiTiff( train_fn );
% train_im = downSampleGaussianImglib( im, [2 2 1], [0.5 0.5 0], [0.5 0.5 0]);
% clear im;
% mskFun = @(im)( true(size(im)) );
% dsFor2dDict = 0;

% imcode = 'rick_0512';
% imcode = 'microtubulesD5';
imcode = 'ds_medulla_train';

[ train_fn, dsFactor_im, train_im, mskFun  ] = imageData( imcode, 1 );

%% params

N = 4e6;
% N = 1e3;

% dsFactor = 3;
% patchSize = [5 5 dsFactor];

% dsFactor = 5;
% patchSize = [9 9 dsFactor];
% nSlices = 1;

dsFactor = 3;
patchSize = [9 9 dsFactor];
dsFor2dDict = 1;
nSlices = 3;

% dsFactor = 5;
% dsFor2dDict = 0;
% patchSize = [15 15 dsFactor];
% nSlices = 1;

% % test parameters
% param.K = 400;  % dictionary size
% param.lambda = 0.4;
% param.numThreads = 15; % number of threads
% param.verbose = 1;
% param.iter = 200;  % num iterations.

% % dictionary building parameters
param.K = 4000;  % dictionary size
param.lambda = 0.4;
param.numThreads = 4; % number of threads
param.verbose = 0;
param.iter = 200;  % num iterations.
tid_big_iters = 3;

normData = 1;

%%

% iniDictCode = '';
% iniDict = [];

iniDictCode = 'ds_medulla,2d,sz9,K=4k,N=4m,lam=p4';
[ iniDict, sz, dsFactor ] = dictionaryData( iniDictCode );

iniDict = iniDict( 1:param.K, : );

%% get the data

if( ~exist( 'Dfile', 'var') || isempty( Dfile ))
    
    if( ~exist( 'train_im', 'var') )
        if( strcmp( train_fn(end-2:end), 'png' ))
            train_im = readImageStack( train_fn );
        else
            train_im = readMultiTiff( train_fn );
        end
    end
    
    if( dsFor2dDict && (dsFactor > 1 ))
        [X_trn, coords] = grabPatchesSimple( train_im, patchSize, N, [], mskFun( train_im ) );
        fprintf('downsampling observations before building 2d dict\n');
        downsampler = Tid.getDownsampler3dzRs();
        X_trn = [downsampler( X_trn', patchSize, dsFactor )]';
    else
        fprintf('taking observations directly before building 2d dict \n');
        [X_trn, coords] = grabPatchesSimple( train_im, [patchSize(1:2) nSlices ], N, [], mskFun( train_im ) );
        size( X_trn )
    end
    
    if( normData )
        fprintf('normalizing observations\n');
        X_trn = bsxfun( @rdivide, X_trn, sqrt(sum(X_trn.^2,2)));
    end
    
     save( fullfile( SAVEPATH, [SAVEPREFIX,'_coords']), 'coords');
    
    %% build the dictionary
    
    tid = Tid( X_trn, patchSize(1:2), patchSize(1:2)-2, param );
    tid.setComparator( 'ncc' );
    tid.bigIters = tid_big_iters;
    tid.nJobs = 20;
    tid.distTol = 0.02
    
    varlist = { 'D2d', 'param', 'tid', 'D2dList', 'simRiList' };
    [D2d, D2dList, simRiList] = tid.buildTranslationInvariantDictionary( iniDict, 1 );

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

