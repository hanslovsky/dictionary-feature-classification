%% simBuild3dDictFrom2d
% dbstop if error; run_script('simBuild3dDictFrom2d', 'test bubble data, f=1');
% dbstop if error; run_script('simBuild3dDictFrom2d', 'test bubble data, NMF, N=1m, f=3');
% dbstop if error; run_script('simBuild3dDictFrom2d', 'test bubble data, N=1m, K=500, f=3');
% dbstop if error; run_script('simBuild3dDictFrom2d', 'test bubble data, N=1m, K=500, f=5');
%
% dbstop if error; run_script('simBuild3dDictFrom2d', 'test bubble data, N=1m, K=500, f=3, noOL, D3dini, cleaned 2d dict');
% dbstop if error; run_script('simBuild3dDictFrom2d', 'test bubble data, N=1m, K=500, f=5, D3dini');

global SAVEPATH
global SAVEPREFIX

%% 
srcdir = '/groups/saalfeld/home/bogovicj/projects/dictionary/dict2dTo3d/sim/bubbles';

train_fn = fullfile( srcdir, 'Bubbles_r3_b1_snr10_1.tif');
test_fn  = fullfile( srcdir, 'Bubbles_r3_b1_snr10_2.tif');

%% params

N = 1000000;
dsFactor = 3;
patchSize = [9 9 dsFactor];
patchSize3d = [9 9 9];

% dsFactor = 5;
% patchSize = [15 15 dsFactor];
% patchSize3d = [15 15 15];

% dictionary building parameters
param.K = 500;  % dictionary size
param.lambda=0.1;
param.numThreads=1; % number of threads
param.verbose=1;
param.iter = 500;  % num iterations.

dictSize3d = 200;
minD3dDiff = -1;
overlappingPatches = 0;
pe_batch = 1;

% Dfile = [];

% 9x9 dict
Dfile = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0172_simBuild3dDictFrom2d/exp0172_dict2d.mat';

% 15x15 dict
% Dfile = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0174_simBuild3dDictFrom2d/exp0174_dict2d.mat';


% D3dfile = [];

% 9x9x3 dict
D3dfile = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0176_simBuild3dDictFrom2d/exp0176_dict3d.mat';



useNMF = 0;

% iniOpts = 'rand';
iniOpts = 'dict';

%% build 2d dictionary

if( isempty( Dfile ))
    
    train_im = readMultiTiff( train_fn );
    X_trn = grabPatchesSimple( train_im, patchSize, N );
    
    if( dsFactor > 1 )
        downsampler = Tid.getDownsampler3dzRs();
        X_trn = [downsampler( X_trn', patchSize, dsFactor )]';
    end
    
    varlist = { 'D2d', 'param', 'useNMF' };
    if( useNMF )
        [D2d, V] = nmf( X_trn', param );
        varlist{4} = 'V';
    else
        D2d = mexTrainDL( X_trn', param );
    end
    
    save( fullfile( SAVEPATH, [SAVEPREFIX,'_dict2d']), varlist{:});
    
else
	fprintf('loading dictionary from file: %s\n', Dfile);
    load( Dfile );
end

%% build 3d dictionary if requested

if( strcmp( iniOpts, 'dict' ))
    if( isempty( D3dfile ))
        fprintf( 'learning 3d dictionary\n' );
        if( ~exist( 'train_im', 'var'))
            train_im = readMultiTiff( train_fn );
        end
        
        X_trn_3d = grabPatchesSimple( train_im, patchSize3d, N );
        if( dsFactor > 1 )
            downsampler = Tid.getDownsampler3dzRs();
            X_trn_3d = [downsampler( X_trn_3d', patchSize3d, dsFactor )]';
        end
        
        varlist = { 'D3d', 'param', 'useNMF' };
        if( useNMF )
            [D3d, V3d] = nmf( X_trn_3d', param );
            varlist{4} = 'V3d';
        else
            D3d = mexTrainDL( X_trn_3d', param );
        end
        size( D3d )
        save( fullfile( SAVEPATH, [SAVEPREFIX,'_dict3d']), varlist{:});
    else
        fprintf( 'loading 3d dictionary from file: %s\n', D3dfile );
        load( D3dfile );
    end
end

%% vis
% if( useNMF )
% 	D2d_nz = D2d(:, ((sum( D2d, 1 ) > 0 )));
% else
%     D2d_nz = D2d;
% end
% 
%
% for i = 1:size(D2d_nz,2)
%     i
%     imdisp( reshape( D2d_nz(:,i), patchSize(1:2)));
%     pause;
%     close all;
%     
% end

%% build 3d

d23 = Dict2dTo3dSampler( D2d', patchSize(1), dsFactor, overlappingPatches );
d23
d23.pc

d23.minDictElemDiff = minD3dDiff;
if( strcmp( iniOpts, 'dict' ))
    d23.Dini = D3d';
    d23.iniLocs = d23.pc.locXyzDim2Idx( [3 3 3], [1 4 7]);
end

% [ D3d ] = d23.build3dDictionary( dictSize3d );
[patches, out2 ] = d23.build3dDictionaryDist( dictSize3d, pe_batch, iniOpts );

system( sprintf('rm -v %s', d23.obj_fn ));

save( fullfile( SAVEPATH, [SAVEPREFIX,'_d23']), 'd23', 'param');
