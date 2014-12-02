%% buildDict3d
%
% dbstop if error; run_script('buildDict3d', 'rick-0512-6, K2d=1000, p2d[9], K3d=200, mxiter=1k' );
% dbstop if error; run_script('buildDict3d', 'rick-0512-6, TEST' );

global SAVEPATH
global SAVEPREFIX

%% 

dsFactor = 5;

train_fn = '/groups/saalfeld/saalfeldlab/data/rick/0512-6_Stack4ECS_Segmentation.corrected.0.02-8-5-0.001-0.001-0.001.tif';


% rick_0512-6,  K=500
% Dfile = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0231_buildDict2d/exp0231_dict2d.mat';
% sz = 9;

% rick_0512-6, K=500
% Dfile = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0232_buildDict2d/exp0232_dict2d.mat';
% sz = 15;

% rick_0512-6, K=1000
Dfile = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0235_buildDict2d/exp0235_dict2d.mat'; 
sz = 9;

% rick_0512-6, K=1000
% Dfile = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0234_buildDict2d/exp0234_dict2d.mat';
% sz = 15;


%% d23 params
 
% dictSize3d = 5;
dictSize3d = 200;
% dictSize3d = 1000;

minD3dDiff = -1;
overlappingPatches = 1;
scaleDictElems     = 0;
scaleByOverlap     = 0;
chooseBestAtIter   = 0;
% intXfmModelType    = 'poly1';
intXfmModelType    = '';
pe_batch = 1;
iniOpts = 'hr';
iters = 1000;
d23Record = 0;

% convEps = 1e-5;

convEps = -1;
convIters = 2*iters;

%% build 3d dictionary

load(Dfile);

if( exist( 'Drand', 'var'))
    D2d = Drand';
elseif( exist( 'Dkmeans', 'var'))
    D2d = Dkmeans';
elseif( exist( 'DkmExemplar', 'var'))
    D2d = DkmExemplar';
end

d23 = Dict2dTo3dSampler( D2d', sz, dsFactor, overlappingPatches, scaleByOverlap );
% d23.cleanDict2d();

d23.chooseBestAtIter = chooseBestAtIter;
d23.intXfmModelType = intXfmModelType;
d23.maxIters = iters;
d23.minDictElemDiff = minD3dDiff;
d23.recordParamsOverIters = d23Record;

d23.convEps = convEps;
d23.convIters = convIters;

d23 %#ok<NOPTS>
d23.pc

if( strcmp( iniOpts, 'dict' ))
    
    d23.Dini = D3d';
    d23.iniLocs = d23.pc.locXyzDim2Idx( [3 3 3], [1 4 7]);
    
elseif( strcmp( iniOpts, 'hr' ))
    
    im = readMultiTiff( train_fn );
    train_im = downSampleGaussianImglib( im, [2 2 1], [0.5 0.5 0], [0.5 0.5 0]);
    clear im;
    
    half = (sz-1)./2;
    zHalf = ceil(half/dsFactor);
    
    szIni = [sz sz (1 + 2*zHalf)];
    DiniLR = grabPatchesSimple( train_im, szIni, dictSize3d );
    
    Dini = Dict2dTo3dSampler.upsampleInitialObservations( DiniLR, szIni, [sz sz sz], dsFactor, {} );
    
%     % debug
%     figure; imdisp3d( reshape( DiniLR, szIni ));
%     figure; imdisp3d( reshape( Dini, [sz sz sz] ));
    
%     d23.Dini = Dini;
    myini = mat2cell( Dini, ones(size(Dini,1),1), size(Dini,2) );
    clear train_im;
end

%%

% [patches, out2 ] = d23.build3dDictionaryDist( dictSize3d, pe_batch, iniOpts );
% 
% % save results
% system( sprintf('rm -v %s', d23.obj_fn ));
% save( fullfile( SAVEPATH, [SAVEPREFIX,'_d23']), 'd23', 'out2');

%%

% figure;
% nrmList = zeros( size( D2d, 2 ), 1 );
% 
% % for i = 1:d23.numDict
% for i = 1:size( D2d, 2 )
%     
%     fprintf( 'i: %d\n', i );
%     
% %     imdisp3d( reshape( d23.D2d( i, :), [sz sz] ));
% 
%     imdisp( reshape( D2d( :, i ), [ sz sz ] ));
%     pause;
%     clf;
% 
% 	nrmList(i) = norm( D2d( :, i ));
% 
% end
% 
% hist( nrmList );

%% test fit dists
myini
[patches, out2 ] = d23.fitParamsToHR_dist( myini );

% save results
system( sprintf('rm -v %s', d23.obj_fn ));
save( fullfile( SAVEPATH, [SAVEPREFIX,'_d23']), 'd23', 'patches');



%%

% iniParams = randi( d23.numDict, d23.pc.numLocs, 1);
% 
% %% test out dict building
% 
% [ patchParams, pv, iteration, costs, patchesByIter, otherParams ] =  d23.build3dPatch( iniParams );
% 
% iterNumUpdated = find( costs > -1 );
% plot( iterNumUpdated, costs( iterNumUpdated ));
% 
% imdisp3d( reshape( pv, d23.sz3d ));

