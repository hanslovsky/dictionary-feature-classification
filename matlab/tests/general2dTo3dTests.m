% general2dTo3dTests

%% patches for testing

patch1 = zeros( 9, 9 );
patch1(5,:) = 1;

patch2 = zeros( 9, 9 );
patch2(5,:) = repmat([0.5 2 0.5], 1, 3);

patch3 = zeros( 9, 9 );
patch3(:,5) = 1;

patch4 = zeros( 9, 9 );
patch4(:,5) = repmat([0.5 2 0.5]', 3, 1);

patch5 = (1./3).*ones( 9, 9 );

patch6 = zeros( 9, 9 );
patch6(1,:) = 1;

patch7 = zeros( 9, 9 );
patch7(:,1) = 1;

patch8 = zeros( 9, 9 );
patch8(end,:) = 1;

patch9 = zeros( 9, 9 );
patch9(:,end) = 1;

patch10 = zeros( 9, 9 );
patch10(1:3,1) = 1;

patch11 = zeros( 9, 9 );
patch11(1:3,end) = 1;

patch12 = zeros( 9, 9 );

patch13 = zeros( 9, 9 );
patch13([32 41 50]) = 1;

patch14 = zeros( 9, 9 );
patch14([1 10 19]) = 1;
 
patch15 = zeros( 9, 9 );
patch15([5 14 23]) = 1;

patch16 = zeros( 9, 9 );
patch16(1:27) = 1/3;
%patch16(1:27) = 1;

D = [   patch1(:), patch2(:), patch3(:), patch4(:), patch5(:), patch6(:), patch7(:), ...
        patch8(:), patch9(:), patch10(:), patch11(:), patch12(:), patch13(:), patch14(:),...
        patch15(:), patch16(:)]';

%% generate a known patch from consistent constraints

patchParams = [ 3 1  1; ... %set
                3 4 12; ... %set
                3 7 12; ... %set
                2 1 15; ... %set
                2 4 15; ... %set
                2 7 15; ... %set
                1 1 12; ... %
                1 4 16; ... % questionable
                1 7 12 ];   %

[rc, lc] = rootConstraintFromArray( patchParams ); 

rc.printAsTree()

dsFactor = 3;
patchSize = [ 9 9 ] ;

clear d23;
d23 = Dict2dTo3d( D, patchSize(1), dsFactor );
% d23 = Dict2dTo3dConstr( D, patchSize(1), dsFactor );

[ pv, patch, cmtx, b ] = d23.patchFromParams( lc );

figure; imdisp3d( patch );

% figure; imdisp( cmtx )

%%

rb = cmtx * pv;

%% look at all the constraints

f = 3; 
xyz = 1;
dim = 3;

for i = 1:size( patchParams, 1)
    i
    xyz = patchParams(i,2)
    dim = patchParams(i,1)
    patch = reshape( D(patchParams(i,3),:), [9 9]);
    
    msk = Dict2dTo3d.planeMaskF( [9 9 9], xyz, dim, f );
    pim = Dict2dTo3d.fill3dWith2d( [9 9 9], dim, xyz, f, patch);
    
    figure; imdisp3d( msk );
    figure; imdisp3d( pim );
    
    pause;
    close all;
end

%%

imidx = zeros( [ 9 9 9 ] );
imidx(1:end) = 1:numel(imidx);
figure; imdisp3d(imidx);

%% test serialization

% d23.allSimilaritiesFast();
% [out,it] = d23.build3dPatch();
% save('~/yourface.mat', 'out')
% out
% load('~/yourface.mat')
% out


%% test compiled method

% d23.allSimilaritiesFast();
% obj_fn = d23.save()

% obj_fn = '/groups/saalfeld/home/bogovicj/projects/dictionary/dict2dTo3d/test/myd23.mat';
% save( obj_fn, 'd23' );

%% test 

% bundle.f =  @run_obj_method_dist;
% bundle.args = {d23.obj_fn  'this'  'build3dPatch'  [2]};
% 
% save('/nobackup/saalfeld/john/tmp/job_mytest_1.mat', 'bundle');
% 
% ns = load('/nobackup/saalfeld/john/tmp/job_mytest/job_mytest_1.mat');


%%

% output = run_obj_method_dist(obj_fn, 'd23', ...
%                                      'build3dPatch', 2, ...
%                                      {});
                                 
%%

% vaSingle = {obj_fn, 'this', 'build3dPatch', 2, {}};
% output = run_obj_method_dist( vaSingle{:} )

%% try a larger patch size

patchSize = 15;
dsFactor = 5;

numDictElems = 25;

% build random dictionary
D = rand( numDictElems, patchSize*patchSize );

d23 = Dict2dTo3d( D, patchSize, dsFactor );
dsFactor = 5;

d23.allSimilaritiesFast();
d23.build3dDictionary( 5 );

                                 
%% test distributed 

% [patches, out2 ] = d23.build3dDictionaryDist( 20, 1 );

%%

% use_gpu = 0;
% pe_batch_size = 1;
% run_script = '/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification/bin/my_run_runObjMethod.sh';
% 
% qsubParams1{1} = obj_fn;
% qsubParams2{1} = 'this';
% qsubParams3{1} = 'build3dPatch';
% qsubParams4{1} = 2;
% 
% qsubParams1{2} = obj_fn;
% qsubParams2{2} = 'this';
% qsubParams3{2} = 'build3dPatch';
% qsubParams4{2} = 2;
% 
% vargin = {qsubParams1, qsubParams2, qsubParams3, qsubParams4 };

% varargin = { obj_fn, 'this', 'build3dPatch', 2 };
% varargin = { {obj_fn}, {'this'}, {'build3dPatch'}, {2} };
%                          
% f_out = qsub_dist( @run_obj_method_dist, pe_batch_size, use_gpu, ...
%                            [], [], run_script, ...
%                            vargin{:});
% f_out

%% reshape f_out to something nice 

% dict3dSize = 2;
% sz3d = [9 9 9];
% 
% n = 0;
% patches3d = zeros( dict3dSize, prod( sz3d ));
% for i = 1:size(f_out,1)
%     
%     patchParams = f_out{i}{1}{1};
%     pv = d23.patchFromParams( patchParams );
%     
%     if( ~isempty( pv ))
%         n = n + 1;
%         patches3d(n,:) = pv; 
%     end
% end

%%
% f = 3;
% 
% dim = 3;
% xyz = 1;
% patch = patch1;
% 
% % xyz dim
% pm1 = Dict2dTo3d.planeMaskF( [9 9 9], xyz, dim, f ); 
% imdisp3d( pm1 );
% 
% % dim xyz
% I2 = Dict2dTo3d.fill3dWith2d( [9 9 9], dim, xyz, f, patch );
% figure; imdisp3d( I2 );
% 
% %% what is the cost of this patch? 
% 
% % expect it to be correct?
% 
% p1 = squeeze(patch(5,:,:));
% dim1 = 1;
% xyz1 = 4;
% 
% % p2 = patch5;
% % dim2 = 1;
% % xyz2 = 4;
% 
% p2 = patch5;
% dim2 = 1;
% xyz2 = 4;
% 
% [ sim, x, cmtx, b, pm1, pm2, overlap ] = d23.patchSimilarity( p1, xyz1, dim1, p2, xyz2, dim2 );
% 

