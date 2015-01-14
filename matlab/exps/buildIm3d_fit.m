% buildIm3d_fit
% dbstop if error; run_script('buildIm3d_fit', 'fitting for testing');
% dbstop if error; buildIm3d_fit;

global SAVEPATH;
global SAVEPREFIX;

%% load image data

im_codeString = 'microtubulesHR';
[ im_fn, dsFactor, im, mskFun  ] = imageData( im_codeString, 1 );
im = mean( im, 4 );

%%

sz = [ 5 5 5 ];
sz2d = [3 3 ];
f = 2;

scaleIni           = 1;
overlappingPatches = 1;
scaleByOverlap     = 0;
doInvCmtx          = 0;

%%
rng( 42 );

N = 1;
DiniHR = grabPatchesSimple( im, sz, N, [], mskFun(im) );
if( scaleIni )
    DiniHR  = bsxfun( @rdivide, DiniHR , sqrt(sum(DiniHR .^2,2)));
end

% isequal( DiniHR, DiniHR2 )
%%
% imdisp3d( reshape( DiniHR, sz ));
% DiniHR2 = DiniHR;

%%

M = 50;
D2d = rand( M, prod( sz2d ));
if( scaleIni )
    D2d = bsxfun( @rdivide, D2d, sqrt(sum(D2d.^2,2)));
end

d23g = Dict2dTo3dSamplerGen( D2d, sz2d, sz, f, overlappingPatches, scaleByOverlap, doInvCmtx );

% i = 1;
% [ idx, curdist, model ] = d23g.fitIdxAndModel(  i, DiniHR, 1 );

[ patchParams, dists, modelList ] = d23g.fidIdxAndModel_dist( DiniHR, 1 );

system( sprintf('rm -v %s', d23g.obj_fn ));
save( fullfile( SAVEPATH, [SAVEPREFIX,'_d23']), 'd23g', 'patchParams', 'dists', 'DiniHR');

