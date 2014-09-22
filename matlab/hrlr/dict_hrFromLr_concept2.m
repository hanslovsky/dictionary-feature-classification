%% dict_hrFromLr_tests

%%

sz = 3;
dorep = 1;

rng = linspace( 0, 1, sz );
[x,y] = ndgrid( rng, rng );
D = [ x(:), ones(numel(x),1), cos(pi.*x(:)), sin(pi.*x(:)) ]';

%% include rotations of each element

clear tid;

tid = Tid_from2D( D', 3 );
tid.buildDictionary();

D = tid.D3d';


%% visualize

% for n = 1:size(Dbig,1)
for n = 1:size(D,1)
    
    figure;
    imdisp( permute(reshape( D(n,:), tid.patchSize3d), [1 2 4 3]), 'border', 0.1 );
%     imdisp( reshape( Dbig(n,:), repmat(sz, 1, 2)) );
   
%     pause;
%     close all;
end
