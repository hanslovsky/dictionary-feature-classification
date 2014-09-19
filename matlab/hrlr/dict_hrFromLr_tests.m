%% dict_hrFromLr_tests

%%

sz = 3;
dorep = 1;

rng = linspace( 0, 1, sz );
[x,y] = ndgrid( rng, rng );
D = [ x(:), ones(numel(x),1), cos(pi.*x(:)), sin(pi.*x(:)) ]';

%% include rotations of each element
rxi = cell2mat( xfmToIdx( squareSymmetry(), [sz sz], 0 ));
rxi = rxi';

nRxi = size(rxi,1);
nD   = size(D, 1);
Dbig = zeros( nD.*nRxi, size(D,2));

k = 1;
numTot = 0;
for i = 1:nD
    dd = D(i,:);
    dd = dd(rxi);

    pd = pdist( dd );
    dd
    [n,m] = ind2subTri( size(dd,1), find( pd < 0.1 ));
%     [n ; m] 
    
    [~, ki ] = dissimilarVecs( dd, 0.1 );
%     dsv
    ki
    
    num2add = length( ki );
    
    Dbig( k:k+num2add-1, : ) = dd(ki,:); 
    Dbig

    k = k + num2add;
end

Dbig = Dbig( 1:(k-1), :);
Dbig

%%
clear tid2;

tid2 = Tid_from2D( Dbig', 3 );
tid2.genVectorTransformations();
tid2.makeDictRotInv( 2 );

Dbig2 = tid2.D';
% return;

%%
N = size(D,1);

i = reshape( 1:numel(x), [sz sz] );
irepz = repmat( i, [1 1 sz]);
irot = repmat( permute( i, [ 3 1 2]), [sz 1 1 ]); 

jrep = irepz(:);
jrot = irot(:);

%% 

if( dorep )
    [ii,jj] = ndgrid( 1:N, 1:N );
    ii = ii(:);
    jj = jj(:);
else
    ii = zeros( N.*(N-1)./2, 1);
    jj = zeros( N.*(N-1)./2, 1);
    n = 1;
    for a = 1:N
        for b = a+1:N
            ii(n) = a;
            jj(n) = b; 
            n = n + 1;
        end
    end
end

Dcomb=D(ii,jrep).*D(jj,jrot);
Dcomb

%%
clear tid3;
tid3 = Tid_from2D( Dcomb', 3 );
tid3.genVectorTransformations();
tid3.makeDictRotInv( 3 );

Dcomb = tid3.D';

%% visualize

for n = 1:size(Dcomb,1)
    
    figure;
    imdisp( permute(reshape( Dcomb(n,:), repmat(sz, 1, 3)), [1 2 4 3]), 'border', 0.1 );
   
%     pause;
%     close all;
end


%% another old test
%
% a = rand(sz);
% b = rand(sz);
% 
% 
% a(jrep).*b(jrot);
% 
% ar = a(:);
% br = b(:);
% 
% ar(jrep).*b(jrot);

%% old 1d examples
% f = linspace( 0, 1, 11);
% D = [ f; 0.5.*ones( 1, 11); cos(pi.*f); sin(pi.*f); ]; 

%%
%aa = repmat( a, [1 1 sz] )
%bb = repmat( permute( b, [ 3 1 2]), [sz 1 1]);
%
%cc = aa .* bb;
%cc(:)

%% possibly faster way to compute the above

%[x,y]= ndgrid(1:2, 1:2)
%xrepz = repmat( x, [1 1 2]);
%yrepz = repmat( y, [1 1 2]);
%
%xrot = repmat( permute( y, [ 3 1 2]), [2 1 1]);
%yrot = repmat( permute( y, [ 3 1 2]), [2 1 1]);
