%% test_hrlrTheory 

%% pseudo-inverse of constraint matrices

Cr = [ 1 1 0 0; 0 0 1 1 ]
Cc = [ 1 0 1 0; 0 1 0 1 ]

Cri = pinv( Cr )
Cci = pinv( Cc )

% Ccr = Cr + Cc
% Ccri = pinv( Ccr )

Ccr = [Cr; Cc]
Ccri = pinv( Ccr )


%% inconsistent case a,b

A = zeros( 6 + 9, 18);
A( 1, 1:3 ) = 1;
A( 2, 4:6 ) = 1;
A( 3, 7:9 ) = 1;

A( 4:6, 10:end ) = repmat( eye(3), 1, 3 );

A( 7:end, 1:9) = eye( 9 );
A( 7:end, 10:end) = -1 .* eye( 9 );

b = zeros( size(A,1), 1 );
b(1:6) = [ 1 13 1 5 6 5 ]'

x = pinv(A) * b


%% consistent case x

A = zeros(6, 9);
A( 1, 1:3 ) = 1;
A( 2, 4:6 ) = 1;
A( 3, 7:9 ) = 1;
A( 4:6, : ) = repmat( eye(3), 1, 3 );

b = [ 1 13 1 5 5 5]'

x = pinv(A) * b

%% consistent case a,b

A = zeros( 6 + 9, 18);
A( 1, 1:3 ) = 1;
A( 2, 4:6 ) = 1;
A( 3, 7:9 ) = 1;

A( 4:6, 10:end ) = repmat( eye(3), 1, 3 );

A( 7:end, 1:9) = eye( 9 );
A( 7:end, 10:end) = -1 .* eye( 9 );

b = zeros( size(A,1), 1 );
b(1:6) = [ 1 13 1 5 5 5 ]'

x = pinv(A) * b

%% consistent case a,b quad prog 2 x 2

% H = [  ...
%      1  0  0  0 -1  0  0  0 ; ...
%     -1  0  0  0  1  0  0  0 ; ...
%      0  1  0  0  0 -1  0  0 ; ...
%      0 -1  0  0  0  1  0  0 ; ...
%      0  0  1  0  0  0 -1  0 ; ...
%      0  0 -1  0  0  0  1  0 ; ...
%      0  0  0  1  0  0  0 -1 ; ...
%      0  0  0 -1  0  0  0  1 ]

% H = [  ...
%      1  0  0  0 -1  0  0  0 ; ...
%      0  1  0  0  0 -1  0  0 ; ...
%      0  0  1  0  0  0 -1  0 ; ...
%      0  0  0  1  0  0  0 -1 ; ...
%     -1  0  0  0  1  0  0  0 ; ...
%      0 -1  0  0  0  1  0  0 ; ...
%      0  0 -1  0  0  0  1  0 ; ...
%      0  0  0 -1  0  0  0  1 ]
 
H = [eye(4) -eye(4); -eye(4) eye(4)]
 
f = zeros( 8, 1 )

Aeq = [ ...
    1 1 0 0 0 0 0 0; ...
    0 0 1 1 0 0 0 0; ...
    0 0 0 0 1 0 1 0; ...
    0 0 0 0 0 1 0 1  ]

% beq = [ 1.5; 1.5 ; 1; 2 ]   % consistent constraints
beq = [ 1.5; 1.5 ; 1; 1 ]   % inconsistent constraints
% beq = [ 2; 2; 1; 1 ]

% lb = zeros( 8, 1);
% ub = ones( 8, 1);

lb = [];
ub = [];

[X,fval,exitflag,output,lambda] = quadprog(H, f, [], [], Aeq, beq, lb, ub);
exitflag
X

fval

% err= sqrt(norm(beq(1:2) - beq(3:4)));
% err

% 0.5 .* norm( null(A)' * beq ).^2
% norm( pinv(A) * beq )

%% inconsistent case x

A = [ 1 1 0 0; 0 0 1 1; 1 0 1 0; 0 1 0 1 ];

x = pinv(A) * beq;
x
mean([X(1:4) X(5:8)],2)

xb = [ x; x ]
% norm( xb' * H * xb )

norm( A*x - beq )
norm(X(1:4)-X(5:end))

%%
nh = size(   H, 1 );
na = size( Aeq, 1 );

M = zeros( size(H,1) + size(Aeq,1) );
M( 1:nh, 1:nh ) = H;
M( 1:nh, nh+1:end ) = Aeq';
M( nh+1:end, 1:nh ) = Aeq;

bb = zeros( size(M,1), 1);
bb(nh+1:end) = beq

x2 = pinv(M) * bb 

est = M * x2

%% see how close the min of the unconstrained min gets to the constrained

x = pinv(A) * beq


%% test constrained vs unconstained on tons of examples

% setup
f = zeros( 8, 1 );
H = [eye(4) -eye(4); -eye(4) eye(4)];

Aeq = [ ...
    1 1 0 0 0 0 0 0; ...
    0 0 1 1 0 0 0 0; ...
    0 0 0 0 1 0 1 0; ...
    0 0 0 0 0 1 0 1  ]

A = [ 1 1 0 0; 0 0 1 1; 1 0 1 0; 0 1 0 1 ];


lb = [];
ub = [];

N = 1000;
mul = 1;

wrongCount = 0 ;
rightCount = 0 ;

opts = optimoptions('quadprog');
opts.Display   = 'off';
opts.Algorithm = 'interior-point-convex';

for i = 1:N
    
    if( mod(i,10)==0 )
        fprintf('%d\n',i);
    end
    
    beq = mul .* rand( 4, 1);
    
    [ xc,fval,exitflag,output,lambda] = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], opts);
    
    xuc =  pinv(A) * beq;
    
    simUC = norm( A*xuc - beq );
    simC = norm( xc(1:4)-xc(5:end));
   
    if( abs( simUC - simC ) > 0.000001 )
        wrongCount = wrongCount + 1;
    else
        rightCount = rightCount + 1;
    end
    
end

fprintf('results were the same for %d of %d cases\n', rightCount, N )

%% consistent case a,b quad prog 3 x 3

N = 27;
 
H = [eye(N) -eye(N); -eye(N) eye(N)]
 
f = zeros( 2*N, 1 )

Aeq = zeros( 2*9, 2*27 );
k = 1;
for i = 1:9
   
   Aeq( i, k:k+2 ) = 1;
   k = k + 3;
   
end

l = 1;
for i = 1:3:9
    
%     Aeq( i:i+2, (27+l):(27+l+8) ) = repmat(eye(3), 1, 3 );
%     Aeq( (9+i):(9+i+2), (l):(l+8) ) = repmat(eye(3), 1, 3 );
    Aeq( (9+i):(9+i+2), (27+l):(27+l+8) ) = repmat(eye(3), 1, 3 );

	l = l + 9;
end

Aeq
imagesc( Aeq )

%%

% beq = [ones(9,1); ones(9,1)]    % consistent constraints
% beq = [ones(9,1); 1.5.*ones(9,1)];   % inconsistent constraints
beq = [ones(9,1); 0; 0; 0; 1; 1; 1; 0; 0; 0];   % inconsistent constraints 2

lb = [];
ub = [];

[X,fval,exitflag,output,lambda] = quadprog(H, f, [], [], Aeq, beq, lb, ub);
exitflag
X

fval


%%

A = zeros( 2*9, 27 );
k = 1;
for i = 1:9
   
   A( i, k:k+2 ) = 1;
   k = k + 3;
   
end

l = 1;
for i = 1:3:9
    
    A( (9+i):(9+i+2), (l):(l+8) ) = repmat(eye(3), 1, 3 );

	l = l + 9;
end

% imagesc(A)

x = pinv(A)*beq;

[ mean( [ X(1:27) X(28:end) ], 2) x ]
% x

%% test downsampling factor of 10

simUCLast = [];
simCLast  = [];

sz = 30;
f  = 10;

sz3 = [sz sz sz];


xyz1 = 1;
xyz2 = 1;

dim1 = 1;
dim2 = 2;

pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, dim1, f );
pm2 = Dict2dTo3d.planeMaskF( sz3, xyz2, dim2, f );
overlap = (pm1 > 0) & (pm2 > 0);

numSame       = 0;
numConsistent = 0;
numTransitive = 0;

N = 251;
Xuc = zeros( N, 1 );
Yc = zeros( N, 1 );

for n = 1:N
    
    patch1 = rand( sz, sz );
    patch2 = rand( sz, sz );
    patch3 = rand( sz, sz );

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % similarity with unconstrained min
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [cmtx1, b1] = Dict2dTo3d.contraintsFromMask( patch1, overlap, pm1 );
    [cmtx2, b2] = Dict2dTo3d.contraintsFromMask( patch2, overlap, pm2 );
    [cmtx3, b3] = Dict2dTo3d.contraintsFromMask( patch3, overlap, pm2 );
    
    cmtx = [ cmtx1; cmtx2 ] ;
    b    = [ b1; b2 ];
    
    
    x = pinv(cmtx) * b;
    sim = sum((cmtx*x - b).^2);
    simUC = norm( cmtx*x - b );
    
    bCompare   = [ b1; b3 ];
    xCompare = pinv(cmtx) * bCompare;
    simUCCompare = norm( cmtx*xCompare - bCompare );
    %
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % similarity with constrained min
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    N = nnz( overlap );
    numConstraintsPerPatch = length(unique( pm1(overlap) ));
    
    H = [eye(N) -eye(N); -eye(N) eye(N)];
    f = zeros( 2*N, 1 );
    
    Aeq = [  cmtx1, zeros( numConstraintsPerPatch, N ); ...
        zeros( numConstraintsPerPatch, N ), cmtx2 ];
    
    lb = [];
    ub = [];
    
    opts = optimoptions('quadprog');
    ops.Display = 'off';
    opts.Algorithm = 'interior-point-convex';
    
    [ xc,fval,exitflag,output,lambda] = quadprog( H, f, [], [], Aeq, b, lb, ub, [], opts );
    
    
    % simC = sum((xc(1:N) - xc(N+1:end)).^2)
    simC = norm( xc(1:N)-xc(N+1:end));
    
    [ xcCompare ] = quadprog( H, f, [], [], Aeq, bCompare, lb, ub, [], opts );
    simCCompare = norm( xcCompare(1:N) - xcCompare(N+1:end));
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % are they the same ?
%     simC
%     simUC

    Xuc( n ) = simUC;
    Yc ( n ) = simC;

    if ( (simC > simCCompare && simUC > simUCCompare) || ...
         (simC < simCCompare && simUC < simUCCompare) )
       
        numTransitive = numTransitive + 1;
    end
       
    if( n > 1 )
        if( simC > simUC && simCLast > simUCLast )
            numConsistent = numConsistent + 1;
        elseif( simC < simUC && simCLast < simUCLast )
            numConsistent = numConsistent + 1;
        end
    end
    
    if( abs(simC - simUC ) < 0.0001 )
        numSame = numSame + 1;
    end
        
    simCLast  = simC;
    simUCLast = simUC;
    
    
end

numSame
numConsistent
numTransitive

%%

plot( Xuc, Yc, '.' );
xlabel('unconstrained min');
ylabel('constrained min');

xmin = min( Xuc );
xmax = max( Xuc );

hold on;
plot( [xmin xmax], [xmin xmax], '--k');

% export_fig('~/unconstrainedVsConstrainedPatchSim_s9f3.png', '-m2' )
% export_fig('~/unconstrainedVsConstrainedPatchSim_s15f5.png', '-m2' )
% export_fig('~/unconstrainedVsConstrainedPatchSim_s30f10.png', '-m2' )

%%

figure;  
imdisp( patch1 );
title('patch1')

figure;  
imdisp( patch2 );
title('patch2');

figure; 
imdisp( permute(reshape( x, sz3 ),[1 2 4 3]), 'border', 0.02 );
title('opt UC');

xmn = mean( [xc(1:N), xc(N+1:end)],2); 
figure; 
imdisp( permute(reshape( xmn, sz3 ),[1 2 4 3]), 'border', 0.02 );
title('opt C');
