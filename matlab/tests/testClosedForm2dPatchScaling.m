%% testClosedForm2dPatchScaling

%% its a quadratic program

Q = [ 1 1 0 0; 0 0 1 1; 1 0 1 0; 0 1 0 1 ];
Q

costFun = @(x,b)( norm(Q * x(1:4) - x([5 5 6 6]) .* b') );

QtQ = Q'* Q

genFun = @(x, beta, b)( diag( [beta(1) beta(1) beta(2) beta(2)]) * Q * x );
obfun = @(x, b)( diag( [beta(1) beta(1) beta(2) beta(2)]) * Q * x );

%% test 3d

clear pc;
pc = PatchConstraints( 2, 2, 0, 0 );
pc.buildCmtx();
Q = pc.cmtx;

costFun = @(x,b)( norm(pc.cmtx * x(pc.numLocs + 1 : end) ...
                  - x(reshape(repmat(1:pc.numLocs, prod(pc.sz2d), 1),1,[])) .* b) );

% b = repmat([ 4 6 3 7 ], 1, 3 );
b = 2 ./prod(pc.sz3d) .* ones( pc.numConstraints, 1 );
b(1:4) = b(1:4) .* 2;
b(5:8) = b(5:8) .* 4;
b(9:end) = b(9:end) .* 6;

Xtrue = ones( pc.numLocs + prod(pc.sz3d), 1 );
Xtrue( pc.numLocs + 1 : end ) = 1./prod(pc.sz3d);

Xtrue( 1:3 ) = [ 1./2 1./4 1./6 ]

costFun( Xtrue, b )

%%
[ H3d, f3d, A3d, b3d, Aeq3d, beq3d ] = pc.buildQuadProg( b );
H3d
A3d
Aeq3d

% [Xest,fval,exitflag] = quadprog( H3d, f3d, A3d, b3d, Aeq3d, beq3d,...
%                           [], [], [], 'interior-point-convex' );
                      
[Xest,fval,exitflag] = quadprog( H3d, f3d, [], [], Aeq3d, beq3d,...
                          [], [], [], 'interior-point-convex' );
Xest
fval
exitflag

fprintf( 'Xest: \n');
fprintf( 'Do the inequality constraints hold: %d\n', all(A3d * Xest < b3d ));
fprintf( 'Do the equality constraints hold: %d\n', all(abs(Aeq3d * Xest) - beq3d ) < 0.0001);
fprintf( 'Real Cost: %d\n', costFun( Xest, b) );
fprintf( 'Quad Cost: %d\n', Xest' * H3d * Xest );
fprintf( '\n\n' );

fprintf( 'Xtrue: \n');
fprintf( 'Do the inequality constraints hold: %d\n', all(A3d * Xtrue < b3d ));
fprintf( 'Do the equality constraints hold: %d\n', all(abs(Aeq3d * Xtrue) - beq3d ) < 0.0001);
fprintf( 'Real Cost: %d\n', costFun( Xtrue, b ));
fprintf( 'Quad Cost: %d\n', Xtrue' * H3d * Xtrue );

%%

% Atemplate = [ 1 1 -1; 1 1 -1; -1 -1 -2 ];
% mbTempalte = ( Atemplate == -1 );  % minus b
% bsTempalte = ( Atemplate == -2 );  % b squared
% H = zeros( 6, 6 );
% 
% % x1 x2 - b1
% % x3 x4 - b2
% % |  | 
% % b3 b4
% 
% b = [ 4 6 3 7];
% alphaTrue = [ 0.5 2 ];
% betaTrue  = 1./alphaTrue;
% b = b .* betaTrue([ 1 1 2 2])
% 
% % b=1
% Atmp = Atemplate;
% Atmp( mbTempalte ) = -b(1); 
% Atmp( bsTempalte ) = b(1)*b(1); 
% 
% rng = [ 1 2 5 ]; % x1, x2, alpha1(5)
% H(rng,rng) = Atmp;
% 
% % b=2
% Atmp = Atemplate;
% Atmp( mbTempalte ) = -b(2); 
% Atmp( bsTempalte ) = b(2)*b(2); 
% 
% rng = [ 3 4 5 ]; % x1, x2, alpha1(5)
% H(rng,rng) = H(rng,rng) + Atmp;
% 
% % b=3
% Atmp = Atemplate;
% Atmp( mbTempalte ) = -b(3); 
% Atmp( bsTempalte ) = b(3)*b(3); 
% 
% rng = [ 1 3 6 ]; % x1, x2, alpha1(5)
% H(rng,rng) = H(rng,rng) + Atmp;
% 
% % b=4
% Atmp = Atemplate;
% Atmp( mbTempalte ) = -b(4); 
% Atmp( bsTempalte ) = b(4)*b(4); 
% 
% rng = [ 2 4 6 ]; % x1, x2, alpha1(5)
% H(rng,rng) = H(rng,rng) + Atmp
% 
% %%
% f = [];
% 
% % constrain scales to be positive
% A = -[ 0 0 0 0 1 0; 0 0 0 0 0 1 ];
% bieq =  [ 0; 0 ];
% 
% % Ceq = [ones( 1, 4 ) 0 0];
% % beq  = 10;
% Ceq = [];
% beq = [];
% 
% [Xest,fval,exitflag] = quadprog( H, f, A, bieq, Ceq, beq );
% Xest
% fval
% exitflag
% 
% fprintf( 'Xest: \n');
% fprintf( 'Do the inequality constraints hold: %d\n', all(A * Xest < bieq ));
% fprintf( 'Do the equality constraints hold: %d\n', all((Ceq * Xest) == beq ));
% fprintf( 'Real Cost: %d\n', costFun( Xest, b) );
% fprintf( 'Quad Cost: %d\n', Xest' * H * Xest );
% fprintf( '\n\n' );
% 
% %% check xtrue
% 
% Xtrue = [ 1 3 2 4 alphaTrue]';
% 
% fprintf( 'Xtrue: \n');
% fprintf( 'Do the inequality constraints hold: %d\n', all(A * Xtrue < bieq ));
% fprintf( 'Do the equality constraints hold: %d\n', all(Ceq * Xtrue == beq ));
% fprintf( 'Real Cost: %d\n', costFun( Xtrue, b ));
% fprintf( 'Quad Cost: %d\n', Xtrue' * H * Xtrue );


