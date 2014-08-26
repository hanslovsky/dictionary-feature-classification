% testDictSymSample


patchSize = [3 3 1];
%D(:,1) = reshape(eye(3),[],1);
%D(:,2) = reshape(ones(3,3),[],1);
%D(:,3) = reshape([1 1 1; 0 0 0; 0 0 0],[],1);
D = DictionarySymSampler.standardDictionary2d();

[xs,ys,zs] = ndgrid( -1:1, -1:1, 0 );
shiftRng = [ xs(:), ys(:), zs(:) ];

dss = DictionarySymSampler( D, [], patchSize, [], [], shiftRng );

samples = dss.sample( 1 )
