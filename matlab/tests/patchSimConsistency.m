%% patchSimConsistency

%% setup
Dfile = [];
Dfile = '/groups/saalfeld/home/bogovicj/reseach/exp/saved_exp/exp0163_segWith2dTo3dDict/exp0163_dict.mat';

% d23File = '/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification/matlab/tests/d23Test.mat';
% d23File = '/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification/matlab/tests/d23Test_small.mat';
d23File = '/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification/matlab/tests/d23Test_mid.mat';

if( ~isempty( Dfile ))
    load( Dfile );
end

D = D(:,1:500);

if( exist( d23File, 'file'))
    load( d23File );
else
    sz = 9;
    f = 3;
    d23 = Dict2dTo3d( D', sz, f );
    d23.allSimilaritiesFast();
    save(d23File,'d23');
end

%%
clear d23c thisobj;
d23c = Dict2dTo3dConsistencyTester([],[],[]);
d23c.clone( d23 );
d23c.NbestTest = d23c.numDict;
d23c

%%

d23FastList = cell( numInis, 1 );
for n = 1:numInis
    thisobj = d23.copy();
    thisobj.build3dPatch();
    d23FastList{n} = thisobj;
end

%%
targetDepth = 8;

nodeList = cell( numInis, 1 );
bestList = cell( numInis, 1 );
for n = 1:numInis
    thisobj = d23FastList{n};
    
    bestList{n} = thisobj.currentBestPatchParams();
    
    it = thisobj.p2dFill3d.getSet().iterator();
    listout = java.util.ArrayList();
    while( it.hasNext())
        listout.add( it.next());
    end
    
    % shuffle the nodes 
    java.util.Collections.shuffle( listout );

    % pick a matching node matching the target depth
    it = listout.iterator();
    while( it.hasNext() )
        nextNode = it.next();
        if( nextNode.getDepth() == targetDepth )
            nodeList{n} = nextNode;
            break;
        end
        
    end
end

nodeList

%% compare 'correct' with fast patch similarities

numInis = 10;
numDict = d23c.numDict;
numLoc  = 9;

d23List = cell( numInis, 1 );

for n = 1:numInis
    n
    thisobj = d23c.copy();
    
    [xyzn, iniNode] = thisobj.initializeBuild( );
    
    for d = 1:numLoc-1
        d
        nextNode = thisobj.pickNextNode( d );
        nextNode
        thisobj.build3dPatchIter( xyzn, nextNode, d );
    end
    
    d23List{n} = thisobj;
end


%% visualize results

% d23ListFile = '/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification/matlab/tests/d23TestList.mat';
% d23ListFile = '/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification/matlab/tests/d23TestList_small.mat';
d23ListFile = '/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification/matlab/tests/d23TestList_4.mat';


if( exist( d23ListFile, 'file' ))
    load( d23ListFile );
    thisobj = Dict2dTo3dConsistencyTester([],[],[]);
    thisobj.clone( d23List{1} );
    thisobj.fillObj = d23List{1}.fillObj;
    thisobj.fillObj2 = d23List{1}.fillObj2;
else
    save(d23ListFile,'d23List');
end

%%
figdir = '/groups/saalfeld/home/bogovicj/projects/dictionary/dict2dTo3d/patchSim';

for n = 1:numInis
    
    thisobj = d23List{n};
    
    leafNode1  = thisobj.pickNextNode( numLoc-1 );
    costsByDepth1 = Dict2dTo3dConsistencyTester.getCostsByDepth( leafNode1, numInis, numLoc );
    costsByDepth1
    
    %
    leafNode2  = thisobj.pickCorrespondingNode( leafNode1 );
    costsByDepth2 = Dict2dTo3dConsistencyTester.getCostsByDepth( leafNode2, numInis, numLoc );
    costsByDepth2
    
    %
    cmap = diverging_map( linspace(0,1,8), [0.2 0.2 1], [1 0.1 0.1]);
    for i = 1:8
        
        plot( costsByDepth1{i}, costsByDepth2{i}, 'o', 'MarkerEdgeColor', cmap(i,:));
        hold on;
    end
    set(gca,'color',[0.2 0.2 0.2])
    xlabel( 'patch similarity (correct)');
    ylabel( 'patch similarity (fast)');
    
    % plot identity line
    axbnds = axis;
    axisX = axbnds(1:2);
    axisY = axbnds(3:4);
    
    minBnd = min([axisX(1) axisY(1)]);
    maxBnd = min([axisX(2) axisY(2)]);
    
    hold on;
    plot( [minBnd maxBnd], [minBnd maxBnd], '--w' );
    pause;
    
    % export_fig('similarityComparison_small.png', '-m2');
    export_fig(sprintf('%s/similarityComparison_mid_%d.png', figdir, n), '-m2');
    close all;
    
end

%%


%%

load('/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification/matlab/tests/d23TestList_4.mat');

numLoc = 9
n = 9

%%

wrong = 0;
right = 0;
N = 10;
M = size(thisobj.dimXyzList,1);

rightRankList = zeros( N*M, 1 );
nn = 1;

testNewer = 1;
testBest  = 0;

for n = 1:N
    
    if( testNewer )
        thisobj = d23List{n};
        leafNode  = thisobj.pickNextNode( numLoc-1 );
    else
        if( testBest )
            leafNode = bestList{n};
        else
            leafNode = nodeList{n};
        end
    end
    
    pathToRoot = leafNode.pathToRoot();
    [ pv, patch, cmtx, b ] = patchFromParams( thisobj, leafNode );
    
%     figure; imdisp3d( patch );
    
    % see if the 2d-LR projections correspond to patches
    patchIdxVals = thisobj.leafNodeToPatchIdxs( leafNode );
    

    for i = 1:M
        
        dim = thisobj.dimXyzList( i, 1 );
        xyz = thisobj.dimXyzList( i, 2 );
        
        msk = Dict2dTo3d.planeMaskF( thisobj.sz3d, xyz, dim, thisobj.f );
        projection = mean( ...
            reshape(patch( (msk>0) ), thisobj.sz3d ./ [1 1 thisobj.f]),...
            3);
        
        j = patchIdxVals( i, 1 );
        selectedPatch = reshape(thisobj.D2d( j, : ), thisobj.sz2d);
        diff = norm( selectedPatch(:) - projection(:) );
        
        distancesToDict = pdist2( projection(:)', thisobj.D2d );
        
        [~,k] = min(distancesToDict);
        [ ~, sortedCoords ] = sort( distancesToDict );
        
        if( k ~= j )
            wrong = wrong + 1;
        else
            right = right + 1;
        end
        
        rightRankList(nn) = find( sortedCoords == j );
        nn = nn + 1;
    end
end

fprintf('matched %d / %d\n', right, (right+wrong));
fprintf('in top 5 %d / %d\n', nnz(rightRankList <= 5), (right+wrong));