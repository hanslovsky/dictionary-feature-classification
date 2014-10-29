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

D = D(:,1:200);

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


%% compare 'correct' with fast patch similarities

numInis = 1;
numDict = d23c.numDict;
numLoc  = 9;

d23List = cell( numInis, 1 );

for n = 1:numInis
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
d23ListFile = '/groups/saalfeld/home/bogovicj/dev/main/dictionary-feature-classification/matlab/tests/d23TestList_mid.mat';


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

leafNode1  = thisobj.pickNextNode( numLoc-1 );
costsByDepth1 = Dict2dTo3dConsistencyTester.getCostsByDepth( leafNode1, numInis, numLoc );
costsByDepth1

%
leafNode2  = thisobj.pickCorrespondingNode( leafNode1 );
costsByDepth2 = Dict2dTo3dConsistencyTester.getCostsByDepth( leafNode2, numInis, numLoc );
costsByDepth2

%%
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

% export_fig('similarityComparison_small.png', '-m2');
export_fig('similarityComparison_mid.png', '-m2');

%%
% for n = 1:numInis
%     thisobj = d23List{n};
%     
%     filler1 = thisobj.fillObj;
%     filler2 = thisobj.fillObj2;
%     
%     someNode  = thisobj.pickNextNode( numLoc-1 );
%     someNode2 = thisobj.pickCorrespondingNode( someNode );
%     
%     pathToRoot  = someNode.pathToRoot();
%     pathToRoot2 = someNode2.pathToRoot();
%     
%     for d = 1:numLoc
%         baseNode = pathToRoot.get( d-1 );
%         children = baseNode.getChildren();
%         
%         baseNode2 = pathToRoot2.get( d-1 );
%         children2 = baseNode2.getChildren();
%         
%         
%     end
% end
