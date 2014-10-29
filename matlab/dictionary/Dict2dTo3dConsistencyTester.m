classdef Dict2dTo3dConsistencyTester < Dict2dTo3d
    properties
        
        NbestTest = 1000;
        
        fillObj;
        fillObj2;
        
        iteration;
        
    end
    
    methods
        
        function this = Dict2dTo3dConsistencyTester( D2d, sz, f )
            this = this@Dict2dTo3d( D2d, sz, f );
            this.NbestTest = min( this.NbestTest, this.numDict );
        end
        
        function obj = copy(this)
            obj = Dict2dTo3dConsistencyTester( this.D2d, this.sz2d(1), this.f );
            obj.numDict = size( obj.D2d, 1 );
            obj.sz3d = this.sz3d;
            
            obj.D3d = this.D3d;
            obj.numDict3d = this.numDict3d;
            
            obj.allSims = this.allSims;
            obj.Nbest = this.Nbest;
            obj.maxItersPerPatch = this.maxItersPerPatch;
            obj.verbose = this.verbose;
            
            obj.dimXyzList = this.dimXyzList;
            
            % specific to this class
            obj.NbestTest = this.NbestTest;
            obj.fillObj   = this.fillObj;
            obj.iteration = this.iteration;
            
        end
       
        function [xyzn, iniNode] = initializeBuild( this )
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            xyzn = this.randomFillOrder();
            
            % fprintf('initializing with random 2d patch\n');
            ii = xyzn( 1 );
            dimIni = this.dimXyzList( ii, 1 );
            xyzIni = this.dimXyzList( ii, 2 );
            
            
            initialPatchIdx = randi( this.numDict );
            rootSpl = SubPatch2dLocation( dimIni, xyzIni, initialPatchIdx, 0 );
            iniNode = SortedTreeNode( rootSpl );
            iniNode2 = SortedTreeNode( rootSpl ); % another copy
            
            this.fillObj  =  Patch2dFill3d( iniNode );
            this.fillObj2 =  Patch2dFill3d( iniNode2 );
        end
        
%        function initializeBuild( this )
%             import net.imglib2.algorithms.patch.*;
%             import net.imglib2.algorithms.opt.astar.*;
%             
%             xyzn = this.randomFillOrder();
%             build3dPatchSetup@Dict2dTo3d( this );
%             
%             rootNode = this.build3dPatchSetup( xyzn );
%             this.fillObj = Patch2dFill3d( rootNode );
%             this.iteration = 1;
%        end

%        function [node] = firstNode( this, xyzn )
%             dimThis = this.dimXyzList( xyzn(1), 1 );
%             xyzThis = this.dimXyzList( xyzn(1), 2 );
%             idxThis = randi( this.numDict, 1, 1 );
%             
%             node = SortedTreeNode(  ...
%                         SubPatch2dLocation( ...
%                             dimThis, xyzThis, idxThis, 0 ), ...
%                         [] );
%        end
       
        function nextNode = pickNextNode( this, depth )
            fprintf('  picking the next node\n');
            it = this.fillObj.getSet().iterator();
            while( it.hasNext() )
                nextNode = it.next();
                if( nextNode.getDepth() == depth )
                    return;
                end
            end   
        end

        function nextNode = pickCorrespondingNode( this, node )
            fprintf('  picking the next node\n');
            it = this.fillObj2.getSet().iterator();
            dat = node.getData();
            while( it.hasNext() )
                nextNode = it.next();
                nextDat  = nextNode.getData();
                if( dat.dim == nextDat.dim && ...
                    dat.xyz == nextDat.xyz && ...
                    dat.idx== nextDat.idx ...
                    )
                        return;
                end
            end
        end
        
        function [patchParams, iteration] = build3dPatchIter( this, xyzn, prevNode, iterin )
            
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            N = size( this.dimXyzList, 1 );
            
            depth = prevNode.getDepth();
            if( depth == N-1 )
                % the next node is best and describes the % parameters for the locally optimal patch
                patchParams = prevNode;
            end
            
            prevNode2 = this.pickCorrespondingNode( prevNode );
            
            iiThis = xyzn( depth + 2 ); % +1 to go from zero to one indexing
                                        % +1 to get the NEXT value
            dimThis = this.dimXyzList( iiThis, 1 );
            xyzThis = this.dimXyzList( iiThis, 2 );
            
            tmpNode = SortedTreeNode(  ...
                SubPatch2dLocation( dimThis, xyzThis, -1, -1 ), ...
                prevNode );
            
           	tmpNode2 = SortedTreeNode(  ...
                SubPatch2dLocation( dimThis, xyzThis, -1, -1 ), ...
                prevNode2 );
            
            % compute all costs
            costs = this.computeCandidateCosts( xyzThis, dimThis, prevNode );
            [ sortedCosts, sortedCostIdxs ] = sort( costs );
            
            % add top NbestTest to a list of candidates
            candList = java.util.ArrayList( this.NbestTest );
            for nn = 1:this.NbestTest
                val = sortedCosts(nn);
                spl = SubPatch2dLocation( dimThis, xyzThis, sortedCostIdxs(nn), val );
                candList.add( spl );
            end
            
            % costs 2
            [ xsectList ] = Dict2dTo3d.intersectingParents( tmpNode2 );
            if( isempty( xsectList ))
                
            else
                costs2 = this.patchCosts( iiThis, xsectList );
                [ sortedCosts2, sortedCostIdxs2 ] = sort( costs2 );
            end
            
            candList2 = java.util.ArrayList( this.NbestTest );
            for nn = 1:this.NbestTest
                val2 = sortedCosts2(nn);
                spl2 = SubPatch2dLocation( dimThis, xyzThis, sortedCostIdxs2(nn), val2 );
                candList2.add( spl2 );
            end
            
            
            prevNode.removeChild( tmpNode );
            prevNode2.removeChild( tmpNode2 );
            
            this.fillObj.addFromNode(   prevNode, candList );
            this.fillObj2.addFromNode( prevNode2, candList2 );
            
            iteration = iterin + 1;
            
        end
       
    end
    
    methods( Static )
        
        function costsByDepth = getCostsByDepth( leafNode, numInis, numLoc )
            costsByDepth = cell( numInis, 1 );
            
            for n = 1:numInis
                
                pathToRoot  = leafNode.pathToRoot();
                
                for d = 1:numLoc
                    baseNode = pathToRoot.get( d-1 );
                    children = baseNode.getChildren();
                    vals = nan( children.size(), 1);
                    
                    it = children.iterator();
                    while( it.hasNext() )
                        dat= it.next().getData();
                        if(dat.idx <= 0 )
                            continue;
                        end
                        vals( dat.idx ) = dat.val;
                    end
                    vals = vals( ~isnan(vals));
                    costsByDepth{d} =  vals;
                    
                end
            end
            
        end
    end
end