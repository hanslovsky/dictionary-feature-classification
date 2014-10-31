classdef PatchConstraintsTests < matlab.unittest.TestCase
   
    properties
        pc9;
        pc15;
    end
    
    methods(TestMethodSetup)
        function createInstance(tc)
            tc.pc9  = PatchConstraints(  9, 3 );
            tc.pc15 = PatchConstraints( 15, 5 );
        end
    end
	methods(TestMethodSetup)
        function destroyInstance(tc)
            tc.pc9  = [];
            tc.pc15 = [];
        end
    end
    
    methods (Test)
        
        function testCmtxConstraintsSingle( tc )
            import net.imglib2.algorithms.opt.astar.*;
            import net.imglib2.algorithms.patch.*;
            
            objlist = { tc.pc9, tc.pc15 };
            
            for i = 1:length(objlist)
                
                thisobj = objlist{i};
                thisobj.buildCmtx();
                
                dxl = thisobj.dimXyzList;
                for j = 1:size(dxl,1)
                    stn = SortedTreeNode( ...
                            SubPatch2dLocation( ...
                                dxl(j,1), dxl(j,2), 4, 3 ...
                             ));
                         
                    cmtx = thisobj.cmtxFromConstraints( stn );
                    
                    % make sure size is correct
                    tc.assertSize( cmtx, [ prod(thisobj.sz2d), prod(thisobj.sz3d) ]);
                    
                    % make sure contents is correct
                    msk = Dict2dTo3d.planeMaskF( thisobj.sz3d, dxl(j,2), dxl(j,1), thisobj.f );
                    cmtxTrue = Dict2dTo3d.contraintsMtx( true(thisobj.sz3d), msk );
                    
                    % figure; imagesc( cmtx );
                    % figure; imagesc( cmtxTrue );
                    
                    tc.assertEqual( cmtx, cmtxTrue );
                end
                
            end  
        end
        
        function testCmtxConstraintsSubset( tc )
            
        end
        
    end
end
