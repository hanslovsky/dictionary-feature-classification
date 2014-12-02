classdef Dict2dTo3d < handle
    % Dict2dTo3d 
    %
    % This class contains logic that builds a high-res 3d dictionary from 
    % 2d patches at a lower resolution.
    %
    % John Bogovic
    % HHMI
    % September 2014
    
    properties ( SetAccess = protected )

        obj_fn;   % the file where this object is stored
        
        pc;       % PatchConstraints object
        
        D2d;      % The dictionary elements
        numDict;  % Number of 2d dictionary elements
        sums2d;   % the dictionary element sums
        
        D3d;        % the 3d dictionary
        numDict3d;  % number of 3d dictionary elements
        paramList;  % the parameters that gave rise to the 3d dict
        costList;   % the costs of each patch
        
        sz2d;     % the size of 2d patches
        sz3d;     % the size of 3d patches
        
        f;        % the downsampling factor ( resolution in z of patches )
        
        summer2Dxy; % downsampling by sum function
        
        % Constraint matrices
        constraintMtxX;
        constraintMtxY;
        constraintMtxZ;
        
        % 
        p2dFill3d;
        pairLocRng;
        dimXyzList;
        
        allSims;
        
        scaleByOverlap = 0;
        
        ndims = 3; 
    end
    
    properties 
        % optimization options
        Nbest = 5;
        maxItersPerPatch = 5000;
        minDictElemDiff = 0.001;
        maxDictCost     = -1;
        
        
        
        % patch configuration options
        overlappingPatches = 1;
        saveParams = 1;
        saveCosts = 1;
        
        % initialization options
        Dini;     % dictionary or observations optionally used for initialization
        iniLocs;  % patch locations that are initialized
        numIniPatches; 
        
        verbose = 1;
    end
    
    methods

        function this = Dict2dTo3d( D2d, sz, f, overlappingPatches, scaleByOverlap )
        % Constructor
        % D2d - 
        % sz  - size of 2d patches
        % f   - downsampling factor  
            
%             import net.imglib2.algorithms.opt.astar.AStarMax;
            import net.imglib2.algorithms.patch.*;
            import java.util.*;

            this.D2d     = D2d;
            if( isempty( this.D2d ))
                return;
            end
            this.numDict = size( this.D2d, 1 );
            this.numDict3d = this.numDict;
           
            this.f    = f;
            
            if( length( sz ) > 2 )
               error('sz must be a 2-vector or scalar');
            end
            if( length( sz ) == 2 )
                if( sz(1) ~= sz(2) )
                    error('sz must be a 2-vector');
                end
                this.sz2d = sz;
                this.sz3d = [ sz sz(1) ];
            end
            if( isscalar( sz ))
                this.sz2d = [ sz sz ];
                this.sz3d = [ sz sz sz ];
            end
                        
%             this.sums2d = Dict2dTo3d.allSums( this.D2d, this.sz2d, this.f);
%             this.summer2Dxy = Tid.sum2dxy();

            half = (this.f - 1)./2;
%             this.pairLocRng = (1+half) : this.f : this.sz3d(1);
            this.pairLocRng = (1) : this.f : this.sz3d(1);
            
            % TODO - remove these eventually since the functionallity
            % is now the in the pc object
%             [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
%             this.dimXyzList = [ dList(:), xyzList(:) ];
            
            if( exist('overlappingPatches','var'))
                this.overlappingPatches = overlappingPatches;
            end
            if (exist('scaleByOverlap','var') && ~isempty(scaleByOverlap))
                this.scaleByOverlap = scaleByOverlap;
            end
            
            % remove low variance (ie uninteresting) dictionary
%             this.cleanDict2d();
            
            % use a PatchConstraints object to compute
            % the constraint matrix once
            this.iniPatchConstraints();

            %this.buildIniLocs();
        end
        
        function setScaleByOverlap( this, scaleByOl )
            this.scaleByOverlap = scaleByOl;
            this.iniPatchConstraints();
        end
        
        function setD3d( this, dict3d )
            this.D3d = dict3d;
        end
       
        function obj = copy(this)
            obj = Dict2dTo3d( this.D2d, this.sz2d(1), this.f );
            obj.sz3d = this.sz3d;
            
            obj.D3d = this.D3d;
            obj.numDict3d = this.numDict3d;
            
            obj.clone( this );
        end
        
        function clone(this, other)
            %
            this.f = other.f;
            
            this.D2d  = other.D2d;
            this.numDict = size( other.D2d, 1 );
            this.sz2d = other.sz2d;
            
            this.sz3d = other.sz3d;
            this.D3d  = other.D3d;
            this.numDict3d = other.numDict3d;
            
            this.paramList = other.paramList;
            this.costList  = other.costList;
            
            this.allSims = other.allSims;
            this.Nbest = other.Nbest;
            this.maxItersPerPatch = other.maxItersPerPatch;
            this.verbose = other.verbose;
            
            this.dimXyzList = other.dimXyzList;
        end
        
        function idx = locXyzDim2Idx( this, dim, xyz )
            idx = ( this.dimXyzList(:,1) == dim & ...
                    this.dimXyzList(:,2) == xyz );
        end
        
        function [dim,xyz] = locIdx2XyzDim( this, idx )
            dim = this.dimXyzList( idx, 1 );
            xyz = this.dimXyzList( idx, 2 );
        end
        
        function obj_fn = save( this, fn )
            global DFEVAL_DIR;
            if( exist( 'fn', 'var') && ~isempty( fn ) )
                this.obj_fn = fn;
            elseif( isempty( this.obj_fn ))
                this.obj_fn = fullfile( DFEVAL_DIR, sprintf('Dict2dTo3d_%s.mat', ...
                    datestr( now, 30 ) ));
            end
            
            obj_fn = this.obj_fn;
            if( ~exist( obj_fn, 'file'))
                save( this.obj_fn, 'this' ); 
            end
        end
        
        function numRemoved = cleanDict2d( this, minVar )
            if( this.verbose)
               fprintf('Cleaning 2d dict\n'); 
            end
            
            if( ~exist( 'minVar','var') || isempty( minVar ))
                minVar = 0.005;
            end
            vlist = zeros( this.numDict, 1);
            for i = 1:this.numDict
                vlist(i) = var( this.D2d(i,:));
            end
            
            goodDictElems = ( vlist > minVar );
            this.D2d = this.D2d( goodDictElems, : );
            
            newSize = nnz( goodDictElems );
            numRemoved = this.numDict - newSize;
            this.numDict = newSize;
        end
        
        function trimDict2d( this )
        % removes "extra" elements from the 2d dictionary
        % that might arise from initializations
            this.D2d = this.D2d( 1:this.numDict, : );
        end
        
        function patch = sample2dTo3d( this )
            import net.imglib2.algorithms.patch.*;
            
            szSm = this.sz2d ./ this.f;
            
            % stores the order in which patches are added to particular bins 
            % (hard-coded '3' corresponds to 3 spatial-dimensions)
            param = zeros( [ 3 max(szSm) ] ); 
            
            % stores which observed 2d patch goes into which 3d position / orientation "bin"
            patchIdxs = zeros( [ 3 max(szSm) ] );
            
            % patch Loc and Idx
            % [ dim xyz idx ]
            patchLocIdx = zeros( 3.*max(szSm), 3 );
            
            % how many bins are there
            mskN = numel( param );
            
            % randomize the order in which the bins will be filled
            xyzn = randperm( mskN );
            
            for ii = 1:mskN
                fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
                
                dim = this.dimXyzList(ii,1);
                xyz = this.dimXyzList(ii,2);
                dxyzi1 = xyzn(ii);
                
                if( ii > 1 )
                    
%                     xsectList = Dict2dTo3d.findIntersections( param, [dim xyz]);
%                     xsectList

                    [ xsectCoords, ~ ] = Dict2dTo3d.findIntersectionsList( patchLocIdx, dim,  xyz);
                    % xsectCoords
                    
                    % get the next best candidate for patch 
                    
                    
                    for nn = 1:size( xsectCoords, 1 )
                        
                        fixedDim = xsectCoords(nn,1);
                        fixedXyz = xsectCoords(nn,2);
                        fixedIdx = xsectCoords(nn,3);
                        
                        dxyzi2 = ( this.dimXyzList(:,1) == fixedDim & ...
                                   this.dimXyzList(:,2) == fixedXyz);
                                
%                         dxyzi2 = find( this.dimXyzList(:,1) == fixedDim & ...
%                                        this.dimXyzList(:,2) == fixedXyz);
%                         tsim = this.allSims(  fixedIdx, :, dxyzi1, dxyzi2 );

                        % get similarities for all patches and the given orientations
                        tsim = this.allSims(  fixedIdx, :, dxyzi1, dxyzi2 )

                        [tsimSort, sortIdx ] = sort( tsim );
                      
                        for mm = 1:this.Nbest
                            aa=1;
                        end
                        
                    end
                    
                else
                    rootSpl = SubPatch2dLocation( dim, xyz, 0 );
                    this.astar.setRoot( rootSpl );
                    patchIdxs(ii) = randi( this.numDict, 1 ); 
                    patchLocIdx(ii,:) = [ dim, xyz, patchIdxs(ii) ];
                end
                
                param( xyzn(ii) ) = ii;
                patchLocIdx(ii,:) = [ dim, xyz, randi( this.numDict, 1 ) ];
                
            end
            
            patch = [];
        end
    
        function [ patches3d, f_out ] = build3dDictionaryDist( this, dict3dSize, pe_batch_size, inioutmethod )
            global DICTPATH;

            f_out = [];
            if( ~exist( 'pe_batch_size', 'var' ) || isempty( pe_batch_size ))
                pe_batch_size = 1;
            end
            
            if( ~exist( 'inioutmethod', 'var' ) || isempty( inioutmethod ))
                inioutmethod = '';
            end
            
%             if( strcmp(inioutmethod, 'dict') )
%                 buildName = 'build3dPatchIni';
%             else
                buildName = 'build3dPatch';
%             end
            
            % params
            this.numDict3d = dict3dSize;
            fun = @run_obj_method_dist;
            use_gpu = 0;
            run_script = fullfile( DICTPATH, 'bin', 'my_run_runObjMethod.sh');
            
            n = 0;
            patches3d = zeros( dict3dSize, prod( this.sz3d ));
            
            % the cell array of parameters
            num_out_args = 6; % build3dPatch has 4 output args
            
%             varargin = {  repmat( {this.obj_fn}, dict3dSize-n, 1), ...
%                           {'this'}, {'build3dPatch'}, {num_out_args} };
            this.save();
            
            varargin = {  repmat( {this.obj_fn}, dict3dSize-n, 1), ...
                          repmat( {'this'}, dict3dSize-n, 1), ...
                          repmat( {buildName}, dict3dSize-n, 1), ...
                          repmat( {num_out_args}, dict3dSize-n, 1), ...
                       };
                      
            if( inioutmethod )
                iniParams = this.iniParamsDist( dict3dSize, inioutmethod );
                varargin{5} = iniParams;
                
                if( strcmp(inioutmethod, 'dict') )
                    varargin{6} = {this.iniLocs};
                end
            end
            
            % optionally save the parameters along with the final recon patch
            if( this.saveParams )
                this.paramList = cell( dict3dSize, 1);
            end
            if( this.saveCosts )
                this.paramList = cell( dict3dSize, 1);
            end
            
            iter = 1;
            while( n < dict3dSize )
                if( this.verbose )
                    fprintf('iter %d, n = %d\n', iter, n );
                    fprintf('submitting %d jobs\n', (dict3dSize-n) );
                end
                varargin{1} = repmat( {this.obj_fn}, dict3dSize-n, 1);
                
                % run everything
                f_out = qsub_dist(fun, pe_batch_size, use_gpu, ...
                            [], [], run_script, ...
                            varargin{:} );
                
                goodToGo = true( size(f_out,1), 1 );
                for i = 1:size(f_out,1)
                    
                    patchParams = f_out{i}{1}{1};
                    pv          = f_out{i}{1}{2};
                    
                    % pv = this.patchFromParams( patchParams );
                    pv = vecToRowCol( pv, 'row' );
                    
                    if( ~isempty( pv ))
                        if( n > 0 && this.minDictElemDiff > 0 )
                            mindiff = min(pdist2( pv, patches3d(1:n,:)));
                            if( mindiff < this.minDictElemDiff )
                                goodToGo(i) = false;
                                continue;
                            end
                        end
                        
                        n = n + 1;
                        if( this.saveParams )
                            this.paramList{n} = patchParams;
                        end
                        if( this.saveCosts )
                            this.costList{n} = f_out{i}{1}{4};
                        end
                        
                        patches3d(n,:) = pv;
                        
                    end
                end
                
                this.saveSupplemental( f_out(goodToGo) );
                
                iter = iter + 1;
            end
            
            this.D3d = patches3d;
        end
        
        function saveSupplemental( this, f_out )
            % do nothing, 
            % but subclasses may use.
        end
        
        function [ patches3d ] = build3dDictionary( this, dict3dSize, inioutmethod )
            this.numDict3d = dict3dSize;
            
            haveinis = 0;
            if( exist( 'inioutmethod', 'var' ) && ~isempty( inioutmethod ) )
                iniParams = this.iniParamsDist( dict3dSize, inioutmethod );
                haveinis = 1;
            end
            
            n = 0;
            i = 1;
            patches3d = zeros( dict3dSize, prod( this.sz3d ));
            while( n < this.numDict3d )
                if( this.verbose )
                    fprintf('building 3d dictionary element %d\n', n);
                end
                if( haveinis )
                    patchParams = this.build3dPatch( iniParams{i} );
                    i = i + 1;
                else
                    patchParams = this.build3dPatch();
                end
                
                pv = this.patchFromParams( patchParams );
                pv = vecToRowCol( pv, 'row' );
                
                if( ~isempty( pv ) )
                    if( n > 0 && this.minDictElemDiff > 0 )
                        mindiff = min(pdist2( pv, patches3d(1:n,:)));
                        if( mindiff < this.minDictElemDiff )
                           continue; 
                        end
                    end
                    
                    n = n + 1;
                    patches3d(n,:) = pv;    
                    
                end
            end
            
            this.D3d = patches3d;
        end
        
        function [ constrList, mskList ] = renderConstraints( this, paramsIn )
            if( isa( paramsIn, 'net.imglib2.algorithms.opt.astar.SortedTreeNode'))
                params = patchParamNodeToArray( paramNode );
            else 
                params = paramsIn;
            end
            
            if( size(params,2) == 3 )
                idxs = params(:,3);
                dims = params(:,1);
                xyzs = params(:,2);
            else
                idxs = params;
                dims = this.pc.dimXyzList(:,1);
                xyzs = this.pc.dimXyzList(:,2);
            end
            
            numConstraints = size( params, 1 );
            
            domask = 0;
            mskList = [];
            if( nargout > 1 )
                domask = 1;
                mskList = cell( numConstraints, 1);
            end
            
            constrList     = cell( numConstraints, 1);
            for n = 1:numConstraints 
                 dim = dims(n);
                 xyz = xyzs(n);
                 msk = Dict2dTo3d.planeMaskF( this.sz3d, xyz, dim, this.f, 0 );
                 
                 thisConstraint = zeros( size(msk) );
                 dictElem = this.D2d( idxs(n), : );
                 mskvals = ( msk > 0 );
                 thisConstraint( mskvals ) = dictElem( msk( mskvals ));
                 constrList{n} = thisConstraint;
                 
                 if( domask )
                    mskList{ n } = msk;
                 end
            end
        end
        
        function [ pvHR, patchHR ] = reconPatch( this, patchLR, doNorm )
            
            % doNormalization by default
            if( ~exist('doNorm','var') || isempty(doNorm))
                doNorm = 1;
            end
            
            patchDim = length( size(patchLR));
            if( patchDim == 3 )
                patch = zeros( size(patchLR,3), numel( patchLR(:,:,1)));
                for z = 1:size(patchLR,3)
                    patch(z,:) = reshape( patchLR(:,:,z), 1, [] );
                end
                % normalize the observed patches
                if( doNorm )
                    patch = patch ./ repmat( sqrt(sum(patch.^2, 2)), 1, size(patch,2));
                end
            else
                patch = patchLR;
            end
            
            numAdded = this.addTemporaryPatchesToDict( patch );
            iniIdx   = (this.numDict-numAdded+1) : (this.numDict);
            
            paramList = [ 3 3 3; 1 4 7; iniIdx]';
            
            [rc, lc] = rootConstraintFromArray( paramList );
            
            this.addToAllSims( numAdded );
            
            patchParams = this.build3dPatch( lc, numAdded );
            
            [pvHR, patchHR] = this.patchFromParams( patchParams );
            
            % fprintf('finished\n');
            
            this.removeTemporaryPatches( numAdded );
            this.remFromAllSims( numAdded );
            
        end

        function [ pv, patch, cmtx, b ] = patchFromParams( this, nodeOrList )
            if( ~isempty(this.pc))
                [ pv, patch, cmtx, b ] = patchFromParamsPre( this, nodeOrList);
            else
                [ pv, patch, cmtx, b ] = patchFromParamsComp( this, nodeOrList);
            end
        end
        
        function [ pv, patch, cmtxi, b ] = patchFromParamsPre( this, splNodeIn )
            % Compute patch from parameters using a precomputed constraint
            % matrix from the PatchConstraints object
            
            fprintf('patch from precomputed constraints\n');
            if( isempty( splNodeIn ))
                patch = [];
                pv = [];
                return;
            end
            
            docell = 0;
            if( ~iscell( splNodeIn ))
                splNodeList = { splNodeIn };
            else
                docell = 1;
                splNodeList = splNodeIn;
                pv    = cell( length(splNodeList), 1 );
                patch = cell( length(splNodeList), 1 );
            end
            
            cmtxi = this.pc.cmtxInv;
            
            for n = 1:length(splNodeList)
                splNode = splNodeList{ n };
                b = this.constraintValue( splNode );
                
                if( docell )
                    pv{n} = cmtxi * b;
                    if( nargout > 1 )
                        patch{n} = reshape( pv{n}, this.sz3d );
                    end
                else
                    pv = cmtxi * b;
                    if( nargout > 1 )
                        patch = reshape( pv, this.sz3d );
                    end
                end
                
            end
        end

        function [ pv, patch, cmtx, b ] = patchFromParamsComp( this, splNodeIn )
        % splNode must be a sortedTreeNode< SubPatch2dLocation >   
            if( isempty( splNodeIn ))
                patch = [];
                pv = [];
                return;
                
            end
            docell = 0;
            if( ~iscell( splNodeIn ))
                splNodeList = { splNodeIn };
            else
                splNodeList = splNodeIn;
                pv = cell( length(splNodeList) );
                docell = 1;
            end
            
            for n = 1:length(splNodeList)
                splNode = splNodeList{ n };
                
                
                % check validity of input
                if( isa( splNode, 'net.imglib2.algorithms.opt.astar.SortedTreeNode'))
                    if( ~ isa( splNode.getData(), 'net.imglib2.algorithms.patch.SubPatch2dLocation' ))
                        error( 'splNode must contain a ''net.imglib2.algorithms.patch.SubPatch2dLocation'' ');
                    end
                else
                    % here
                    try
                        splNode = Dict2dTo3d.patchParamArrayToNode(splNode);
                    catch e
                        error('tried to make a node from an array but failed');
                    end
                    %                 error( 'splNode must be a ''net.imglib2.algorithms.opt.astar.SortedTreeNode'' ');
                end
                
                paramNode = splNode;
                param     = paramNode.getData();
                
                numConstraints = splNode.getDepth() + 1;
                
                N = numConstraints .* prod( this.sz2d );  % the number of constraints
                M = prod( this.sz3d ); % the number of elements in the HR patch
                
                cmtx = zeros( N, M );
                b    = zeros( N, 1 );
                
                k = 1;
                for i = 1:numConstraints
                    
                    dim = param.dim;
                    xyz = param.xyz;
                    idx = param.idx;
                    %pause;
                    
                    msk = Dict2dTo3d.planeMaskF( this.sz3d, xyz, dim, this.f);
                    
                    for j = 1:prod(this.sz2d)
                        cmtx( k, (msk == j) ) = 1;
                        b( k ) = this.D2d( idx, j );
                        k = k + 1;
                    end
                    
                    if( ~paramNode.isRoot())
                        paramNode = paramNode.getParent();
                        param     = paramNode.getData();
                    end
                end
                
                if( docell )
                    pv = pinv( cmtx ) * b;
                     if( nargout > 1 )
                        patch = reshape( pv, this.sz3d );
                    end
                else
                    pv{n} = pinv( cmtx ) * b;
                    if( nargout > 1 )
                        patch{n} = reshape( pv{n}, this.sz3d );
                    end
                end
                
            end % loop over nodes
        end
       
        function xyzn = randomFillOrderBlah( this )
            N = size( this.dimXyzList, 1 );

            isgood = 0;
            xyzn = [];
            while( ~isgood )
                xyzn = randperm( N );
                isgood = this.validateFillOrder( xyzn );
            end
        end

        function xyzn = randomFillOrder( this )
            N = size( this.dimXyzList, 1 );
            i1 = randi( N, 1, 1 );

            % pick a second plane that intersects the first
            validOpts = find(this.dimXyzList(:,1) ~= this.dimXyzList(i1,1));
            j = randi( length(validOpts), 1, 1);
           
            i2 = validOpts(j);

            % randomly choose the order of the remaining 
            % plane orientations
            remi = setdiff( 1:N, [i1 i2] );
            xyzn = [i1 i2 remi(randperm(N-2))];

        end

        function isgood = validateFillOrder( this, xyzn )
        % xyzn is a fill order - determines if it is 'valid'
        % where valid means that for i > 1, the plane defined by
        % there exists a j such that 
        %   this.dimXyzList(i) intersectts with this.dimXyzList(j) 
        
            %TODO - determine if we need this and implement if so
            %fillOrder = this.dimXyzList( xyzn, : );
            isgood = 1;
        end

        function [patchParams, pv, iteration, costs, paramsByIter, otherParams] = build3dPatch( ...
                        this, iniPatchFill, num2exclude )
            
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            paramsByIter = {};
            otherParams  = {};
            
            if( ~exist('num2exclude','var') || isempty(num2exclude))
                num2exclude = 0;
            end
            
            N = size( this.dimXyzList, 1 );
            
            % randomize the order in which the patch bins will be filled
            % xyzn = randperm( N );
            xyzn = this.randomFillOrder();
            
            haveIni = 0;
            if( exist( 'iniPatchFill', 'var' ) && ~isempty( iniPatchFill ))
                
                % fprintf('initializing with the input\n');
                haveIni = 1;
                
                this.p2dFill3d = Patch2dFill3d( iniPatchFill );
            
                iniParams = Dict2dTo3d.patchParamNodeToArray( iniPatchFill );
                
                exclude = find( ismember(this.dimXyzList(:,1), iniParams(:,1)) & ...
                                ismember(this.dimXyzList(:,2), iniParams(:,2)));
                
                xyzn = [exclude', setdiff( xyzn, exclude, 'stable')];
                
            end
            
            done = false;
            iteration = 1;
            while( ~done )
                if( this.verbose && mod(iteration, 100 ) == 0 )
                      fprintf('iteration %d\n', iteration); 
                end
                if ( iteration > this.maxItersPerPatch )
                    
                    patchParams = this.currentBestPatchParams();
                    return;
                    
                elseif( iteration == 1 && ~haveIni )
                   
                    rootNode = this.build3dPatchSetup( xyzn );
                    this.p2dFill3d = Patch2dFill3d( rootNode );
                    
                else
                    
                    prevNode = this.p2dFill3d.next();
                    depth = prevNode.getDepth();
                    if( depth == N-1 )
                        
                        % the next node is best and describes the % parameters for the locally optimal patch
                        patchParams = prevNode;
                        done = true;
                        break;
                    end
                    
                    % iiPrev = xyzn( depth + 1 );
                    iiThis = xyzn( depth + 2 );
                    
                    dimThis = this.dimXyzList( iiThis, 1 );
                    xyzThis = this.dimXyzList( iiThis, 2 );
                    
                    tmpNode = SortedTreeNode(  ...
                                SubPatch2dLocation( dimThis, xyzThis, -1, -1 ), ...
                                prevNode );
                    [ xsectList ] = Dict2dTo3d.intersectingParents( tmpNode );
                   
                    if( isempty( xsectList ))
                        
%                         fprintf('no intersections...picking randomly\n');
                        % pick 'Nbest' patches randomly and give them equal cost 
                        randomPatchIndexes = randi( this.numDict, 1, this.Nbest);
                        costs = ones( this.numDict, 1 );
                        
                        % TODO - what cost should be given?
                        costs( randomPatchIndexes ) = 0;
                    else
%                         fprintf('computing costs\n');
                        if( isempty( this.allSims ))
                            %costs = this.patchCostsAllSlow( xyzThis, dimThis, xsectList );
                            costs = this.computeCandidateCosts( xyzThis, dimThis, prevNode );
                        else
                            costs = this.patchCosts( iiThis, xsectList );
                        end
                        
                    end
                    prevNode.removeChild( tmpNode );

                    costs = costs( 1: length(costs)-num2exclude);
                    
                    [ sortedCosts, sortedCostIdxs ] = sort( costs );
                    
                    
                    candList = java.util.ArrayList( this.Nbest );
                    for nn = 1:this.Nbest
                       val = sortedCosts(nn);
                       spl = SubPatch2dLocation( dimThis, xyzThis, sortedCostIdxs(nn), val );
                       candList.add( spl );
                    end
                    
                    this.p2dFill3d.addFromNode( prevNode, candList );
                    
                end
                iteration = iteration + 1;
            end
            
            [pv] = this.patchFromParams( patchParams );
        end
         
        function costs = computeCandidateCosts( this, xyz, dim, config )
            fprintf('USING NEW COST METHOD!\n');
            costs = this.allPatchConfigCosts( dim, xyz, config );
        end
        
        function currentBest = currentBestPatchParams( this )
            N = size( this.dimXyzList, 1 );
            it = this.p2dFill3d.getSet().iterator();
            while( it.hasNext())
               thisone = it.next(); 
               if( thisone.getDepth() == (N-1))
                  currentBest = thisone;
                  return;
               end
            end
            currentBest = [];
        end
        
        function constraints = getSumConstraints( this, patchIdxs )
            
            szdown = this.sz3d ./ this.f;
            constraints = nan( szdown );
            
            [d,xyz] = ind2sub( size( patchIdxs ), find( patchIdxs > 0 ));
            for i = 1:length(xyz)
                patch = reshape( this.D2d( patchIdxs(d(i),xyz(i)), : ), this.sz2d );
                constraints( Dict2dTo3d.planeMask( szdown, xyz(i), d(i) )) = this.summer2Dxy( patch, this.f );
            end
        end
        
        function c = collectConstraints( this, constraints, xyz, d )
            szSm3 = this.sz3d ./ this.f;
            szSm2 = this.sz2d ./ this.f;
            planeMask = Dict2dTo3d.planeMask( szSm3, xyz, d );
            c = reshape( constraints( planeMask ), szSm2 );
        end
        
        function buildConstraintMatrix( this )
           
            [x,y,z] = ndgrid( 1:this.f, 1:this.f, 1:this.f );
            N = numel(x);
            
            this.constraintMtxX = zeros( this.f, N );
            this.constraintMtxY = zeros( this.f, N );
            this.constraintMtxZ = zeros( this.f, N );
            
            for i = 1:this.f
                % rows
                this.constraintMtxX(i,:) = reshape((x==i),1,[]);
                
                % colums
                this.constraintMtxY(i,:) = reshape((y==i),1,[]);
                
                % slices
                this.constraintMtxZ(i,:) = reshape((z==i),1,[]);
                
            end
        end

        function [ costs ] = patchCosts( this, dxyzi1, intersectionList )
            
            N = size(intersectionList,1);
            constraintCosts = zeros( N, this.numDict );
            
            for nn = 1:N
                
                fixedDim = intersectionList(nn,1);
                fixedXyz = intersectionList(nn,2);
                fixedIdx = intersectionList(nn,3);
                
                dxyzi2 = ( this.dimXyzList(:,1) == fixedDim & ...
                           this.dimXyzList(:,2) == fixedXyz);
                               
                constraintCosts( nn,:) = this.allSims(  fixedIdx, :, dxyzi1, dxyzi2 );
            end
            
%             fprintf('constraint costs\n');
%             size( constraintCosts )
            
            % take max over constraint costs for final cost
            costs = max( constraintCosts, [], 1 );
            
%             fprintf('new patch costs\n');
%             size(costs)
        end
        
        function costs = allPatchConfigCosts( this, dim, xyz, config )
            
            costs = zeros( this.numDict, 1 );
            for i = 1:this.numDict
                costs( i ) = this.patchConfigCost( i, dim, xyz, config );
            end
        end
        
        function [cost, cmtx, b] = patchConfigCost( this, patchIdx, dim, xyz, config )
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            
            tmpNode = SortedTreeNode( ...
                        SubPatch2dLocation( dim, xyz, patchIdx, 9999 ), ...
                        config );
                    
            
            % build constraint mtx
            patchNumElem   = prod( this.sz2d );
            numConstraints = ( tmpNode.getDepth() + 1 ) * patchNumElem;
            numVariables   = prod( this.sz3d );
            
            cmtx = zeros( numConstraints, numVariables );
            b    = zeros( numConstraints, 1 );

            thisNode = tmpNode;
            
            k = 1;
            for i = 1 : ( tmpNode.getDepth() + 1 )
                thisdim = thisNode.getData().dim;
                thisxyz = thisNode.getData().xyz;
                thisidx = thisNode.getData().idx;
                
                % TODO - dont recompute this msk every time
                msk = Dict2dTo3d.planeMaskF( this.sz3d, thisxyz, thisdim, this.f );
            
                for j = 1:patchNumElem
                    
                    cmtx( k, (msk==j) ) = 1;
                    b(k) = this.D2d( thisidx, j );
                    
                    k = k + 1;
                end
                
                thisNode = thisNode.getParent();
            end
            
            % for debug purposes
%             figure; imdisp( cmtx );
%             aaassddeeww=1;
            
            % solve and compute the cost
            cmtxi = pinv( cmtx );
            x = cmtxi * b;
            cost = norm( cmtx * x - b );
            
            config.removeChild( tmpNode );
        end
        
        function [ cost ] = patchCostsAllSlow( this, xyz1, n1, intersectionList )
            N = size(intersectionList,1);
            similarityList = zeros( this.numDict, N);
            for j = 1:this.numDict
                p1 = reshape( this.D2d( j, :), this.sz2d );
                for i = 1:N
                    p2 = reshape( this.D2d( intersectionList(i,3), :), this.sz2d );
                    similarityList(j,i) = this.patchSimilarity( ...
                        p1, xyz1, n1, ...
                        p2, intersectionList(i,2), intersectionList(i,1));
                    
                end
            end
            cost = max( similarityList, [], 2 );
            
        end
        
        function [ cost ] = patchCostsSlow( this, p1, xyz1, n1, intersectionList )
            N = size(intersectionList,1);
            cost = zeros( N, 1);
            for i = 1:N
                p2 = reshape( this.D2d( intersectionList(i,3), :), this.sz2d );
                cost(i) = this.patchSimilarity( ...
                            p1, xyz1, n1, ...
                            p2, intersectionList(i,2), intersectionList(i,1));
                
            end
            cost = max( similarityList );
        end
       
            
        % TODO make this methods abstract when more possibilities are available
        % this version will use unconstrained linear least squares 
        function [ sim, x, cmtx, b, pm1, pm2, overlap ] = patchSimilarity( this, p1, xyz1, n1, p2, xyz2, n2 )
            
            % get the value of the sum-contraint each patch describes over
            % the overlapping area
            sz3 = this.sz3d;

            pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, n1, this.f );
            pm2 = Dict2dTo3d.planeMaskF( sz3, xyz2, n2, this.f );
            overlap = (pm1 > 0) & (pm2 > 0);
           
            [cmtx1, b1] = Dict2dTo3d.contraintsFromMask( p1, overlap, pm1 );
            [cmtx2, b2] = Dict2dTo3d.contraintsFromMask( p2, overlap, pm2 );

            cmtx = [ cmtx1; cmtx2 ] ;
            b    = [ b1; b2 ];
            x = pinv(cmtx) * b;
            
            sim = sum((cmtx*x - b).^2);
        end
        
        function allSims = allSimilaritiesFast( this )
        % returns a 4d array describing the similarities between all pairs
        % of patches in all positions / orientations
        % 
        % dimensions are : ( patch1Index patch2Index patch1XyzDim patch2XyzDim )
        % the this.dimXyzList property stores the patch orientations.
        %
        % As an example:
        % allSims( i, j, k, l ) stores the similarities for patches i and j
        %   for orientations 
        %   this.dimXyzList( k ) and
        %   this.dimXyzList( l ) 
        %
            numP = size( this.D2d, 1 );
            
            sz3 = this.sz3d;
            [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
            N = numel(dList);
            
            % this.dimXyzList = [ dList(:), xyzList(:) ];
            
            allSims = 9999.*ones( numP, numP, N, N );
            for k = 1:N
                xyz1 = xyzList(k);
                n1   = dList(k);
                
                for l = 1:N
                    
                    % don't bother computing similarities
                    % for identical orientations - they'll
                    % never be used
                    if( l == k )
                        continue;
                    end
                    
                    xyz2 = xyzList(l);
                    n2   = dList(l);
                    
                    % skip when the constraint planes are 'parallel'
                    % because this will never be used.
                    if( n1 == n2 )
                       continue; 
                    end
                    
                    pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, n1, this.f );
                    pm2 = Dict2dTo3d.planeMaskF( sz3, xyz2, n2, this.f );
                    overlap = (pm1 > 0) & (pm2 > 0);
                    
                    [ cmtx1, idxs1 ] = Dict2dTo3d.contraintsMtx( overlap, pm1 );
                    [ cmtx2, idxs2 ] = Dict2dTo3d.contraintsMtx( overlap, pm2 );
                    
                    cmtx = [ cmtx1; cmtx2 ];
                    cmtxi = pinv( cmtx );
                    
                    for i = 1:numP
                        b1 = this.D2d( i, idxs1 );
                        
                        for j = i:numP
                            
                            b2 = this.D2d( j, idxs2 );
                            
                            b    = [ b1'; b2' ];
                            x = cmtxi * b;
                            sim = sum(( cmtx*x - b ).^2);

                            allSims( i, j, k, l ) = sim;
                            allSims( j, i, k, l ) = sim;

                        end
                    end
                end
%                 endtime = toc;
%                 fprintf('time: %fs\n', endtime );
            end
            this.allSims = allSims;
        end
        
        
        function allSims = specSimilaritiesFast( this, iprange, jprange )
            
            asym = 0;
            if( ~exist( 'iprange', 'var') || isempty( iprange ) )
               fprintf('default iprange\n');
               iprange =  1 : size( this.D2d, 1 );
               asym = 1;
            end
            
            if( ~exist( 'jprange', 'var') || isempty( jprange ) )
               fprintf('default jprange\n');
               jprange =  1 : size( this.D2d, 1 );
               asym = 1;
            end
            
            numIp = length( iprange );
            numJp = length( jprange );
            
            sz3 = this.sz3d;
            [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
            N = numel(dList);
            
            % this.dimXyzList = [ dList(:), xyzList(:) ];
            
            allSims = 9999.*ones( numIp, numJp, N, N );
            for k = 1:N
                xyz1 = xyzList(k);
                n1   = dList(k);
                
                for l = 1:N
                    
                    % don't bother computing similarities
                    % for identical orientations - they'll
                    % never be used
                    if( l == k )
                        continue;
                    end
                    
                    xyz2 = xyzList(l);
                    n2   = dList(l);
                    
                    % skip when the constraint planes are 'parallel'
                    % because this will never be used.
                    if( n1 == n2 )
                       continue; 
                    end
                    
                    pm1 = Dict2dTo3d.planeMaskF( sz3, xyz1, n1, this.f );
                    pm2 = Dict2dTo3d.planeMaskF( sz3, xyz2, n2, this.f );
                    overlap = (pm1 > 0) & (pm2 > 0);
                    
                    [ cmtx1, idxs1 ] = Dict2dTo3d.contraintsMtx( overlap, pm1 );
                    [ cmtx2, idxs2 ] = Dict2dTo3d.contraintsMtx( overlap, pm2 );
                    
                    cmtx = [ cmtx1; cmtx2 ];
                    cmtxi = pinv( cmtx );
                    
                    for i = 1:numIp
                        b1 = this.D2d( iprange(i), idxs1 );
                        
                        for j = 1:numJp
                            
                            b2 = this.D2d( jprange(j), idxs2 );
                            
                            b    = [ b1'; b2' ];
                            x = cmtxi * b;
                            sim = sum(( cmtx*x - b ).^2);

                            allSims( i, j, k, l ) = sim;
                            
                        end
                    end
                end
%                 endtime = toc;
%                 fprintf('time: %fs\n', endtime );
            end
            
        end
        
        function allSims = allSimilarities( this )
        % returns a 4d array describing the similarities between all pairs
        % of patches in all positions / orientations
        % 
        % dimensions are : ( patch1Index patch2Index patch1XyzDim patch2XyzDim )
        % the this.dimXyzList property stores the patch orientations.
        %
        % As an example:
        % allSims( i, j, k, l ) stores the similarities for patches i and j
        %   for orientations 
        %   this.dimXyzList( k ) and
        %   this.dimXyzList( l ) 
        %
            [dList, xyzList] = ndgrid( 1:3, this.pairLocRng);
            N = numel(dList);
            
%             this.dimXyzList = [ dList(:), xyzList(:) ];
            
            allSims = 9999.*ones( this.numDict, this.numDict, N, N );
            for i = 1:this.numDict
                fprintf('i: %d of %d\n', i, this.numDict );
                p1 = reshape( this.D2d(i,:), this.sz2d );
                for j = i:this.numDict
                    p2 = reshape( this.D2d(j,:), this.sz2d );
                    for k = 1:N
                        xyz1 = xyzList(k);
                        n1   = dList(k);

                        for l = 1:N

                            % don't bother computing similarities
                            % for identical orientations - they'll
                            % never be used
                            if( l == k )
                                continue;
                            end

                            xyz2 = xyzList(l);
                            n2   = dList(l);
                            sim = this.patchSimilarity( p1, xyz1, n1, p2, xyz2, n2 );
                            
                            allSims( i, j, k, l ) = sim;
                            allSims( j, i, k, l ) = sim;
                        end
                    end
                end
%                 fprintf('time: %fs\n', endtime );
            end
            
        end
        
        function b = constraintValue( this, obj )
            if( isa( obj, 'net.imglib2.algorithms.opt.astar.SortedTreeNode'))
                b = this.pc.constraintValueNode( this.D2d, obj );
            else
                b = this.pc.constraintValueList( this.D2d, obj );
            end
        end
    end
    
    methods ( Access = protected )
        function numAdded = addTemporaryPatchesToDict( this, tempPatches, updateNumDict )
            if( ~exist('updateNumDict','var') || isempty(updateNumDict))
                updateNumDict = 0;
            end
            
            numAdded = size( tempPatches, 1 );
            this.D2d = [ this.D2d; tempPatches ];
            if( updateNumDict )
                this.numDict = size( this.D2d, 1 );
            end
        end
        
        function removeTemporaryPatches( this, numToRemove )
            this.D2d = this.D2d( 1: size(this.D2d,1)-numToRemove, : );
            this.numDict = size( this.D2d, 1 );
        end
        
        function newAllSims = remFromAllSims( this, numRem )

            sz1 = size( this.allSims, 1 );
            sz2 = size( this.allSims, 2 );
            newAllSims = this.allSims( 1:sz1-numRem, 1:sz2-numRem, :, :);
            this.allSims = newAllSims;
        end
        
        function newAllSims = addToAllSims( this, numNew ) 
            
            oldAllSims = this.allSims;
            numOld = size(oldAllSims,1);
            
            irng = 1:numOld ;
            jrng = (numOld +1):(numOld +numNew);
            
            newSims = this.specSimilaritiesFast( jrng, irng );
            
            newSimsSize = size( oldAllSims );
            newSimsSize(1:2) = newSimsSize(1:2) + numNew ;
            
            newAllSims = zeros( newSimsSize );
            
            newAllSims( irng, irng, :, : ) = oldAllSims;
            newAllSims( jrng, irng, :, : ) = newSims;
            newAllSims( irng, jrng, :, : ) = permute( newSims, [2 1 3 4]);
            newAllSims( jrng, jrng, :, : ) = 9999; % make this stupid high
            this.allSims = newAllSims;

        end

        function initialConstraintNode = build3dPatchSetup( this, xyzn )
            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;
            
            % fprintf('initializing with random 2d patch\n');
            ii = xyzn( 1 );
            dimIni = this.dimXyzList( ii, 1 );
            xyzIni = this.dimXyzList( ii, 2 );
            
            % dxyzi1 = xyzn(ii);
            
            initialPatchIdx = randi( this.numDict );
            rootSpl = SubPatch2dLocation( dimIni, xyzIni, initialPatchIdx, 0 );
            initialConstraintNode = SortedTreeNode( rootSpl );
            
        end
        
        function iniPatchConstraints( this )
%             fprintf('scaleByOverlap: %d\n', this.scaleByOverlap);
            this.pc = PatchConstraints( this.sz3d(1), this.f, this.overlappingPatches, this.scaleByOverlap );
            
            fprintf('Computing constraint matrix inverse...');
            this.pc.buildCmtx();
            fprintf('.done\n');
        end
    end 
    
    methods( Static )
        
        function slc = slice( patch, i, isRow )
        % isRow - if true, sample ith row, 
        %         else sample ith column
            if( isRow )
                slc = patch(i,:);
            else
                slc = patch(:,i);
            end
        end
        
        function slc = slicePadded( patch, i, pad, isRow )
        % isRow - if true, sample ith row, 
        %         else sample ith column
            if( isRow )
                slc = patch(i-pad:i+pad,:);
            else
                slc = patch(:,i-pad:i+pad);
            end
        end
        
        function ds = downsampleVec( vec, f )
        % downsample vector 
            ds = vecToRowCol( vec, 'row' );
            ds = mean(reshape( ds, f, [] ));
        end
        
        function similarity = BADpatchConsistency2d( p1, i1, d1, p2, i2, d2, f )
            
            slc1 = Dict2dTo3d.slice( p1, i1, d1 );
            slc2 = Dict2dTo3d.slice( p2, i2, d2 );
            
            slc1ds = Dict2dTo3d.downsampleVec( slc1, f );
            slc2ds = Dict2dTo3d.downsampleVec( slc2, f );
            
            similarity = norm( slc1ds - slc2ds );
        end
        
        function similarity = patchConsistency2d( p1, i1, d1, p2, i2, d2, f )
            
            half = (f-1)./2;
            slc1 = Dict2dTo3d.slicePadded( p1, i1, half, d1 );
            slc2 = Dict2dTo3d.slicePadded( p2, i2, half, d2 );
            
            dsFun = Tid.sum2dxy();
            
            slc1ds = sum(dsFun( slc1, f ));
            slc2ds = sum(dsFun( slc2, f ));
            
            similarity = norm( slc1ds - slc2ds );
        end
        
        function similarity = patchConsistencyConstraints( constraints, X, sz, f )
           
            N = size( X, 1 );
            similarity = zeros( N, 1 );
            dsFun = Tid.sum2dxy();
            
            for i = 1:N
                
                p = reshape( X(i,:), sz );
                pd = dsFun( p, f );
                notNanInds = ~isnan(constraints);
%                 constraints(notNanInds)
%                 pd(notNanInds)
%                 constraints(notNanInds) - pd(notNanInds)
%                 fprintf('\n\n')
                similarity(i) = norm(constraints(notNanInds) - pd(notNanInds));
                
            end
        end
        
        function psv = patchSumVecs( p, i, isRow, downsampleFactor )
            half = (downsampleFactor-1)./2;
            slc = Dict2dTo3d.slicePadded( p, i, half, isRow );
            dsFun = Tid.sum2dxy();
            psv = dsFun( slc, downsampleFactor );
        end
        
        function sumMtx = allSums( X, sz, f )
            [N,M] = size( X );
            if( prod(sz) ~= M )
                error('size inconsistent with patch matrix X')
            end
            
            half = (f-1)./2;
            numValid = sz - ( 2 * half );
            L = sum( numValid );
            
            sumSz = max(sz./f);
            
            sumMtx = zeros( N*L, sumSz );
            
            k = 1;
            for n = 1:N
                for d = 0:1
                    for ri = 2 : sz(d+1)-1;
                        sumMtx( k, : ) = Dict2dTo3d.patchSumVecs( reshape( X(n,:), sz ), ri, d, f );
                        k = k + 1;
                    end
                end
            end
        end
        
        function intersection = findIntersections( msk, v )
            j = setdiff( 1:3, v(1));
            m = msk( j, : );
            intersection = m( m > 0 );
        end
        
        function [coords, logicalIdx] = findIntersectionsList( xyzDimList, dim, xyz )
            logicalIdx = (xyzDimList(:,2) > 0 ) & (xyzDimList(:,1) ~= dim);
            coords = xyzDimList( logicalIdx, : );
        end
                
        % xyz{1,2} are 3-vectors giving a point planes 1 and 2
        % n{1,2} are 3-vectors giving the normal vectors 
        function [msk,d] = planeIntersections( sz, xyz1, n1, xyz2, n2 )
        %function [msk,d,ix,iy,iz] = planeIntersections( sz, xyz1, n1, xyz2, n2 )
            
            msk1 = Dict2dTo3d.planeMask( sz, xyz1, n1);
            msk2 = Dict2dTo3d.planeMask( sz, xyz2, n2);
            msk = msk1 & msk2;

            if( nargout == 1)
                return;
            end
            
            if( n1 == n2 )
                d = 0;
            else
                d = setdiff( 1:3, [n1 n2] );
            end
            
            %if( nargout > 2 )
            %    [x,y,z] = meshgrid(1:sz(1), 1:sz(2), 1:sz(3));
            %    ix = x(msk);
            %    iy = y(msk);
            %    iz = z(msk);
            %end
            
        end

        function line = patchIntersectionFromMask( patch, xyz, n, msk )
            sz = size( msk );
            pmsk = Dict2dTo3d.planeMask( sz, xyz, n );
            line = patch( msk( pmsk )); 
        end
        
        function I = fill3dWith2d( sz3d, dim, xyz, f, patch)
            I = zeros( sz3d );
            msk = Dict2dTo3d.planeMaskF( sz3d, xyz, dim, f);
            I( msk > 0 ) = patch( msk(msk>0) );
        end
        
        % n is {1,2,3}
        function msk = planeMask( sz, xyz, n )
            msk = false( sz );
            
            if( isscalar(xyz) )
                val = xyz;
            else
                switch n
                    case 1
                        val = xyz(1);
                    case 2
                        val = xyz(2);
                    case 3
                        val = xyz(3);
                    otherwise
                        error('invalid normal direction');
                end
            end
            
            switch n
               case 1
                   msk( val, :, : ) = true;
               case 2
                   msk( :, val, : ) = true;
               case 3
                   msk( :, :, val ) = true;
               otherwise 
                   error('invalid normal direction');
            end
        end
        
        % n is {1,2,3}
        function msk = planeMaskF( sz, xyz, n, f, centered )
            
            msk = zeros( sz );
            
            if( ~exist('centered','var') || isempty( centered ))
               centered = false; 
            end
            
            if( isscalar(xyz) )
                val = xyz;
            else
                switch n
                    case 1
                        val = xyz(1);
                    case 2
                        val = xyz(2);
                    case 3
                        val = xyz(3);
                    otherwise
                        error('invalid normal direction');
                end
            end
            
            if( centered )
                half = (f-1)./2;
                rng = val-half : val+half;
            else
                rng = val : val + f - 1;
            end
            
            valid = (rng>0) & (rng<=sz(1));
            rng = rng( valid );

            N = prod(sz(1:2));
            
            v = repmat( reshape(1:N, sz(1:2)), [1 1 nnz(valid)]);
            
            switch n
                case 1
                    msk( rng, :, : ) = permute(v, [3 1 2]);
                case 2
                    msk( :, rng, : ) = permute(v, [1 3 2]);
                case 3
                    msk( :, :, rng ) = v;
                otherwise
                    error('invalid normal direction');
            end
        end
        
        % Deprecated
        function intersection = findIntersectionsBig( msk, v )
           n = v(4);
           switch n
               case 1
                   m = msk( v(1), :, :, [2 3] );
               case 2
                   m = msk( :, v(2), :, [1 3] );
               case 3
                   m = msk( :, :, v(3), [1 2] );
               otherwise 
                   error('invalid normal vector');
           end
           intersection = m( m > 0 );
        end
        
        function simMtx = allSimilaritiesSums( X, sz ) 
        % X must be (N x M)
        % where N is the number of patches and
        % M is prod( sz )
        % and sz is a vector describing the patch size
            [N,M] = size( X );
            if( prod(sz) ~= M )
                error('size inconsistent with patch matrix X')
            end
            half = (f-1)./2;
            numValid = sz - ( 2 * half );
            L = sum( numValid );
            
            simMtx = zeros( N , L );
        end
        
        function [cmtx, idxs ] = contraintsMtx( constraintPtMask, pm )
            
            pm = pm( constraintPtMask );
            
            % remove zero
            % since it indicates points not in the patch
            patchIdxsInMask = setdiff( unique( pm ), 0 );
            
            N = length( patchIdxsInMask );  % the number of constraints
            M = nnz( constraintPtMask ); % the number of elements in the HR patch
            cmtx = zeros( N , M );
            
            for i = 1:N
                cmtx( i,  (pm == patchIdxsInMask(i)) ) = 1;
            end
            idxs = pm( patchIdxsInMask );
        end
        
        function [cmtx, b] = contraintsFromMask( patch, constraintPtMask, pm )
            
            pm = pm( constraintPtMask );
             
            % remove zero 
            % since it indicates points not in the patch
            patchIdxsInMask = setdiff( unique( pm ), 0 );
            
            N = length( patchIdxsInMask );  % the number of constraints
            M = nnz( constraintPtMask ); % the number of elements in the HR patch
            cmtx = zeros( N , M );
            
            for i = 1:N
                cmtx( i,  (pm == patchIdxsInMask(i)) ) = 1;
            end
            
            b = patch(  pm( patchIdxsInMask ));
        end
        
        function [ xsectParentList ] = intersectingParents( dimXyzNode )
            
            depth = dimXyzNode.getDepth();
            xsectParentList = zeros( depth, 3 );
            parent = dimXyzNode.getParent();
            n = 0;
            for i = 1:depth
                
                if( dimXyzNode.getData().dim ~= parent.getData().dim )
                    n = n + 1;
                    xsectParentList( n, : ) = [ parent.getData().dim, ...
                                                parent.getData().xyz,...
                                                parent.getData().idx ];
                end
                
                parent = parent.getParent();
            end
            
            xsectParentList = xsectParentList(1:n,:);
            
        end

        function [node, root ] = patchParamArrayToNode( paramArray )

            import net.imglib2.algorithms.patch.*;
            import net.imglib2.algorithms.opt.astar.*;

            parentNode = []; % equivalent to java null

            for i = 1:size(paramArray,1)

                node = SortedTreeNode( ...
                        SubPatch2dLocation( ...
                            paramArray(i,1), paramArray(i,2),...
                            paramArray(i,3), 0 ), ...
                        parentNode);
                
                parentNode = node;
                
                if( i == 1 )
                    root = node;
                end
            end
            
        end
        
        function params = patchParamNodeToArray( paramNode )
            
            params = [  paramNode.getData().dim, ...
                        paramNode.getData().xyz, ...
                        paramNode.getData().idx ];
                    
            next = paramNode;
            while( ~next.isRoot())
                next = next.getParent();
                params = [ params; ...
                    next.getData().dim, ...
                    next.getData().xyz, ...
                    next.getData().idx ];
                
            end
        end
        
        function iterCounts = iterationCounts_dist( f_out )
        % outputs the final iteration count for each distributed buildPatch3d job
            N = length( f_out );
            iterCounts = zeros( N, 1);
            for i = 1:N
                iterCounts(i) = find( ( f_out{i}{1}{4} > -1 ), 1, 'last');
            end
        end
        
    end % static methods
    
end % class
