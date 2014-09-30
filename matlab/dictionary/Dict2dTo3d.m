classdef Dict2dTo3d < handle
    % Tid - transformationally invariant dictionary
    %
    % John Bogovic
    % HHMI
    % August 2014
    
    properties ( SetAccess = private )
        
        D2d;      % The dictionary elements
        numDict;  % Number of 2d dictionary elements
        sums2d;   % the dictionary element sums
        
        sz2d;     % the size of 2d patches
        sz3d;     % the size of 3d patches
        
        f;        % the downsampling factor ( resolution in z of patches )
        
        summer2Dxy; % downsampling by sum function
        
        ndims = 3; 
    end
    
    methods
        
        % Constructor
        % D2d - 
        % sz  - size of 2d patches
        % f   - downsampling factor
        function this = Dict2dTo3d( D2d, sz, f )
            
            this.D2d     = D2d;
            this.numDict = size( this.D2d, 1 );
           
            this.f    = f;
            
            if( length(sz ) > 2 )
               error('sz must be a 2-vector or scalar');
            end
            if( length(sz ) == 2 )
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
            
            this.sums2d = Dict2dTo3d.allSums( this.D2d, this.sz2d, this.f);
            
            this.summer2Dxy = Tid.sum2dxy();
        end
        
        function patch = sample2dTo3d( this )
            
%             [N,M] = size( this.D2d );
            szSm = this.sz2d ./ this.f;
            
            % hard-coded '3' corresponds to 3 spatial-dimensions
            msk = zeros( [ 3 max(szSm) ] );
            mskN = numel(msk);
            xyzn = randperm( mskN );
            
            patchIdxs = zeros( [ 3 max(szSm) ] );
            
            for ii = 1:mskN
                
                [d,n] = ind2sub( size(msk), xyzn(ii));
                
                if( ii > 1 )
                    xsectList = Dict2dTo3d.findIntersections( msk, [d n]);
                    xsectList
                else
                   patchIdxs(ii) = randi( this.numDict, 1 ); 
                end
                
                msk( xyzn(ii) ) = ii
                patchIdxs(ii)
                
                pause;
            end
            
            patch = [];
        end
        
        function constraints = getSumConstraints( this, patchIdxs )
            
            szdown = this.sz3d ./ this.f;
            constraints = nan( szdown );
            
            [d,xyz] = ind2sub( size( patchIdxs ), find( patchIdxs > 0 ));
            for i = 1:length(xyz)
                patch = reshape( this.D2d( patchIdxs(d(i),xyz(i)), : ), this.sz2d );
                constraints( Dict2dTo3d.planeMask( szdown, xyz(i), d(i) )) = this.summer2Dxy( patch, this.f );
%                 constraints
%                 fprintf('\n\n');
            end
        end
        
        function c = collectConstraints( this, constraints, xyz, d )
            szSm3 = this.sz3d ./ this.f;
            szSm2 = this.sz2d ./ this.f;
            planeMask = Dict2dTo3d.planeMask( szSm3, xyz, d );
%             planeMask
            c = reshape( constraints( planeMask ), szSm2 );
        end
        
    end
    
    methods( Static )
        
        % isRow - if true, sample ith row, 
        %         else sample ith column
        function slc = slice( patch, i, isRow )
            if( isRow )
                slc = patch(i,:);
            else
                slc = patch(:,i);
            end
        end
        
        % isRow - if true, sample ith row, 
        %         else sample ith column
        function slc = slicePadded( patch, i, pad, isRow )
            if( isRow )
                slc = patch(i-pad:i+pad,:);
            else
                slc = patch(:,i-pad:i+pad);
            end
        end
        
        % downsample vector 
        function ds = downsampleVec( vec, f )
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
                
        % xyz{1,2} are 3-vectors giving a point planes 1 and 2
        % n{1,2} are 3-vectors giving the normal vectors 
        function [msk,d] = planeIntersections( sz, xyz1, n1, xyz2, n2 )
        %function [msk,d,ix,iy,iz] = planeIntersections( sz, xyz1, n1, xyz2, n2 )
            
            msk1 = Dict2dTo3d.planeMask( sz, xyz1, n1);
            msk2 = Dict2dTo3d.planeMask( sz, xyz2, n2);
            msk = msk1 & msk2;

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
        
        % X must be (N x M)
        % where N is the number of patches and
        % M is prod( sz )
        % and sz is a vector describing the patch size
        function simMtx = allSimilarities( X, sz ) 
            [N,M] = size( X );
            if( prod(sz) ~= M )
                error('size inconsistent with patch matrix X')
            end
            half = (f-1)./2;
            numValid = sz - ( 2 * half );
            L = sum( numValid );
            
            simMtx = zeros( N , L );
        end
        
    end
    
end
