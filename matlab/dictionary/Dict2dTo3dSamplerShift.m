classdef Dict2dTo3dSamplerShift < Dict2dTo3dSampler
    % Dict2dTo3dSamplerShift
    %
    % John Bogovic
    % HHMI
    % January 2015
    
    properties( SetAccess = protected )
        D2d_shifted;
        D2d_shifted_downsampled;
    end
    
    methods
        
        function this = Dict2dTo3dSamplerShift( D2d, sz2d, f, pc, scaleByOverlap, comparator )
            
            this = this@Dict2dTo3dSampler( D2d, sz2d, f, pc, scaleByOverlap, comparator );
            if( ~isa( this.comparator, 'PatchCompareShift' ))
                error( 'comparator must be of type ''PatchCompareShift''');
            end
            
            if( this.scaleByOverlap )
                this.paramScales = this.pc.overlapFraction;
            end
            
            this.D2d_shifted = this.pc.buildShiftedDictionary( this.D2d, this.comparator );
            
            this.D2d_shifted_downsampled = cell( this.comparator.numShifts, 1 );
            for i = 1:this.comparator.numShifts
                this.D2d_shifted_downsampled{i} = this.pc.downsample2dDictByDimSlow( this.D2d_shifted{i} );
            end
            this.D2d_downsampled = this.D2d_shifted_downsampled(all( this.comparator.shifts == 0, 2 ));
            
%             this.D2d_downsampled = this.pc.downsample2dDictByDimSlow( this.D2d );
%             this.D2d_shifted_downsampled = this.pc.buildShiftedDownsampledDictionary( this.D2d, this.comparator );
                
        end % constructor
                
        function [ dictIdxs, dictCosts, shifts, shiftIdxs, models, x, pmsk, dmtx ] = bestKdictsShift( this, xin, i, K, isLR  )
            
            if( ~exist('isLR', 'var'))
                isLR = [];
            end
            
            models = {};
            [ x, pmsk, dmtx, omsk ] = this.pc.getSubPatchI( xin, i, isLR );
            x = x(:)';
            
            % here D is cell array where each cell contains
            % the values of the dictionary for a given shift
%             if( numel( x ) == prod( this.pc.sz2d ) )
%                 D = this.D2d_shifted;
%             else
%                 D = this.D2d_shifted_downsampled;
%             end

            Dexp = cell2mat( this.D2d_shifted );
            if( ~isempty( dmtx ))
                disp('here');
                Dexp = Dexp(:, pmsk(:)) * dmtx';
            end

            % check size of K
            K = min( K, size( Dexp, 1 ));
            
            % TODO deal with models if they exist
            
            % compute the distance
            [ dists ] =  this.comparator.distances( x, Dexp );
            [ dictsSorted, is ] = sort( dists );

            dictIdxsWSort = is( 1:K ) % this has the sort index mixed in
            dictCosts = dictsSorted(1:K);
            
            dictIdxs = mod( dictIdxsWSort, this.numDict );
            dictIdxs( dictIdxs == 0 ) = this.numDict;
            
            shiftIdxs = 1 + ( (dictIdxsWSort - dictIdxs)./ this.numDict );
            shifts = this.comparator.shifts( shiftIdxs, : );
        end
        
        function [ dictIdxs, dictCosts, shifts, shiftIdxs, models, x, msk, didTp ] = bestKdictsShift1( this, xin, i, K, isLR  )
            
            if( ~exist('isLR', 'var'))
                isLR = [];
            end
            
            x = vecToRowCol( xin, 'row');
            
            models = {};
            
            [x, ~, msk, ~, didTp ] = this.pc.getSubImage( x, i, isLR );
            outSz = size( x );
%             [ x, pmsk, dmtx, omsk ] = this.pc.getSubPatchI( x, i, isLR );
            x = x(:)';
            
            % here D is cell array where each cell contains
            % the values of the dictionary for a given shift
            if( numel( x ) == prod( this.pc.sz2d ) )
                D = this.D2d_shifted;
            else
                D = this.D2d_shifted_downsampled;
            end
            Dexp = cell2mat( D );
            
            % check size of K
            K = min( K, size( Dexp, 1 ));
            
            % TODO deal with models if they exist
            
            % compute the distance
            [ dists ] =  this.comparator.distances( x, Dexp );
            [ dictsSorted, is ] = sort( dists );

            dictIdxsWSort = is( 1:K ) % this has the sort index mixed in
            dictCosts = dictsSorted(1:K);
            
            dictIdxs = mod( dictIdxsWSort, this.numDict );
            dictIdxs( dictIdxs == 0 ) = this.numDict;
            
            shiftIdxs = 1 + ( (dictIdxsWSort - dictIdxs)./ this.numDict );
            shifts = this.comparator.shifts( shiftIdxs, : );
        end
        
        function [ dictParams, xhat, dictCost, xhatup, shift ] = dictParamsLasso( this, xin, loci, params, isLR )
            % uses the spams toolbox's implementation of lasso to estimate
            % a sparse coefficient vector that explains the observation 'x'
            % at location 'loci' using the 2d dictionary
            
            %             msk = this.pc.planeMaskLRI( loci );
            %             [x,~,didTp] = PatchConstraints.downsampleByMaskDim( x(:), msk );
            
            [ x, pmsk, dmtx, omsk ] = this.pc.getSubPatchI( xin, loci, isLR );
            outSz = size( x );
            x = x(:);
            
            DshiftList = this.D2d_shifted;
            
            dictCost   = inf;
            dictParams = [];
            shift      = [];
            shiftIdx   = 0;
            
            for i = 1 : this.comparator.numShifts
                
                if( isempty( dmtx ))
                    D = DshiftList{ i }(:,pmsk(:));
                else
                    D = DshiftList{ i }(:,pmsk(:)) * dmtx';
                end
                dictParamsTmp = mexLasso( x , D', params );
                
                % evaluate
                xhatTmp = D' * dictParamsTmp;
                dictCostTmp =  norm( x - xhatTmp );
                
                if( dictCostTmp < dictCost )
                    dictCost   = dictCostTmp;
                    dictParams = dictParamsTmp;
                    xhat       = xhatTmp;
                    shift      = this.comparator.shifts( i, : );
                    shiftIdx = i;
                end
            end % loop over shifts
            
            if( nargout > 3 )
                if( isLR )
                    xhatup =  (DshiftList{ i }(:,pmsk(:)))' * dictParams;
                else
                    % if the observation is high res,
                    % then the current xhat is already high res
                    xhatup = xhat;
                end
            end
            
        end
        
        function [ dictParams, xhat, dictCost, xhatup, shift ] = dictParamsLasso1( this, xin, loci, params, isLR )
            % uses the spams toolbox's implementation of lasso to estimate
            % a sparse coefficient vector that explains the observation 'x'
            % at location 'loci' using the 2d dictionary
            
            x = vecToRowCol( xin, 'row');
            %             msk = this.pc.planeMaskLRI( loci );
            %             [x,~,didTp] = PatchConstraints.downsampleByMaskDim( x(:), msk );
            [ x, sz, ~, ~, ~ ] = this.pc.getSubImage( x, loci, isLR );
            outSz = size( x );
            x = x(:);
            
            isLowRes = false;
            if( numel( x ) == prod( this.pc.sz2d ) )
                DshiftList = this.D2d_shifted;
            else
                DshiftList = this.D2d_shifted_downsampled;
                isLowRes = true;
            end
            
            dictCost   = inf;
            dictParams = [];
            shift      = [];
            
            for i = 1 : this.comparator.numShifts
                
                D = DshiftList{ i };
                dictParamsTmp = mexLasso( x , D', params );
                
                % evaluate
                xhatTmp = D' * dictParamsTmp;
                dictCostTmp =  norm( x - xhatTmp );
                
                if( dictCostTmp < dictCost )
                    dictCost   = dictCostTmp;
                    dictParams = dictParamsTmp;
                    xhat       = xhatTmp;
                    shift      = this.comparator.shifts( i, : );
                end
            end % loop over shifts
            
            if( nargout > 3 )
                if( isLowRes )
                    xhatup = this.D2d' * dictParams;
                else
                    % if the observation is high res,
                    % then the current xhat is already high res
                    xhatup = xhat;
                end
            end
            
        end % dictParamsLasso
    
        
        
    end
end
