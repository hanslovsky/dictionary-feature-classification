function out = reshapePy( in, varargin )
% Reshapepy reshapes in the manner of python with the last dimension
% varying fastest

paramrev = varargin(end:-1:1);
if( length( paramrev ) == 1 )
   paramrev{1} =  paramrev{1}(end:-1:1);
end

inrs = reshape( in, paramrev{:} );   % reshape 
out = permute( inrs, length(size(inrs)):-1:1);   % switch dimensions
