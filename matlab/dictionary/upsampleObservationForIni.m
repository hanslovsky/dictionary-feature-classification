function im_re = upsampleObservationForIni( im, sz, usFactorZ, interpOpts )
% Usage:
%   im_re = upsampleObservationForIni( im, sz, usFactorZ, interpOpts )

if( ~exist('interOpts','var') || ~iscell(interOpts))
    interpOpts = {};
end

if( isscalar( sz ))
    half = (sz-1)./2;
    xyrng = -half:half;
    
    zHalf = ceil(half/usFactorZ);
    zrng = -zHalf : zHalf;
    
    xrng = xyrng;
    yrng = xyrng;
else
    xhalf = (sz(1)-1)./2;
    yhalf = (sz(2)-1)./2;
    half = xhalf;
        
    xrng = -xhalf:xhalf; %#ok<*BDSCI>
    yrng = -yhalf:yhalf;
    
    zHalf = (sz(3)-1)./2;
    zrng = -zHalf : zHalf;
end


[ xo, yo, zo ] = meshgrid( xrng, yrng, zrng * half );
[ xn, yn, zn ] = meshgrid( xrng, yrng, xrng );

im_re = interp3( xo, yo, zo, im, xn, yn, zn, interpOpts{:} );