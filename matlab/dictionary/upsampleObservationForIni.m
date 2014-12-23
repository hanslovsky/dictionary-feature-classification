function im_re = upsampleObservationForIni( im, sz, usFactorZ, interpOpts )
% Usage:
%   im_re = upsampleObservationForIni( im, sz, usFactorZ, interpOpts )

if( ~exist('interpOpts','var'))
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

if( ~isempty( interpOpts ) && any(strcmp( interpOpts, 'nearest' )))
    fprintf('UPSAMPLING NEAREST-NEIGHBOR');
	zSamples = repmat(1:sz(3), usFactorZ , 1);
    im_re = im( :, :, zSamples );
else
    [ xo, yo, zo ] = meshgrid( xrng, yrng, zrng * half );
    [ xn, yn, zn ] = meshgrid( xrng, yrng, xrng );
    im_re = interp3( xo, yo, zo, im, xn, yn, zn, interpOpts{:} );
end