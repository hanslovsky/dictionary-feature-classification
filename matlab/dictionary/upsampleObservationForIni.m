function im_re = upsampleObservationForIni( im, sz, dsFactorZ )
% Usage:
%   im_re = upsampleObservationForIni( im, sz, dsFactorZ )

if( ~isscalar( sz ))
    error( 'sz input must be scalar');
end

half = (sz-1)./2;
xyrng = -half:half; %#ok<BDSCI>

zHalf = ceil(half/dsFactorZ);
zrngSample = -zHalf : zHalf;

[ xo, yo, zo ] = meshgrid( xyrng, xyrng, zrngSample * dsFactorZ );
[ xn, yn, zn ] = meshgrid( xyrng, xyrng, xyrng );

im_re = interp3( xo, yo, zo, im, xn, yn, zn );