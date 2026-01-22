% Upsample and perform direct quadrature of stresslet kernel
function u = stresslet_with_upsamp( xt, density, body, fac)
% Upsampling
fgrid = upsample(body, fac);
% Interpolation of density
fdensity = fgrid.interpolate(density);
% Quadrature weight and normal combined
fnorm = bsxfun(@times, fgrid.n, fgrid.w);
% Compute via direct summation
u = stresslet(xt, fgrid.x, fnorm, fdensity);
end
