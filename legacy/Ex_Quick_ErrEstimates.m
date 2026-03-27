% Ex_Quick_ErrEstimates.m: Script to showcase the use of precomputation and
% use of quick error estimates for classifying target points for 
% computing the stresslet kernel. Components of code previously developed 
% by Ludvig af Klinteberg and Chiara Sorgentone
%
% Authors: 
%   Pritpal 'Pip' Matharu, David Krantz
%
% Max Planck Institute for Mathematics in the Sciences
% KTH Royal Institute of Technology
% Date: 2026/01/12
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('src'));
clear all
%% Parameters
% Density flag
% 1: Analytically defined density
% Otherwise, use stresslet identity
dflag = 0;
% Tolerance for method
TOL = 1e-6;
% Upsampling factors
upsamp_fac = 2:3;
% Discretization points
nth = 40;
nph = 60;
% Particle shape
geometry.shape = 'spheroid';
% geometry.shape = 'peanut';

% Interpolate theta root? If not, determine analytically (spheroid) or newton
rflag = true;

FS = 16; % Font size on plots
%% Setup particle
body = init_axsymbody(geometry, nph, nth);

%% Precompute uniform error estimate and create interpolant
errestinterp = uniform_error_estimates(body, upsamp_fac, rflag);

%% Density
qb = density_selection(body, dflag);

%% Evaluation grid
xpts = xgrid_selection(body);

%% Categorize points using quick estimate
disp('Categorize points for upsampling regions')
tic
[precomp, estplt, tmpR, phi0_thstar, intpind_p] = target_precompute( xpts.xout, qb, TOL, body, errestinterp);
toc

%% Compute
disp('Direct and Upsampling Quadrature')
uquad = zeros(size(xpts.xout));
% Direct and Upsampling region
for k=1:numel(precomp.upsampfac)
  % Upsample the density and compute directly
  uu = stresslet_with_upsamp( xpts.xout(precomp.upsamp_mask{k}, :), ...
    qb, body, precomp.upsampfac(k));
  uquad(precomp.upsamp_mask{k},:) = uu;
end
disp('Special Quadrature')
% Special Quadrature Region (computed using "high" upsampling)
uquad(precomp.SQ_mask,:) = stresslet_with_upsamp( xpts.xout(precomp.SQ_mask, :), qb, body, 10);

% Add points back with non-computed values
u = zeros(size(xpts.xt));
u(xpts.mask, :) = uquad;

%% Plotting

% Reference solution
switch dflag
  case 1 % Analytically provided density, use "high" upsampling
    uref = stresslet_with_upsamp( xpts.xt, qb, body, 10);
  otherwise % Stresslet identity, 0 solution
    uref = 0*u;
end

% Error plots
plot_quaderrs(precomp, body, xpts, u, uref,  FS)
% Estimates plots
upfac = 1;
plot_estcontour(xpts, body, qb, estplt, upfac, uref, FS)

return
%% ****************************** FUNCTIONS ******************************

%% xgrid_selection: Set of target points, based on body selected
function xpts = xgrid_selection(body)
switch lower(body.geometry.shape)
  case 'spheroid'
    M = 10; L = 0.2;
    xv = linspace(-L/2, L/2, M);
    zv = linspace(-L, L, M);
    [X,Z] = meshgrid(xv, zv);
    yp = 0.0; Y = yp*ones(size(X));
    xt = [X(:), Y(:), Z(:)];
  case 'peanut'
    M = 200; L = 5;
    xv = linspace(-L/1.5, L/1.5, M);
    zv = linspace(-L, L, M);
    [X,Z] = meshgrid(xv, zv);
    yp = 0.0; Y = yp*ones(size(X));
    xt = [X(:), Y(:), Z(:)];
end

% Determine which points lie outside of particle body
[x_ext, mask_ext] = find_exterior_points(body, xt);
% Outpute
xpts.xt = xt;
xpts.xv = xv;
xpts.zv = zv;
xpts.X = X;
xpts.Z = Z;
xpts.xout = x_ext;
xpts.mask = mask_ext;
end

%% density_selection: Density function, based on inputted density flag
function density = density_selection(body, dflag)
switch dflag
  case 1 % Analytic Density
    sigma_1 = @(theta,phi) 2.1+sin(8*theta)+sin(12*phi);
    sigma_2 = @(theta,phi) 2+sin(6*theta).*cos(6*phi);
    sigma_3 = @(theta,phi) sin(5*theta).*exp(-cos(phi).^2)+1.03;
    q1 = sigma_1(body.theta, body.phi);
    q2 = sigma_2(body.theta, body.phi);
    q3 = sigma_3(body.theta, body.phi);
    density = [q1(:), q2(:), q3(:)];
  otherwise % Stresslet Identity
    % Number of points on particle
    Nb = size(body.x, 1);
    density = [300*ones(Nb,1), 100*ones(Nb,1), 200*ones(Nb,1)];
end

% Simple test to see if density is well resolved
qtmp = density;
errM = 0.0;
for comp = 1:3
  qmat = reshape(qtmp(:, comp), body.nth, body.nph);
  for ii = 1:body.nth
    qhat = fft(qmat(ii, :).')/(body.nph/2);
    qhat = qhat(1:body.nph/2); % Assume density is real value
    qhatM = max(abs(qhat));
    err = max(abs(qhat(end-4:end)./qhatM));
    errM = max(errM, err);
    if err > 1e-12
      fprintf('DENSITY MAY NOT BE RESOLVED!!! Err = %g\n', err)
    end
  end
end
fprintf('Max Err in Density Resolution = %g\n', errM)

end

%% target_precompute: Categorize quadrature required using quick estimate
function [precomp, pltout, tmpR, phi0_thstar, intpind_p] = target_precompute( xt, density, TOL, body, errest)

% Determine closest particle points to each target point
tmpR = knnsearch(body.x, xt);

% Error estimates for determining quadrature method
xycoord = sqrt(xt(:, 1).^2 + xt(:, 2).^2);
zcoord  = abs(xt(:, 3));
% Interpolate estimate, no upsampling (index 1)
esttmp  = cellfun(@(c) 10.^(c(xycoord, abs(zcoord))), {errest.errest_INT{1, :}}, 'UniformOutput',false);
estval  = reshape(cell2mat(esttmp), size(xycoord, 1), 3);
% Modify uniform error estimate with density evaluated at root

if errest.flag % Interpolate root in theta
  th_intp = errest.thRe_INT{1}(xycoord, abs(zcoord))+1i*errest.thIm_INT{1}(xycoord, abs(zcoord));
  [qp, qt, phi0_thstar, intpind_p] = Qintp(body, xt, body.a, body.c, reshape(density, body.nth, body.nph, 3), tmpR, th_intp);
else % Use analytic root expressions if available, otherwise Newton solver
  [qp, qt, phi0_thstar, intpind_p] = Qintp(body, xt, body.a, body.c, reshape(density, body.nth, body.nph, 3), tmpR);
end
% Take max of density in both directions and modify uniform estimate
Mqpqt = max(qp, qt);
Errest = max(estval.*Mqpqt, [], 2);

% For plotting
pltout.Errest_uni    = estval;
pltout.estup_stor{1} = Errest;

% Mask for direct quadrature based on error estimates
precomp.direct_mask = Errest < TOL;         % Points to be computed via direct quadrature
precomp.upsamp_mask{1} = precomp.direct_mask;
available           = ~precomp.direct_mask; % Remaining target points

precomp.upsampfac = errest.upsampfac;
% Upsampled error estimate (use previous density evaluations for speed)
for k=2:numel(errest.upsampfac)
  % Use precomputed error estimates to determine masks
  upsamp    = errest.upsampfac(k);
  % Retrieve precomputed upsampled error estimate
  esttmp    = cellfun(@(c) 10.^(c(xycoord, abs(zcoord))), {errest.errest_INT{upsamp, :}}, 'UniformOutput',false);
  estval    = reshape(cell2mat(esttmp), size(xycoord, 1), 3);
  % Modify error estimate
  errest_up = max(estval.*Mqpqt, [], 2);

  % For plotting
  pltout.estup_stor{k} = errest_up;

  % Determine which points can be computed via upsampling factor
  precomp.upsamp_mask{k} = available & (errest_up < TOL);
  available              = available & (~precomp.upsamp_mask{k}); % Update points remaining to compute
end
% Determine points requiring special quadrature
precomp.SQ_mask = available;        % SQ for points for all other points
% Store tolerance
precomp.tol = TOL;
end


