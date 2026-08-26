% Precomputation of error estimate, using uniform density
function errestinterp = uniform_error_estimates(body, upsampfacs, iflag)

disp('Precomputing Quadrature Error Estimate')

% Set upsampling factors and storage cell for estimate
upfacvals = unique([1, upsampfacs]); % Ensure 1 is included (direct quad)
est_INT = cell(length(upfacvals),3);

% If interpolating root in theta on-the-fly
if iflag
  % Storage for theta roots
  Rerts_INT = cell(length(upfacvals),1);
  Imrts_INT = cell(length(upfacvals),1);
end

%% Setup tabulation grid and estimate values
% Set up evaluation grid
nt = 100;
nb = 100;
% Setup interpolation based on distance from z axis and in z
xy = abs(linspace(0, 2*(body.a+body.c), nt-1));
xy = [-xy(2) xy]; % Add additional point, to avoid issues on symmetry axis
z  = abs(linspace(0, 2*(body.a+body.c), nb-1));
z = [-z(2) z];    % Add additional point, to avoid issues on symmetry axis
[xy_grid, z_grid] = ndgrid(xy, z);
% Create grid for interpolant
tmpxy  = reshape(xy_grid, nt, nb);
tmpz   = reshape(z_grid, nt, nb);

% Set Coordinates (values of x and y must match distance values)
xintp = xy;
y     = 0;
[X, Y, Z] = ndgrid(xintp, y, z);
x_eval    = [X(:), Y(:), Z(:)];

% Determine points outside of spheroid, to compute error est
[x_exteval, evalmask] = find_exterior_points(body, x_eval);

% Determine point on the axis of symmetry outside of particle (y values should all be equal to zero)
xext_xy = find(x_eval(:, 1) == 0 & x_eval(:, 2) == 0 & evalmask);
% ** Find the point along the axis that is closest to particle surface
[~, maxin_ind] = min(x_eval(xext_xy, 3));

%% Loop through all upsampling factor values
fprintf('Precomputing Upsampled Error Estimate, from kappa = %d, ..., %d \n', upfacvals(1), upfacvals(end))
for j = 1:length(upfacvals)
  upfac = upfacvals(j);
  % Store values from body or upsampled body full estimate
  refined_body = upsample(body, upfac);

  fprintf('Precomputing Upsampled Error Estimate, kappa = %d \n', upfac)
  Esttmp = zeros(size(x_exteval, 1), 3); % Values of error estimate outside of particle
  Errest = zeros(size(x_eval, 1), 3);    % Values of error estimate everywhere
  th_tmp = zeros(size(x_exteval, 1), 1); % Values of theta root

  % Indexes of closest points
  ind_R = knnsearch(refined_body.x, x_exteval);
  % Only loop through points outside of particle
  for ii =1:size(x_exteval, 1)
    [Esttmp(ii, :), th_tmp(ii)] = errest_unifdens(x_exteval(ii, :), ind_R(ii), body.at, body.dadt, body.ct, body.dcdt, refined_body);
  end
  % Use stresslet identity density for precomputing error est
  Errest(evalmask, :)   = Esttmp;
  Errest(~evalmask, :)  = 1/eps^2;        % Points inside spheroid, set (arbitrarily) to 1/eps^2. *Note: Magnitude must be large enough to offset potential gaps between surface and interpolation points
  Errest(isnan(Errest)) = 1/eps^2;        % Any NaNs are also set to 1/eps^2
  % **Safety measure since tabulation grid is not adaptive, ensure that we do not miscategorize evaluation
  Errest(xext_xy(maxin_ind), :)  = 1/eps^2;      % Error estimate cannot be trusted on the axis of symmetry, closest to particle
  Errest                = Errest+eps; % Set lower bound on error estimate as eps (avoids value being exactly 0)

  % Loop through each component
  for iii = 1:3
    tmpest_up = reshape(Errest(:,iii), nt, nb);
    % Create gridded interpolant for error estimate
    est_INT{j, iii}     = griddedInterpolant(tmpxy, tmpz, log10(tmpest_up), 'linear');
  end

  % Save theta roots, if interpolating on-the-fly
  if iflag
    theta0 = zeros(size(x_eval, 1), 1);    % Values of error estimate everywhere
    theta0(evalmask, :)   = th_tmp;
    tmprts = reshape(theta0, nt, nb);
    % Create gridded interpolant for roots
    % Separate real and imaginary parts
    Rerts_INT{j} = griddedInterpolant(tmpxy, tmpz, real(tmprts), 'linear');
    Imrts_INT{j} = griddedInterpolant(tmpxy, tmpz, abs(imag(tmprts)), 'linear');
  end
end

%% Store and Save
% Store gridded interpolant
errestinterp.errest_test = Errest;
errestinterp.errest_INT = est_INT;
errestinterp.XY_grid = xy_grid; % Added for plotting figures
errestinterp.Z_grid  = z_grid;  % Added for plotting figures
errestinterp.upsampfac = upfacvals;

% Store theta root interpolation
if iflag
  errestinterp.flag     = iflag;
  errestinterp.thRe_INT = Rerts_INT;
  errestinterp.thIm_INT = Imrts_INT;
else
  errestinterp.flag     = false;
end

end

%% ****************************** FUNCTIONS ******************************
%% errest_unifdens: Determine the full error estimate for axisymmetric body
% Use linear map t=2*theta/pi-1 for theta root
function [est, theta0phistar, t0phistar, phi0tstar, t_star, theta_star, s_star]=...
  errest_unifdens(xt, ind_R, a, dadt, c, dcdt, bodyvals)
% est -------------- combined trapezoidal and Gauss-Legendre error estimate
% theta0phistar ---- root in theta direction
% phi0tstar -------- root in phi direction
% t_star ----------- value of t for closest grid point
% s_star------------ value of phi for closest grid point
%*** We use the linear mapping in t=2*theta/pi-1 when determining the root
% in theta direction, and the value of theta for all calculations and
% estimates otherwise  ***%
[theta_star, s_star, t_star, imap, dfac, phi0tstar, t0phistar, theta0phistar, Gp0t, Gt0p] = ...
  axsym_roots(bodyvals, ind_R, xt, a, c, dadt, dcdt, bodyvals.geometry.shape);
% We use the stresslet
p = 5/2;

%% Find linear approximation
% Value of the surface at the point on the body
srfvec = axsym_vector(theta_star,s_star,a,c);
% Difference between point on surface and evalution point
rvec = srfvec-xt;
% Semi-analytical estimates
dp0fun = semiroot(rvec, axsym_drdphi(theta_star,s_star,a),           axsym_drdtheta(theta_star,s_star,dadt,dcdt)); % phi
dt0fun = semiroot(rvec, axsym_drdtheta(theta_star,s_star,dadt,dcdt), axsym_drdphi(theta_star,s_star,a)); % theta

%% Set up integration in both directions
% Will use Gauss-laguerre quadrature to capture error estimates
[xlag,wlag] = gausslaguerre8();
% estimate ~ exp(-c*t), where

% Trapezoidal rule
grad_perp_TZ = axsym_drdtheta(theta_star,s_star,dadt,dcdt);
grad_TZ      = axsym_drdphi(theta_star,s_star,a);
C_TZ         = bodyvals.nph*norm(grad_perp_TZ)/norm(grad_TZ);
% Gauss-Legendre quadrature
grad_perp_GL = axsym_drdphi(theta_star,s_star,a);
grad_GL      = axsym_drdtheta(theta_star,s_star,dadt,dcdt);
C_GL         = 2*(bodyvals.nth*norm(grad_perp_GL)/norm(grad_GL));


%% Phi direction (Trapezoidal rule)

% Check both positive and negative imaginary part for root, and take smaller error
estTZ = zeros(2, 3);
phi0tmp = [phi0tstar; conj(phi0tstar)];
for ii = 1:2
  phi0tstar = phi0tmp(ii);
  % Determine estimate in each component
  for iii =1:3
    % Use stresslet identity density for precomputing error estimate
    density         = zeros(1, 3);
    density(:, iii) = 1.0;

    % Density evaluated at root in phi with fixed theta
    % When using the SI identity for precomputations, avoid interpolation
    f_s = density;
    % If using non-uniform, discrete density
    % f_s = localintp_est_phi0(bodyvals, density, phi0tstar, theta_star);

    % Evalute kernel at root
    fp0t = kernel_phi0(p, xt, f_s, theta_star, phi0tstar, a, c, dadt, dcdt);

    % Numerical integration quantities
    est_const = abs( 2/gamma(p)*fp0t*Gp0t^(-p) );
    ds        = xlag/C_TZ;
    dtpos = phi0tstar - dp0fun(ds)  + dp0fun(0);
    dtneg = phi0tstar - dp0fun(-ds) + dp0fun(0);
    estintegrand_pos = trapz_errfunc_deriv(dtpos, bodyvals.nph, p-1);
    estintegrand_neg = trapz_errfunc_deriv(dtneg, bodyvals.nph, p-1);
    % Collect quantities to determine error
    estTZ(ii, iii) = est_const * ...
      sum((estintegrand_pos+estintegrand_neg).*exp(xlag).*wlag/C_TZ);
  end
end
% Error estimate from phi direction contribution (trapzoidal integration)
est1 = min(estTZ, [], 1);
est1(isnan(est1)) = 0.0; % If target point at north/south pole, ensure value not NaN

%% Theta direction (Gauss-Legendre rule)

% Check both positive and negative imaginary part for root, and take smaller error
estGL = zeros(2, 3);
t0tmp     = [t0phistar; conj(t0phistar)];
for ii = 1:2
  t0phistar = t0tmp(ii);

  % Determine estimate in each component
  for iii =1:3
    % Use stresslet identity density for precomputing error estimate
    density         = zeros(1, 3);
    density(:, iii) = 1.0;

    % Density evaluated at root in theta with fixed phi
    % When using the SI identity for precomputations, avoid interpolation
    f_t = density;
    % If using non-uniform, discrete density
    % f_t = localintp_th0(bodyvals, density, theta0phistar, theta_star, s_star);

    % Evalute kernel at root
    ft0p = kernel_t0(p, xt, f_t, t0phistar, s_star, a, c, dadt, dcdt, imap, dfac);

    % Numerical integration quantities
    est_const = abs(2/gamma(p)*ft0p*Gt0p^(-p));
    ds        = xlag/(C_GL);
    dtpos = ( t0phistar - ( dt0fun(ds)  - dt0fun(0) ) );
    dtneg = ( t0phistar - ( dt0fun(-ds) - dt0fun(0) ) );
    estintegrand_pos = gl_errfunc_deriv(dtpos, bodyvals.nth, p-1);
    estintegrand_neg = gl_errfunc_deriv(dtneg, bodyvals.nth, p-1);
    % Collect quantities to determine error
    estGL(ii, iii) = est_const * ...
      sum((estintegrand_pos+estintegrand_neg).*exp(xlag).*wlag/(C_GL));
  end
end
% Error estimate from theta direction contribution (Gauss-Legendre integration)
est2 = min(estGL, [], 1);

%% Final combined estimate
est = est1 + est2;

% Inverse map for root, to theta
t0phistar = imap(t0phistar);

end % End of S3Q_fullintpest

%% kernel_t0: Dependent on kernel (single, double, stresslet), t0 at kernel
function ft0p = kernel_t0(p, xt, f_t, t0phistar, s_star, a, c, dadt, dcdt, imap, dfac)
switch p
  case 1/2 % Single layer
    % Surface element at root
    [srfele, ~, ~] = element_normal(t0phistar, s_star, a, dadt, dcdt, imap, dfac);
    % ft0p           = f_t(t0phistar)*srfele; % Multiply by density at root
    ft0p           = f_t*srfele; % Multiply by density at root
  case 3/2 % Double layer
    % Parameterization of surface, evaluated at root
    srfvec          = axsym_vector(t0phistar, s_star, a, c, imap);
    % Normal vector
    [~, ~, srfnorm] = element_normal(t0phistar, s_star, a, dadt, dcdt, imap, dfac);
    % ft0p            = f_t(t0phistar).*sum((srfvec-xt) .* srfnorm, 2); % Evaluate kernel
    ft0p            = f_t.*sum((srfvec-xt) .* srfnorm, 2); % Evaluate kernel
  otherwise % Stresslet
    srfvec    = axsym_vector(t0phistar, s_star, a, c, imap);
    [~, nvec, ~] = element_normal(t0phistar, s_star, a, dadt, dcdt, imap, dfac);
    r    = all_differences(srfvec, xt);
    rq   = custom_scalar_product(r, f_t); % N×m
    rn   = custom_scalar_product(r, nvec); % N×m
    tmp  = (rq .* rn ); % N×m
    f    = bsxfun(@times, r, tmp); % N×m×3
    %             ft0p = -6*(sum(sum(abs(f))));
    ft0p = -6*((max(abs(f))));
end
end

%% kernel_phi0: Dependent on kernel (single, double, stresslet), phi0 at kernel
function fp0t = kernel_phi0(p, xt, f_s, theta_star, phi0tstar, a, c, dadt, dcdt)
switch p
  case 1/2 % Single layer
    % Surface element at root
    [srfele, ~, ~] = element_normal(theta_star, phi0tstar, a, dadt, dcdt);
    % fp0t           = f_s(phi0tstar)*srfele; % Multiply by density at root
    fp0t           = f_s*srfele; % Multiply by density at root
  case 3/2  % Double layer
    % Parameterization of surface, evaluated at root
    srfvec          = axsym_vector(theta_star,phi0tstar,a,c);
    % Normal vector
    [~, ~, srfnorm] = element_normal(theta_star, phi0tstar, a, dadt, dcdt);
    % fp0t            = f_s(phi0tstar).*sum((srfvec-xt) .* srfnorm, 2); % Evaluate kernel
    fp0t            = f_s.*sum((srfvec-xt) .* srfnorm, 2); % Evaluate kernel
  otherwise % Stresslet
    % Stresslet Kernel
    % %f = -6 * r_i * [r_j q_j(y)] * [r_k n_k(y)] / r^5,
    srfvec    = axsym_vector(theta_star,phi0tstar,a,c);
    [~, nvec, ~] = element_normal(theta_star, phi0tstar, a, dadt, dcdt);
    r    = all_differences(srfvec, xt);
    rq   = custom_scalar_product(r, f_s); % N×m
    rn   = custom_scalar_product(r, nvec); % N×m
    tmp  = (rq .* rn ); % N×m
    f    = bsxfun(@times, r, tmp); % N×m×3
    %             fp0t = -6*(sum(sum(abs(f))));
    fp0t = -6*(max(abs(f)));
end
end

%% axsym_roots: Compute roots of axissymmetric particle
function [theta_star, s_star, t_star, imap, dfac, phi0tstar, t0_phistar, ...
  theta0_phistar, Gp0t, Gt0p] = ...
  axsym_roots(bodyvals, ind_R, xeval, a, c, dadt, dcdt, shape)
%% Find index of closest quadrature node (on body) to evaluation point
[I1,I2]    = ind2sub([bodyvals.nth,bodyvals.nph], ind_R);
theta_star = bodyvals.theta(I1);        % theta value of point on body  [pi, 0]
s_star     = bodyvals.phi(I2);          % phi value of point on body    [0, 2pi]
t_star     = (2/pi)*theta_star - 1; % Linear map in t, Full Sphere  [-1, 1]
imap       = @(t) (pi/2)*(t + 1);   % inverse linear mapping, Full Sphere
dfac       = (pi/2);                % factor from derivative chain rule

%% Approx. phi0star (root in phi direction)
% Use newton solver, as it produces slightly better results down to
% numerical precision for the root than the analytic expression
gamma_tilde_S  = @(s)     axsym_vector(theta_star, s, a, c);
dgamma_tilde_S = @(s)     axsym_drdphi(theta_star, s, a);
% We make a small imaginary perturbation of the initial guess to ensure the
% Newton solver finds the complex root
%[phi0tstar, Gp0t] = f_newton(gamma_tilde_s, dgamma_tilde_s, xt, s_star+1i*0.1);

%% Analytically determine phi roots
% Expressions for finding roots analytically
lambda_PR = (1./(2*a(theta_star))).*((a(theta_star).^2 + xeval(:, 1).^2 +xeval(:, 2).^2 +(c(theta_star) - xeval(:, 3)).^2)./sqrt(xeval(:, 1).^2+xeval(:, 2).^2));
% Determine value of theta and phi which provide roots in the complex plane
% *Plus or minus imaginary part, so for simplicity we choose positive WLOG
phi0tstar = mod(atan2(xeval(:, 2), xeval(:, 1)),2*pi)  + 1i*log(lambda_PR - sqrt(lambda_PR.^2 - 1));
% Geometric factor
Gp0t=2*sum((gamma_tilde_S(phi0tstar)-xeval).*dgamma_tilde_S(phi0tstar));


%% Approx. theta0star (root in theta direction)
% Determine roots in theta using Newton solver
gamma_tilde_T  = @(th)     axsym_vector(th, s_star, a, c, imap);
% Note derivative factor to account for mapping
dgamma_tilde_T = @(th)     axsym_drdtheta(th, s_star, dadt, dcdt, imap, dfac);

% If target point lies on z = 0, use Newton solver to determine
% root. Otherwise, use analytic formula for spheroids
if xeval(:, 3) ~= 0 & strcmp(shape, 'spheroid')
  %% Analytically determine phi roots for Prolate Spheroid
  theta0_phistar = find_theta_root_spheroid(xeval, a(pi/2), c(0), s_star);
  t0_phistar     = (2/pi)*theta0_phistar - 1;
else
  % We make a small imaginary perturbation of the initial guess to ensure the
  % Newton solver finds the complex root
  t0_phistar = f_newton(gamma_tilde_T,dgamma_tilde_T,xeval,t_star+1i*0.1);
  theta0_phistar = imap(t0_phistar);
end
% Geometric factor
Gt0p=2*sum((gamma_tilde_T(t0_phistar)-xeval).*dgamma_tilde_T(t0_phistar));

end

%% semiroot: Semi-analytical estimate by numeric integration in second direction
% Notation used in paper (maybe with s<->t) % https://doi.org/10.1016/j.camwa.2022.02.001
function ds0fun = semiroot(rvec, drds, drdt)
r2 = norm(rvec)^2;
r_dot_drdt = dot(rvec,drdt);
drdt2 = norm(drdt)^2;
r_dot_drds = dot(rvec,drds);
drdt_dot_drds = dot(drdt,drds);
drds2 = norm(drds)^2;
%     a = r2 + 2*r_dot_drdt*dt + drdt2*dt.^2;
%     b = 2*r_dot_drds + 2*drdt_dot_drds*dt;
%     c = drds2;
%     % Explicit root for first-order linearization
%     ds0fun = -b./(2*c) + 1i*sqrt(a./c - (b./(2*c)).^2);
aa = @(dt) r2 + 2*r_dot_drdt*dt + drdt2*dt.^2;
bb = @(dt) 2*r_dot_drds + 2*drdt_dot_drds*dt;
cc = @(dt) drds2;
% Explicit root for first-order linearization
ds0fun = @(dt) -bb(dt)./(2*cc(dt)) + 1i*sqrt(aa(dt)./cc(dt) - (bb(dt)./(2*cc(dt))).^2);
end

%% axsym_drdphi: Derivative of parameterization in phi
function fsym = axsym_drdphi(t, phi, a, imap)

if nargin > 3
  % Inverse linear mapping
  theta = imap(t);
else
  theta = t;
end
% Derivative of parameterization wrt phi,
% defined as a function of theta and phi
fsym = [-a(theta).*sin(phi),a(theta).*cos(phi),0*theta];

end

%% element_normal: Element normal
function [W, n, normvec] = element_normal(t, phi, a, da, dc, imap, dfac)

if nargin < 6
  % *** Defined such that t = theta
  % Derivative wrt theta, at root
  drdtheta     = axsym_drdtheta(t, phi, da, dc);
  % Derivative wrt phi, at root
  drdphi       = axsym_drdphi(t, phi, a);
else
  % Derivative wrt theta, at root
  drdtheta     = axsym_drdtheta(t, phi, da, dc, imap, dfac);
  % Derivative wrt phi, at root
  drdphi       = axsym_drdphi(t, phi, a, imap);
end

% Normal direction
n            = cross(drdtheta, drdphi);
% Determinant of normal
W            = norm(n);
% Normal vector
normvec      = n/W;

end

%% all_differences: Computes all possible differences r = x - y
function r = all_differences(x, y)
%
% Computes all possible differences r = x - y
%
% Input:
% - x: M×3 matrix of M points
% - y: N×3 matrix of N points
%
% Output:
% - r: N×M×3 array containing all possible differences x - y.
r = zeros(size(y,1), size(x,1), 3); % N×M×3
for d=1:3
  xd = x(:,d).'; % 1×M
  yd = y(:,d); % N×1
  rd = bsxfun(@minus, xd, yd); % N×M
  r(:,:,d) = rd;
end
end

%% custom_scalar_product: Computes the scalar product 
function ab = custom_scalar_product(a, b)
%
% Input:
% - a: N×m×3 array
% - b: N×3 matrix or 1×3 matrix
%
% Output:
% - ab: N×m matrix
%
% Computes the scalar product of a and b by reshaping b to an
% N×1×3 or 1×1×3 array and then doing the scalar product
% along the third dimension.
b = reshape(b, [size(b,1) 1 3]); % N×1×3 or 1×1×3
ab = sum(bsxfun(@times, a, b), 3); % N×m
end

%% gausslaguerre8: Gauss-Laguerre quadrature points and weights
function [xlag,wlag] = gausslaguerre8()

xlag = [
  0.170279632305101
  0.903701776799380
  2.251086629866131
  4.266700170287659
  7.045905402393466
  10.758516010180996
  15.740678641278004
  22.863131736889265
  ]';
wlag = [
  0.369188589341639
  0.418786780814343
  0.175794986637171
  0.033343492261216
  0.002794536235226
  0.000090765087734
  0.000000848574672
  0.000000001048001
  ]';
end

%% trapz_errfunc_deriv: qth derivative of the trapezoidal error function
function knq = trapz_errfunc_deriv(z, n, q)
% knq = gl_errfunc_deriv(z, n, q)
%
% qth derivative of the trapezoidal error function
% returns absolute value

b   = abs(imag(z));
knq = 2*pi*n^q*exp(-n*b);
end

%% gl_errfunc_deriv: qth derivative of the Gauss-Legendre error function
function knq = gl_errfunc_deriv(z, n, q)
% knq = gl_errfunc_deriv(z, n, q)
%
% qth derivative of the Gauss-Legendre error function,
% returns absolute value

% Transform to first quadrant
z = abs(real(z)) + 1i*(imag(z));


knq = 2*pi./(z + sqrt(z.^2-1)).^(2*n+1);
if q ~= 0
  knq = knq .* (-(2*n+1)./sqrt(z.^2-1)).^q;
end

knq = abs(knq);
end


