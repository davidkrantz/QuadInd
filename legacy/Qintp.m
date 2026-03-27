% Modify density at roots, using local linear interpolation
function [qp, qt, phi0tstar, intpind_p] = Qintp(body, xt, a, c, dens_q, ind_R, thetaroot_int)

%% Find closest quadrature node (on body) to evaluation point
[i1,i2]    = ind2sub([body.nth,body.nph], ind_R);
theta_star = body.theta(i1);        % theta value of point on body  [pi, 0]
s_star     = body.phi(i2);          % phi value of point on body    [0, 2pi]
s_star     = s_star(:);

%% Analytically determine phi roots
% % x, y, z coordinates
xcoord = xt(:, 1);
ycoord = xt(:, 2);
zcoord = xt(:, 3);

% Expressions for finding roots analytically
lambda_pr = (1./(2*a)).*((a^2 + xcoord.^2 +ycoord.^2 +(c - zcoord).^2)./sqrt(xcoord.^2+ycoord.^2));
% Determine value of theta and phi which provide roots in the complex plane
% *Plus or minus imaginary part, so for simplicity we choose positive WLOG
phi0tstar = mod(atan2(ycoord, xcoord), 2*pi)  + 1i*log(lambda_pr - sqrt(lambda_pr.^2 - 1));

% Density evaluated at root in phi with fixed theta
[f_s, intpind_p] = localintp_est_phi0(body, dens_q, phi0tstar, theta_star);
% Error estimate from phi direction contribution (trapzoidal integration)
qp = abs(f_s);

%% Determine theta roots
if nargin < 7
  % If target point lies on z = 0, use Newton solver to determine
  % root. Otherwise, use analytic formula for spheroids
  if xt(:, 3) ~= 0 & strcmp(body.geometry.shape, 'spheroid')
    % Analytically determine phi roots for Prolate Spheroid
    theta0phistar = find_theta_root_spheroid(xt, body.at(pi/2), body.ct(0), s_star);
  else
    % Use Newton solver
    t_star     = (2/pi)*theta_star(:) - 1; % Linear map in t, Full Sphere  [-1, 1]
    imap       = @(t) (pi/2)*(t + 1);   % inverse linear mapping, Full Sphere
    dfac       = (pi/2);                % factor from derivative chain rule
    t0_phistar = zeros(size(xt, 1), 1);
    for i = 1:size(xt, 1)
      gamma_tilde_T  = @(th)     axsym_vector(th, s_star(i), body.at, body.ct, imap);
      % Note derivative factor to account for mapping
      dgamma_tilde_T = @(th)     axsym_drdtheta(th, s_star(i), body.dadt, body.dcdt, imap, dfac);
      % We make a small imaginary perturbation of the initial guess to ensure the
      % Newton solver finds the complex root
      t0_phistar(i) = f_newton(gamma_tilde_T,dgamma_tilde_T,xt(i, :),t_star(i)+1i*0.1);
    end
    theta0phistar = imap(t0_phistar);
  end
else
  % Use interpolated root
  theta0phistar = thetaroot_int;
end

% Density evaluated at root in theta with fixed phi
f_t = localintp_th0(body, dens_q, theta0phistar, theta_star, s_star);
% Error estimate from theta direction contribution (Gauss-Legendre integration)
qt = abs(f_t);


end % End of errestintp


%% localintp_est_phi0: Local interpolation at a complex root in phi
function [q0, intpind_p] = localintp_est_phi0(body, fdens, phi00, thetavals)
%                 f0_loc = fsigma(thmap_loc(tGL),phi0);
% Number of interpolation points
n_intp = 2;

% Theta of base grid
tmpth  = body.theta;

% Loop through each component
q0       = zeros(size(thetavals, 1), 3);
intpind_p = zeros(size(phi00, 1), n_intp);
for ii = 1:length(thetavals)
  % Must do interpolation one target point at a time
  theta_star = thetavals(ii);
  phi_star = real(phi00(ii));

  % Distance to point wrt to theta and phi, sorted
  tmpphi = body.phi;
  [~, iDP]      = sort(abs(tmpphi - phi_star));
  [~, local_it] = min(abs(tmpth - theta_star));
  % Distance to point sorted
  if iDP(1) < n_intp % Create local interpolant with 2 points, accounting for periodicity
    tmpphi(end-(n_intp-1):end) = abs(tmpphi(end-(n_intp-1):end) - 2*pi);
    [~, iDP] = sort(abs(tmpphi - phi_star));
  elseif iDP(1) > length(tmpphi) - n_intp
    tmpphi(1:n_intp) = abs(tmpphi(1:n_intp) + 2*pi);
    [~, iDP] = sort(abs(tmpphi - phi_star));
  end
  % Sorted indices
  local_ip = sort(iDP(1:n_intp));
  % Phi values for local interpolation
  p_pts    = body.phi(local_ip);

  % Store local interpolation index for later use in S3Q
  intpind_p(ii, :) = local_ip;

  % Local interpolation
  q_ppts = fdens(local_it, local_ip, :);
  q0(ii,1) = (q_ppts(1,1,1).*(p_pts(2) - phi00(ii)) + q_ppts(1,2,1).*(phi00(ii) - p_pts(1)))/(p_pts(2) - p_pts(1));
  q0(ii,2) = (q_ppts(1,1,2).*(p_pts(2) - phi00(ii)) + q_ppts(1,2,2).*(phi00(ii) - p_pts(1)))/(p_pts(2) - p_pts(1));
  q0(ii,3) = (q_ppts(1,1,3).*(p_pts(2) - phi00(ii)) + q_ppts(1,2,3).*(phi00(ii) - p_pts(1)))/(p_pts(2) - p_pts(1));

end

end % End of localintp_phi0

%% localintp_th0: Local interpolation at a complex root in phi
function q0 = localintp_th0(body, fdens, theta00, thetavals, phivals)

% Loop through each component
q0 = zeros(length(phivals), 3);
qt = zeros(1,3);
for ii = 1:length(phivals)
  phi_star   = phivals(ii);
  theta_star = thetavals(ii);

  % Number of interpolation points
  n_intp = 2;
  % Distance to point wrt to theta and phi, sorted
  tmpphi = body.phi;
  tmpth  = body.theta;
  [~, local_ip] = min(abs(tmpphi - phi_star));
  [~, local_it] = mink(abs(tmpth - theta_star),n_intp);
  % Theta values for local interpolation
  t_pts    = tmpth(local_it);

  % Local interpolation
  q_tpts = fdens(local_it, local_ip, :);
  % Use 2-point linear interpolation to evaluate at exact theta root
  qt(1) = (q_tpts(1,1,1).*(t_pts(2) - theta00(ii)) + q_tpts(2,1,1).*(theta00(ii) - t_pts(1)))/(t_pts(2) - t_pts(1));
  qt(2) = (q_tpts(1,1,2).*(t_pts(2) - theta00(ii)) + q_tpts(2,1,2).*(theta00(ii) - t_pts(1)))/(t_pts(2) - t_pts(1));
  qt(3) = (q_tpts(1,1,3).*(t_pts(2) - theta00(ii)) + q_tpts(2,1,3).*(theta00(ii) - t_pts(1)))/(t_pts(2) - t_pts(1));
  qt(isnan(qt)) = 0;

  % Modified estimate
  q0(ii, :) = qt;
end

end % End of localintp_th0