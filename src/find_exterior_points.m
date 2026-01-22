% Determine which points lie outside of particle body
function [x_exterior, mask] = find_exterior_points(body, x)

% Spheroid have analytic expression, otherwise iteratively determine
switch lower(body.geometry.shape)
  case 'spheroid'
    mask = ( (x(:,1).^2+x(:,2).^2)/body.a^2 + x(:,3).^2/body.c^2 > 1 );
  otherwise
    % Determine distance from origin
    distx = sum(x.^2, 2);
    theta = zeros(size(x, 1), 1);
    % Since axissymmetric, can determine phi value
    phival = mod(atan2(x(:, 2), x(:, 1)), 2*pi);

    % Loop through and determine closest value that lies on the body
    % *** Can be improved beyond fminbnd
    for i = 1:size(x, 1)
      % fun = @(theta)  (x(i, 3) - body.ct(theta)*cos(theta) );
      gamma  = @(t) axsym_vector(t, phival(i), body.at, body.ct);
      R = @(t) sum((gamma(t)-x(i, :)).^2, 2);
      theta(i) = fminbnd(R, 0, pi);
    end

    % Determine corresponding onsurface value, and determine if distance is
    % greater than point distance
    onsurf = axsym_vector(theta, phival, body.at, body.ct);
    mask = sum(onsurf.^2, 2) < distx;

end

% Output evaluation points outside of particle
x_exterior = x(mask, :);

end
