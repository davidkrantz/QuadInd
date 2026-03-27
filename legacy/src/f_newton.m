% Simple Newton solver
function [t0,Gt0]=f_newton(gamma, dgamma, xt, t_init)

% Functional
R2= @(t) sum((gamma(t)-xt).^2, 2);
% Derivative of functional
R2_t=@(t) 2*sum((gamma(t)-xt).*(dgamma(t)), 2);

% Initialize solver
add=1;
t0=t_init;
k=0; kmax=100;

% Loop through to determine root
while (abs(add)>1e-12 && k<kmax)
  k   = k+1;
  add = R2(t0)/R2_t(t0);
  t0  = t0-add;
end

% If the solver fails, try again with more steps and smaller step size to
% avoid overshooting
if k==kmax || isnan(real(t0)) || isnan(imag(t0))
  disp('not converging! Trying again, smaller stepsize')

  % Initialize again
  add=1;
  t0=t_init;
  k=0; kmax=2000;

  % Loop through to determine root
  while (abs(add)>1e-12 && k<kmax)
    k   = k+1;
    add = R2(t0)/R2_t(t0);
    t0  = t0-0.1*add; % Now taking step size of 0.1
  end
end

% Printed if max iterations is reached after a second time
if k==kmax
  disp('not converging!')
end

% Geometric factor
Gt0=2*dot(gamma(t0)-xt,dgamma(t0));
end
