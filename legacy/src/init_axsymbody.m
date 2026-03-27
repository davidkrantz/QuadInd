% Function to setup axis-symmetric particle body
function body = init_axsymbody(geometry, nph, nth)

% Select axisymmetric geometry
[at, ct, dadth, dcdth] = axsym_parameterization(geometry);

body = struct();

% Uniformly spaced points in periodic direction
phi = linspace(0,2*pi,nph+1);
phi = phi(1:nph);
% Quadrature weights for trapezoidal rule
wphi = (phi(2)-phi(1))*ones(size(phi));

% Gauss-Legendre directly on [0, pi]
[s,w] = lgwt(nth, 0, pi);
theta = s;
wtheta = w;

% x, y, z coordinates of body
X = at(theta)*cos(phi);
Y = at(theta)*sin(phi);
Z = ct(theta)*ones(size(phi));

% Jacobian for body
[ J, N1, N2, N3 ] = axsym_jacobian( at, dadth, dcdth, phi, theta );

% Quadrature weights, accounting for Jacobian
W = J.*(wtheta * wphi);

% Store values of the body grid
body.x = [X(:) Y(:) Z(:)];
body.n = [N1(:) N2(:) N3(:)];
body.w = W(:);
body.phi = phi;
body.theta = theta;
body.nph = nph;
body.nth = nth;

% Store parameterization
body.geometry = geometry;
% Use exact semi-axis values, not max over GL nodes (which gives incorrect values)
switch lower(geometry.shape)
  case 'spheroid'
    % Use defaults if not specified (must match axsym_parameterization defaults)
    if ~isfield(geometry, 'a'); geometry.a = 0.05; end
    if ~isfield(geometry, 'c'); geometry.c = 0.1;  end
    body.a = geometry.a;  % Exact semi-axis in xy-plane
    body.c = geometry.c;  % Exact semi-axis in z
  otherwise
    body.a = max(at(theta));  % For non-spheroid, use max of at(theta)
    body.c = max(ct(theta));
end
body.at = at;
body.ct = ct;
body.dadt = dadth;
body.dcdt = dcdth;

end

%% axsym_parameterization: particle parameterization, based on body selected
function [at,ct,dadt,dcdt] = axsym_parameterization(geometry)
switch lower(geometry.shape)
  case 'spheroid'
    % Default
    if ~isfield(geometry, 'a'); geometry.a = 0.05; end
    if ~isfield(geometry, 'c'); geometry.c = 0.1;  end

    a = @(theta) geometry.a*ones(size(theta));
    c = @(theta) geometry.c*ones(size(theta));
    dadth = @(theta) zeros(size(theta));
    dcdth = @(theta) zeros(size(theta));
  case 'peanut'
    a = @(theta) 2.2 + cos(2*theta);
    c = @(theta) 2.2 + cos(2*theta);
    dadth = @(theta) -2*sin(2*theta);
    dcdth = @(theta) -2*sin(2*theta);
  % case 'star'
  %   a = @(theta) 1 + 0.3*cos(8*theta);
  %   c = @(theta) 1 + 0.3*cos(8*theta);
  %   dadth = @(theta) -0.3*8*sin(8*theta);
  %   dcdth = @(theta) -0.3*8*sin(8*theta);
  % case 'smootherstar'
  %   a = @(theta) 1 + 0.2*cos(5*theta);
  %   c = @(theta) 1 + 0.2*cos(5*theta);
  %   dadth = @(theta) -0.2*5*sin(5*theta);
  %   dcdth = @(theta) -0.2*5*sin(5*theta);
  otherwise
    % warning("GET_PARAMETERIZATION: geometry must be 'spheroid', 'peanut', 'star' or 'smootherstar'");
    warning("GET_PARAMETERIZATION: geometry must be 'spheroid' or 'peanut' ");
end
% Store parameterization and derivatives in theta and phi
at = @(theta) a(theta).*sin(theta);  % radius in x and y
ct = @(theta) c(theta).*cos(theta);  % radius in z
dadt = @(theta) dadth(theta).*sin(theta)+a(theta).*cos(theta);  % derivative wrt theta in x and y
dcdt = @(theta) dcdth(theta).*cos(theta)-c(theta).*sin(theta);  % derivative wrt theta in z
end

%% axsym_jacobian: Jacobian and normals function for particle
function [J, NX, NY, NZ] = axsym_jacobian( a, da, dc, phi, theta )
% a           --- distance to origin in x and y, dependent on theta
% c           --- distance to origin in z, dependent on theta
% da          --- derivative wrt theta
% dc          --- derivative wrt theta
% phi is [0,2pi]  --- polar angle
% theta is [0,pi] --- azimuthal angle

% General axisymmetric rigid body
at   = a(theta);
dadt = da(theta);
dcdt = dc(theta);

% Compute Jacobian determinant components
J1 = -dcdt.*at*cos(phi);
J2 = -dcdt.*at*sin(phi);
J3 = dadt.*at*ones(size(phi));

J = sqrt(J1.^2 + J2.^2 + J3.^2);

NX = J1./J;
NY = J2./J;
NZ = J3./J;

% Jacobian can zero out at poles, causing NaN normals
idx_nan = (isnan(NX) | isnan(NY) | isnan(NZ));
if any(idx_nan)
    idxN = find(idx_nan & abs(theta) < 1e-10);
    idxS = find(idx_nan & abs(theta-pi) < 1e-10);
    NX(idxN) = 0;
    NY(idxN) = 0;
    NZ(idxN) = 1;
    NX(idxS) = 0;
    NY(idxS) = 0;
    NZ(idxS) = -1;
end

if any(isnan(NX(:)) | isnan(NY(:)) | isnan(NZ(:)))
    error('get_jacobian failed (NaN in normals)')
end

end

%% lgwt: Gauss-Legendre quadrature weights
function [x,w]=lgwt(N,a,b)
% lgwt.m
%
% [x,w]=lgwt(N,a,b)
%
% This script is for computing definite integrals using Legendre-Gauss
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
%
% This code is covered by the BSD license (see bottom of file).

N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps


  L(:,1)=1;
  Lp(:,1)=0;

  L(:,2)=y;
  Lp(:,2)=1;

  for k=2:N1
    L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
  end

  Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);

  y0=y;
  y=y0-L(:,N2)./Lp;

end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

% Copyright (c) 2009, Greg von Winckel
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
end
