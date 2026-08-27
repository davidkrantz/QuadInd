function [x, w] = GaussLegendre(N, a, b)
    % GAUSSLEGENDRE Compute Gauss-Legendre quadrature nodes and weights
    %
    %   [x, w] = GaussLegendre(N, a, b) computes N Gauss-Legendre nodes
    %   and weights on the interval [a, b].
    %
    % Inputs:
    %   N - Number of quadrature points
    %   a - Left endpoint of interval
    %   b - Right endpoint of interval
    %
    % Outputs:
    %   x - Column vector of N quadrature nodes
    %   w - Column vector of N quadrature weights
    %
    % Example:
    %   [x, w] = quadind.util.GaussLegendre(10, 0, pi);
    %   integral_approx = sum(sin(x) .* w);  % Approximates integral of sin on [0,pi]
    %
    % Based on code by Greg von Winckel (BSD license, see below)
    
    arguments
        N (1,1) {mustBePositive, mustBeInteger}
        a (1,1) {mustBeReal}
        b (1,1) {mustBeReal}
    end
    
    if a >= b
        error('quadind:GaussLegendre:invalidInterval', ...
            'Interval must satisfy a < b');
    end
    
    N = N - 1;
    N1 = N + 1;
    N2 = N + 2;
    
    xu = linspace(-1, 1, N1)';
    
    % Initial guess
    y = cos((2*(0:N)' + 1) * pi / (2*N + 2)) + (0.27/N1) * sin(pi * xu * N / N2);
    
    % Legendre-Gauss Vandermonde Matrix
    L = zeros(N1, N2);
    
    % Compute zeros of N+1 Legendre polynomial using Newton-Raphson
    y0 = 2;
    
    while max(abs(y - y0)) > eps
        L(:,1) = 1;
        L(:,2) = y;
        
        for k = 2:N1
            L(:,k+1) = ((2*k - 1) * y .* L(:,k) - (k-1) * L(:,k-1)) / k;
        end
        
        Lp = (N2) * (L(:,N1) - y .* L(:,N2)) ./ (1 - y.^2);
        
        y0 = y;
        y = y0 - L(:,N2) ./ Lp;
    end
    
    % Linear map from [-1,1] to [a,b]
    x = (a * (1 - y) + b * (1 + y)) / 2;
    
    % Compute the weights
    w = (b - a) ./ ((1 - y.^2) .* Lp.^2) * (N2/N1)^2;
end

% Copyright (c) 2009, Greg von Winckel
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
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
