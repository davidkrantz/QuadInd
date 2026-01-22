% Function for upsample via fixed factor
function new = upsample(body, k)
% Upsample grid by a factor k in each direction, returns a new grid
new = init_axsymbody(body.geometry, k*body.nph, k*body.nth);

% Upsampler for density
fi = get_upsampler(body, new);
new.interpolate = @(q) [fi(q(:,1)),fi(q(:,2)),fi(q(:,3))];

end

%% get_upsampler: interpolation handle to upsampled grid
function hdl = get_upsampler(grid1, grid2)

[T, B] = upsampling_pre(grid1, grid2);
T = T.';

m = numel(grid1.phi);
n = numel(grid1.theta);

hdl = @(y) reshape(B*reshape(y, n, m)*T, [], 1);

end

%% upsampling_pre: interpolation weights in both directions
function [T, B] = upsampling_pre(grid1, grid2)

T = trig_interp_matrix(grid1.phi(:), grid2.phi(:));

x1 = grid1.theta(:);
x2 = grid2.theta(:);
bclag_pre = bclag_interp_weights(x1);
B = bclag_interp_matrix(x1, bclag_pre, x2);
end

%% trig_interp_matrix: interpolation in periodic direction
function T = trig_interp_matrix(x, xi)
% T = trig_interp_matrix(x, xi)
%
% fi = T*f
% x must be equispaced
% x and xi assumed columns

assert(size(x,2)==1)
assert(size(xi,2)==1)

n = numel(x);
N = numel(xi);

% DFT matrix (precompute possible, but cost is negligible)
D = fft(eye(n)); % Identical to dftmtx(n)
p = ceil(n/2);
D = [D(p+1:end,:);D(1:p,:)]; % fftshift

% Non-uniform inverse DFT matrix
I = complex(zeros(N,n));
q = floor(n/2);
z = exp(1i*xi);

% Main loop, this is what costs
w = z.^-q;
I(:,1) = w;
for l=2:n
  w = w.*z;
  I(:,l) = w;
end
I = I/n;

% % Slower, but more accurate than looped version
% for l=1:n
%     k = l-1-q;
%     I(:,l) = z.^k;
% end
% I = I/n;

% Want real(I*D)
T = real(I)*real(D) - imag(I)*imag(D);
end

%% trig_interp_matrix: interpolation in GL direction
function w = bclag_interp_weights(x)

% Barycentric Lagrange Interpolation.
% Berrut, J.-P., & Trefethen, L. N. (2004).
% SIAM Review, 46(3), 501–517. doi:10.1137/S36144502417716

assert(size(x,2)==1);
n = numel(x);

w = zeros(size(x));
for j=1:n
  w(j) = 1/prod(x(j)-x(1:n~=j));
end
end

function B = bclag_interp_matrix(x, w, xi)

% Barycentric Lagrange Interpolation.
% Berrut, J.-P., & Trefethen, L. N. (2004).
% SIAM Review, 46(3), 501–517. doi:10.1137/S0036144502417715

assert(size(x,2)==1)
assert(size(xi,2)==1)
assert(all(size(x)==size(w)));

n = numel(x);
N = numel(xi);

B = zeros(N,n);
[denom exact] = deal(zeros(size(xi)));

for j=1:n
  xdiff = xi-x(j);
  temp = w(j)./xdiff;
  B(:,j) = temp;
  denom = denom + temp;
  exact(xdiff==0) = j;
end

B = bsxfun(@rdivide,B,denom);
jj = find(exact);
B(jj,:) = 0;
B(jj + N*(exact(jj)-1)) = 1;
end