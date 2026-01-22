% Direct Quadrature of stresslet kernel
function u = stresslet(x, y, n, q)
% u = stresslet(x, y, n(y), q(y))
%
% Compute stresslet potential
% u_j = q_i r_i r_j r_k n_k / r^5

if size(y, 1) == 1
  % Vectorize in x
  r = bsxfun(@minus, x, y);
  r5i = sum(r.^2, 2).^(-5/2);
  qr = r*q';
  rn = r*n';
  u = -6*bsxfun(@times, r, qr.*rn.*r5i);
else
  % Iterate
  u = zeros(size(x));
  for i=1:size(x,1)
    % Vectorize in y
    r = bsxfun(@minus, x(i,:), y);
    r5i = sum(r.^2, 2).^(-5/2);
    qr = sum(q.*r, 2);
    rn = sum(r.*n, 2);
    u(i,:) = -6*sum(bsxfun(@times, r, qr.*rn.*r5i), 1);
  end
end
end
