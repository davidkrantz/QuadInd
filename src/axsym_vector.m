%% axsym_vector: Parameterization of particle body
function fsym = axsym_vector(t, phi, a, c, imap)

if nargin > 4
  % Inverse linear mapping
  theta = imap(t);
else
  theta = t;
end

% Spheroid parameterization, defined as a function of theta and phi
fsym = [a(theta).*cos(phi),a(theta).*sin(phi),c(theta)];

end
