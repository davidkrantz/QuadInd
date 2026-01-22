%% axsym_drdtheta: Derivative of parameterization in theta
function fsym = axsym_drdtheta(t, phi, da, dc, imap, dfac)

if nargin > 4
  % Inverse linear mapping
  theta = imap(t);
else
  theta = t;
  dfac  = 1;
end
% Derivative of spheroid parameterization wrt theta, compensating for map
% defined as a function of theta and phi
fsym = [dfac*da(theta).*cos(phi),dfac*da(theta).*sin(phi),dfac*dc(theta)];

end
