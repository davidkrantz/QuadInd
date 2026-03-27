% Analytic expression for determining roots of spheroid
function [theta0,theta0_all,beta] = find_theta_root_spheroid(xvals, Aa, Ca, phi)
if Aa == Ca
  error('A cannot equal C');
end

theta0_all = zeros(size(xvals, 1), 8);
beta   = zeros(size(xvals, 1), 4);

delta = Aa^2 - Ca^2;
tau   = Ca*xvals(:, 3)+1i*Aa*(xvals(:, 1).*cos(phi)+xvals(:, 2).*sin(phi));
d2    = xvals(:, 1).^2+xvals(:, 2).^2+xvals(:, 3).^2+Ca^2;

A=delta/4;
B=conj(tau);
C=-(delta/2+d2);
D=tau;
E=delta/4;

delta0 = C.^2-3*B.*D+12*A.*E;
delta1 = 2*C.^3-9*B.*C.*D+27*(B.^2).*E+27*A.*(D.^2)-72*A.*C.*E;

p = (8*A.*C-3*B.^2)./(8*A.^2);
q = (B.^3-4*A.*B.*C+8*(A.^2).*D)./(8*A.^3);

Q = ((delta1+sqrt(delta1.^2-4*delta0.^3))/2).^(1/3);
S = sqrt(-(2/3)*p+(Q+delta0./Q)/(3*A))/2;

beta(:, 1) = -B./(4*A)-S+sqrt(-4*S.^2-2*p+q./S)/2;
beta(:, 2) = -B./(4*A)-S-sqrt(-4*S.^2-2*p+q./S)/2;
beta(:, 3) = -B./(4*A)+S+sqrt(-4*S.^2-2*p-q./S)/2;
beta(:, 4) = -B./(4*A)+S-sqrt(-4*S.^2-2*p-q./S)/2;

theta0_all(:, 1)=angle(beta(:, 1))+1i*log(abs(beta(:, 1)));
theta0_all(:, 2)=angle(beta(:, 1))-1i*log(abs(beta(:, 1)));
theta0_all(:, 3)=angle(beta(:, 2))+1i*log(abs(beta(:, 2)));
theta0_all(:, 4)=angle(beta(:, 2))-1i*log(abs(beta(:, 2)));
theta0_all(:, 5)=angle(beta(:, 3))+1i*log(abs(beta(:, 3)));
theta0_all(:, 6)=angle(beta(:, 3))-1i*log(abs(beta(:, 3)));
theta0_all(:, 7)=angle(beta(:, 4))+1i*log(abs(beta(:, 4)));
theta0_all(:, 8)=angle(beta(:, 4))-1i*log(abs(beta(:, 4)));

%     [~,idx] = min(abs(imag(theta0_all)));
[~,idx] = min(abs(imag(theta0_all)), [], 2);
theta0 = zeros(size(theta0_all, 1), 1);
for ii = 1:size(theta0_all, 1)
  theta0(ii, 1) = theta0_all(ii, idx(ii));
end


% gives same as above, but sometimes conj difference, so no practical
% difference
%[theta0new,theta0_allnew,betanew] = find_theta_root_spheroid_manuscript_derivation(xvals, Aa, Ca, phi);
%transpose(theta0_allnew)
end
