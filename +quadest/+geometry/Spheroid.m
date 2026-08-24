classdef Spheroid < quadest.geometry.AxsymGeometry
    % SPHEROID Prolate or oblate spheroid geometry
    %
    % A spheroid is an ellipsoid of revolution, parameterized as:
    %   x = a * sin(theta) * cos(phi)
    %   y = a * sin(theta) * sin(phi)
    %   z = c * cos(theta)
    %
    % where a is the semi-axis in the xy-plane and c is the semi-axis in z.
    %   - Prolate spheroid: c > a (elongated along z)
    %   - Oblate spheroid: a > c (flattened along z)
    %
    % Usage:
    %   geom = quadest.geometry.Spheroid('a', 0.05, 'c', 0.1);  % prolate
    %   geom = quadest.geometry.Spheroid();  % defaults: a=0.05, c=0.1
    %
    % The spheroid has analytic formulas for:
    %   - Interior/exterior test
    %   - Complex theta root (quartic solution)
    %
    % See also: quadest.geometry.AxsymGeometry, quadest.geometry.Peanut
    
    properties (SetAccess = private)
        a   % Semi-axis in xy-plane
        c   % Semi-axis in z direction
    end
    
    methods
        function obj = Spheroid(varargin)
            % SPHEROID Construct a spheroid geometry
            %
            %   geom = Spheroid() creates spheroid with default a=0.05, c=0.1
            %   geom = Spheroid('a', val, 'c', val) specifies semi-axes
            
            % Parse inputs
            p = inputParser;
            addParameter(p, 'a', 0.05, @(x) isscalar(x) && x > 0);
            addParameter(p, 'c', 0.1, @(x) isscalar(x) && x > 0);
            parse(p, varargin{:});
            
            obj.a = p.Results.a;
            obj.c = p.Results.c;
        end
        
        function val = at(obj, theta)
            % AT Radial distance in xy-plane: a * sin(theta)
            val = obj.a * sin(theta);
        end
        
        function val = ct(obj, theta)
            % CT Z-coordinate: c * cos(theta)
            val = obj.c * cos(theta);
        end
        
        function val = dadt(obj, theta)
            % DADT Derivative of at w.r.t. theta: a * cos(theta)
            val = obj.a * cos(theta);
        end
        
        function val = dcdt(obj, theta)
            % DCDT Derivative of ct w.r.t. theta: -c * sin(theta)
            val = -obj.c * sin(theta);
        end
        
        function val = maxRadius(obj)
            % MAXRADIUS Maximum radial extent in xy-plane
            val = obj.a;
        end
        
        function val = maxHeight(obj)
            % MAXHEIGHT Maximum z extent
            val = obj.c;
        end
        
        function mask = isExterior(obj, points)
            % ISEXTERIOR Test if points are outside the spheroid
            %
            %   mask = isExterior(obj, points) returns logical array
            %   where true indicates the point is exterior.
            %
            % Uses the ellipsoid equation: x²/a² + y²/a² + z²/c² > 1
            
            arguments
                obj
                points (:,3) {mustBeNumeric}
            end
            
            x = points(:, 1);
            y = points(:, 2);
            z = points(:, 3);
            
            mask = (x.^2 + y.^2) / obj.a^2 + z.^2 / obj.c^2 > 1;
        end
        
        function theta0 = findThetaRoot(obj, targets, phi)
            % FINDTHETAROOT Find complex theta root using quartic formula
            %
            %   theta0 = findThetaRoot(obj, targets, phi) computes the complex
            %   theta value where the kernel becomes singular.
            %
            % This uses an analytic quartic solution specific to spheroids.
            %
            % Inputs:
            %   targets - [M×3] target points
            %   phi     - [M×1] phi values (typically from closest grid point)
            %
            % Output:
            %   theta0 - [M×1] complex theta roots (smallest |imag| selected)
            
            arguments
                obj
                targets (:,3) {mustBeNumeric}
                phi (:,1) {mustBeNumeric}
            end
            
            delta = obj.a^2 - obj.c^2;
            if abs(delta) <= sqrt(eps) * max(obj.a^2, obj.c^2)
                % The quartic loses its leading coefficient in the spherical
                % limit, use instead Newton solver
                theta0 = findThetaRoot@quadest.geometry.AxsymGeometry( ...
                    obj, targets, phi);
                return;
            end
            
            M = size(targets, 1);
            Aa = obj.a;
            Ca = obj.c;
            
            % Quartic coefficients for beta=exp(i*theta)
            tau = Ca * targets(:,3) + 1i * Aa * (targets(:,1) .* cos(phi) + targets(:,2) .* sin(phi));
            d2 = targets(:,1).^2 + targets(:,2).^2 + targets(:,3).^2 + Ca^2;
            
            A = delta / 4;
            B = conj(tau);
            C = -(delta / 2 + d2);
            D = tau;
            E = delta / 4;
            
            % Solve quartic using Ferrari's method
            delta0 = C.^2 - 3*B.*D + 12*A.*E;
            delta1 = 2*C.^3 - 9*B.*C.*D + 27*(B.^2).*E + 27*A.*(D.^2) - 72*A.*C.*E;
            
            p = (8*A.*C - 3*B.^2) ./ (8*A.^2);
            q = (B.^3 - 4*A.*B.*C + 8*(A.^2).*D) ./ (8*A.^3);
            
            Q = ((delta1 + sqrt(delta1.^2 - 4*delta0.^3)) / 2).^(1/3);
            S = sqrt(-(2/3)*p + (Q + delta0./Q) / (3*A)) / 2;
            
            % Four roots (in beta = exp(i*theta))
            beta = zeros(M, 4);
            beta(:,1) = -B./(4*A) - S + sqrt(-4*S.^2 - 2*p + q./S) / 2;
            beta(:,2) = -B./(4*A) - S - sqrt(-4*S.^2 - 2*p + q./S) / 2;
            beta(:,3) = -B./(4*A) + S + sqrt(-4*S.^2 - 2*p - q./S) / 2;
            beta(:,4) = -B./(4*A) + S - sqrt(-4*S.^2 - 2*p - q./S) / 2;
            
            % Convert to theta via beta = exp(i*theta)
            % theta = angle(beta) + i*log|beta| or angle(beta) - i*log|beta|
            theta_all = zeros(M, 8);
            theta_all(:,1) = angle(beta(:,1)) + 1i*log(abs(beta(:,1)));
            theta_all(:,2) = angle(beta(:,1)) - 1i*log(abs(beta(:,1)));
            theta_all(:,3) = angle(beta(:,2)) + 1i*log(abs(beta(:,2)));
            theta_all(:,4) = angle(beta(:,2)) - 1i*log(abs(beta(:,2)));
            theta_all(:,5) = angle(beta(:,3)) + 1i*log(abs(beta(:,3)));
            theta_all(:,6) = angle(beta(:,3)) - 1i*log(abs(beta(:,3)));
            theta_all(:,7) = angle(beta(:,4)) + 1i*log(abs(beta(:,4)));
            theta_all(:,8) = angle(beta(:,4)) - 1i*log(abs(beta(:,4)));
            
            % Select root with smallest imaginary part
            [~, idx] = min(abs(imag(theta_all)), [], 2);
            theta0 = zeros(M, 1);
            for ii = 1:M
                theta0(ii) = theta_all(ii, idx(ii));
            end

            % Optionally use MATLAB's polynomial root solver instead of
            % explicit Ferrari formula (may be more robust). Convert each
            % beta root with theta=-i*log(beta), include adjacent 2*pi
            % representatives, and choose the root closest to [0,pi].
            % theta0 = zeros(M, 1);
            % for ii = 1:M
            %     beta = roots([A, B(ii), C(ii), D(ii), E]);
            %     thetaPrincipal = -1i * log(beta);
            %     thetaCandidates = [thetaPrincipal - 2*pi; ...
            %                        thetaPrincipal; thetaPrincipal + 2*pi];
            %     re = real(thetaCandidates);
            %     intervalOffset = max(0, -re) + max(0, re - pi);
            %     distance = hypot(intervalOffset, imag(thetaCandidates));

            %     residual = abs(sum((obj.evaluate(thetaCandidates, ...
            %         phi(ii) * ones(size(thetaCandidates))) - targets(ii,:)).^2, 2));
            %     scale = max(sum(abs(obj.evaluate(thetaCandidates, ...
            %         phi(ii) * ones(size(thetaCandidates)))).^2, 2), obj.a^2 + obj.c^2);
            %     normalizedResidual = residual ./ scale;
            %     valid = isfinite(thetaCandidates) & normalizedResidual < 1e-8;
            %     distance(~valid) = Inf;
            %     [bestDistance, idx] = min(distance);
            %     if ~isfinite(bestDistance)
            %         theta0(ii) = findThetaRoot@quadest.geometry.AxsymGeometry( ...
            %             obj, targets(ii,:), phi(ii));
            %     else
            %         theta0(ii) = thetaCandidates(idx);
            %     end
            % end
        end
        
        function disp(obj)
            % DISP Display spheroid parameters
            if obj.c > obj.a
                type = 'prolate';
            elseif obj.a > obj.c
                type = 'oblate';
            else
                type = 'sphere';
            end
            fprintf('Spheroid (%s): a = %.4f, c = %.4f\n', type, obj.a, obj.c);
        end
    end
end
