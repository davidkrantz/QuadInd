classdef RootFinder
    % ROOTFINDER Complex root finding for error estimation
    %
    % Provides static methods for computing complex roots in phi and theta
    % directions, which are needed for error estimation.
    %
    % See also: quadest.errorest.ErrorEstimator
    
    methods (Static)
        function [phi0, Gphi] = computePhiRoot(geometry, targets, theta_star)
            % COMPUTEPHIROOT Compute complex phi roots analytically
            %
            %   [phi0, Gphi] = computePhiRoot(geometry, targets, theta_star)
            %
            % The phi root has an analytic formula for all axisymmetric geometries.
            %
            % Inputs:
            %   geometry    - AxsymGeometry object
            %   targets     - [M×3] target points
            %   theta_star  - [M×1] theta values (typically closest grid point)
            %
            % Outputs:
            %   phi0 - [M×1] complex phi roots
            %   Gphi - [M×1] geometric factors for error estimation
            
            x = targets(:, 1);
            y = targets(:, 2);
            z = targets(:, 3);
            
            % Use the parameterization values at theta_star for the phi root formula
            % NOT the constant semi-axis coefficients!
            % For spheroid: at_val = a*sin(theta_star), ct_val = c*cos(theta_star)
            at_val = geometry.at(theta_star);
            ct_val = geometry.ct(theta_star);
            
            % Analytic formula for phi root
            rxy = sqrt(x.^2 + y.^2);
            lambda = (1 ./ (2 * at_val)) .* ...
                ((at_val.^2 + x.^2 + y.^2 + (ct_val - z).^2) ./ rxy);
            
            phi0 = mod(atan2(y, x), 2*pi) + 1i * log(lambda - sqrt(lambda.^2 - 1));
            
            % Compute geometric factor G = 2 * dot(gamma(phi0) - target, dgamma/dphi(phi0))
            % Note: MATLAB's dot() conjugates the first argument for complex vectors,
            % which matches the legacy implementation.
            if nargout > 1
                gamma_phi0 = geometry.evaluate(theta_star, phi0);
                dgamma_phi0 = geometry.drdphi(theta_star, phi0);
                
                r = gamma_phi0 - targets;
                % Use dot() for each row to match legacy behavior with complex vectors
                Gphi = zeros(size(targets, 1), 1);
                for ii = 1:size(targets, 1)
                    Gphi(ii) = 2 * dot(r(ii,:), dgamma_phi0(ii,:));
                end
            end
        end
        
        function [theta0, Gtheta] = computeThetaRoot(geometry, targets, phi_star, theta_init)
            % COMPUTETHETAROOT Compute complex theta roots
            %
            %   [theta0, Gtheta] = computeThetaRoot(geometry, targets, phi_star, theta_init)
            %
            % For spheroids, uses analytic quartic formula.
            % For other geometries, uses Newton iteration.
            %
            % Inputs:
            %   geometry   - AxsymGeometry object
            %   targets    - [M×3] target points
            %   phi_star   - [M×1] phi values (typically closest grid point)
            %   theta_init - [M×1] initial theta guesses (closest grid point)
            %
            % Outputs:
            %   theta0 - [M×1] complex theta roots
            %   Gtheta - [M×1] geometric factors
            
            % Use geometry's built-in root finder
            theta0 = geometry.findThetaRoot(targets, phi_star);
            
            % Compute geometric factor
            % Note: MATLAB's dot() conjugates the first argument for complex vectors,
            % which matches the legacy implementation.
            if nargout > 1
                % Use linear map t = 2*theta/pi - 1 for geometric factor
                imap = @(t) (pi/2) * (t + 1);
                dfac = pi/2;
                
                t0 = 2*theta0/pi - 1;
                
                gamma_t0 = geometry.evaluate(imap(t0), phi_star);
                dgamma_t0 = dfac * geometry.drdtheta(imap(t0), phi_star);
                
                r = gamma_t0 - targets;
                % Use dot() for each row to match legacy behavior with complex vectors
                Gtheta = zeros(size(targets, 1), 1);
                for ii = 1:size(targets, 1)
                    Gtheta(ii) = 2 * dot(r(ii,:), dgamma_t0(ii,:));
                end
            end
        end
        
        function ds0fun = semiRoot(rvec, drds, drdt)
            % SEMIROOT Semi-analytical root for perpendicular direction integration
            %
            %   ds0fun = semiRoot(rvec, drds, drdt)
            %
            % Returns a function handle ds0fun(dt) that gives the approximate
            % root in the ds direction as a function of perturbation dt.
            %
            % Reference: https://doi.org/10.1016/j.camwa.2022.02.001
            
            r2 = sum(rvec.^2);
            r_dot_drdt = dot(rvec, drdt);
            drdt2 = sum(drdt.^2);
            r_dot_drds = dot(rvec, drds);
            drdt_dot_drds = dot(drdt, drds);
            drds2 = sum(drds.^2);
            
            % Coefficients of quadratic approximation
            aa = @(dt) r2 + 2*r_dot_drdt*dt + drdt2*dt.^2;
            bb = @(dt) 2*r_dot_drds + 2*drdt_dot_drds*dt;
            cc = @(dt) drds2;
            
            % Explicit root for first-order linearization
            ds0fun = @(dt) -bb(dt)./(2*cc(dt)) + 1i*sqrt(aa(dt)./cc(dt) - (bb(dt)./(2*cc(dt))).^2);
        end
        
        function [t0, converged] = newtonSolve(gamma, dgamma, target, t_init, config)
            % NEWTONSOLVE Newton iteration for complex root finding
            %
            %   [t0, converged] = newtonSolve(gamma, dgamma, target, t_init, config)
            %
            % Solves |gamma(t) - target|^2 = 0 for complex t.
            
            if nargin < 5
                config = quadest.util.Config();
            end
            
            R2 = @(t) sum((gamma(t) - target).^2, 2);
            R2_t = @(t) 2 * sum((gamma(t) - target) .* dgamma(t), 2);
            
            % First attempt
            t0 = t_init + config.newtonInitPerturbation;
            tol = config.newtonTol;
            maxIter = config.newtonMaxIter;
            converged = false;
            
            for k = 1:maxIter
                step = R2(t0) / R2_t(t0);
                t0 = t0 - step;
                if abs(step) < tol
                    converged = true;
                    return;
                end
            end
            
            % Fallback with smaller step size
            if ~converged || isnan(real(t0)) || isnan(imag(t0))
                t0 = t_init + config.newtonInitPerturbation;
                maxIter = config.newtonMaxIterFallback;
                stepScale = config.newtonStepSizeFallback;
                
                for k = 1:maxIter
                    step = R2(t0) / R2_t(t0);
                    t0 = t0 - stepScale * step;
                    if abs(step) < tol
                        converged = true;
                        return;
                    end
                end
            end
        end
    end
end
