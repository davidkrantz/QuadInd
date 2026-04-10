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
                % Vectorized row-wise dot product (dot() conjugates first arg
                % for complex vectors, matching legacy behavior)
                Gphi = 2 * sum(conj(r) .* dgamma_phi0, 2);
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
                % Vectorized row-wise dot product (dot() conjugates first arg
                % for complex vectors, matching legacy behavior)
                Gtheta = 2 * sum(conj(r) .* dgamma_t0, 2);
            end
        end
        
        function [t0, converged] = newtonSolve(geometry, target, phi_star, t_init, dfac, config)
            % NEWTONSOLVE Closure-free Newton solver for complex theta root
            %
            %   [t0, converged] = newtonSolve(geometry, target, phi_star, t_init, dfac, config)
            %
            % Solves |gamma(t) - target|^2 = 0 for complex t, where
            % gamma(t) = geometry.evaluate(dfac*(t+1), phi_star).
            %
            % This is the canonical Newton solver for all complex theta root
            % finding in the package. It is closure-free: directly calls
            % geometry.evaluate() and geometry.drdtheta() per step, avoiding
            % function handle overhead.
            %
            % Inputs:
            %   geometry - AxsymGeometry object
            %   target   - [1×3] target point
            %   phi_star - scalar phi value
            %   t_init   - scalar initial guess in mapped t-space
            %   dfac     - derivative scaling factor (typically pi/2)
            %   config   - (optional) Config object for solver parameters
            %
            % Outputs:
            %   t0        - complex root in mapped t-space
            %   converged - logical, true if solver converged
            
            if nargin < 6
                config = quadest.util.Config();
            end
            
            tol = config.newtonTol;
            t0 = t_init + config.newtonInitPerturbation;
            converged = false;
            
            % First attempt: standard Newton
            for k = 1:config.newtonMaxIter
                theta_k = dfac * (t0 + 1);
                gval = geometry.evaluate(theta_k, phi_star);
                dgval = dfac * geometry.drdtheta(theta_k, phi_star);
                rv = gval - target;
                R2v = sum(rv.^2, 2);
                R2tv = 2 * sum(rv .* dgval, 2);
                step = R2v / R2tv;
                t0 = t0 - step;
                if abs(step) < tol
                    converged = true;
                    return;
                end
            end
            
            % Fallback with smaller step size
            if ~converged || isnan(real(t0)) || isnan(imag(t0))
                t0 = t_init + config.newtonInitPerturbation;
                stepScale = config.newtonStepSizeFallback;
                
                for k = 1:config.newtonMaxIterFallback
                    theta_k = dfac * (t0 + 1);
                    gval = geometry.evaluate(theta_k, phi_star);
                    dgval = dfac * geometry.drdtheta(theta_k, phi_star);
                    rv = gval - target;
                    R2v = sum(rv.^2, 2);
                    R2tv = 2 * sum(rv .* dgval, 2);
                    step = R2v / R2tv;
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
