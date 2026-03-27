classdef Peanut < quadest.geometry.AxsymGeometry
    % PEANUT Peanut-shaped axisymmetric geometry
    %
    % A peanut shape with parameterization:
    %   a(theta) = c(theta) = amplitude + cos(2*theta)
    %
    % where amplitude controls the size (default 2.2).
    %
    % Usage:
    %   geom = quadest.geometry.Peanut();           % default amplitude=2.2
    %   geom = quadest.geometry.Peanut('amplitude', 3.0);
    %
    % Unlike the spheroid, the peanut geometry requires Newton iteration
    % to find complex theta roots.
    %
    % See also: quadest.geometry.AxsymGeometry, quadest.geometry.Spheroid
    
    properties (SetAccess = private)
        amplitude   % Base amplitude (default 2.2)
    end
    
    properties (Access = private)
        config      % Configuration for Newton solver
    end
    
    methods
        function obj = Peanut(varargin)
            % PEANUT Construct a peanut geometry
            %
            %   geom = Peanut() creates peanut with default amplitude=2.2
            %   geom = Peanut('amplitude', val) specifies amplitude
            
            % Parse inputs
            p = inputParser;
            addParameter(p, 'amplitude', 2.2, @(x) isscalar(x) && x > 0);
            parse(p, varargin{:});
            
            obj.amplitude = p.Results.amplitude;
            obj.config = quadest.util.Config();
        end
        
        function val = baseRadius(obj, theta)
            % BASERADIUS Base radius function: amplitude + cos(2*theta)
            val = obj.amplitude + cos(2*theta);
        end
        
        function val = dbaseRadius(obj, theta)
            % DBASERADIUS Derivative of base radius: -2*sin(2*theta)
            val = -2 * sin(2*theta);
        end
        
        function val = at(obj, theta)
            % AT Radial distance in xy-plane: baseRadius * sin(theta)
            val = obj.baseRadius(theta) .* sin(theta);
        end
        
        function val = ct(obj, theta)
            % CT Z-coordinate: baseRadius * cos(theta)
            val = obj.baseRadius(theta) .* cos(theta);
        end
        
        function val = dadt(obj, theta)
            % DADT Derivative of at w.r.t. theta
            % d/dtheta [r(theta)*sin(theta)] = r'*sin + r*cos
            r = obj.baseRadius(theta);
            dr = obj.dbaseRadius(theta);
            val = dr .* sin(theta) + r .* cos(theta);
        end
        
        function val = dcdt(obj, theta)
            % DCDT Derivative of ct w.r.t. theta
            % d/dtheta [r(theta)*cos(theta)] = r'*cos - r*sin
            r = obj.baseRadius(theta);
            dr = obj.dbaseRadius(theta);
            val = dr .* cos(theta) - r .* sin(theta);
        end
        
        function val = maxRadius(obj)
            % MAXRADIUS Maximum radial extent in xy-plane
            % Maximum of (amplitude + cos(2*theta)) * sin(theta)
            % Approximate by sampling
            theta_sample = linspace(0, pi, 1000);
            val = max(obj.at(theta_sample));
        end
        
        function val = maxHeight(obj)
            % MAXHEIGHT Maximum z extent
            % Maximum of (amplitude + cos(2*theta)) * cos(theta)
            theta_sample = linspace(0, pi, 1000);
            val = max(abs(obj.ct(theta_sample)));
        end
        
        function mask = isExterior(obj, points)
            % ISEXTERIOR Test if points are outside the peanut
            %
            % Uses iterative search to find closest surface point.
            
            arguments
                obj
                points (:,3) {mustBeNumeric}
            end
            
            M = size(points, 1);
            distX = sum(points.^2, 2);  % Distance from origin
            
            % Determine phi from (x,y) coordinates
            phiVal = mod(atan2(points(:,2), points(:,1)), 2*pi);
            
            % Find closest point on surface for each target
            theta = zeros(M, 1);
            for i = 1:M
                R = @(t) sum((obj.evaluate(t, phiVal(i)) - points(i,:)).^2, 2);
                theta(i) = fminbnd(R, 0, pi);
            end
            
            % Compare surface distance to point distance
            surfPts = obj.evaluate(theta, phiVal);
            surfDist = sum(surfPts.^2, 2);
            
            mask = surfDist < distX;
        end
        
        function theta0 = findThetaRoot(obj, targets, phi)
            % FINDTHETAROOT Find complex theta root using Newton iteration
            %
            %   theta0 = findThetaRoot(obj, targets, phi) computes the complex
            %   theta value where the kernel becomes singular.
            %
            % Uses Newton solver since no analytic formula exists for peanut.
            %
            % Inputs:
            %   targets - [M×3] target points
            %   phi     - [M×1] phi values (typically from closest grid point)
            %
            % Output:
            %   theta0 - [M×1] complex theta roots
            
            arguments
                obj
                targets (:,3) {mustBeNumeric}
                phi (:,1) {mustBeNumeric}
            end
            
            M = size(targets, 1);
            theta0 = zeros(M, 1);
            
            % Use linear map t = 2*theta/pi - 1 for Newton solver
            imap = @(t) (pi/2) * (t + 1);  % inverse map: t -> theta
            dfac = pi/2;                    % derivative factor
            
            for i = 1:M
                % Parameterization in t-space
                gamma = @(t) obj.evaluate(imap(t), phi(i));
                dgamma = @(t) dfac * obj.drdtheta(imap(t), phi(i));
                
                % Find closest real theta first (as initial guess)
                R_real = @(t) sum((obj.evaluate(t, phi(i)) - targets(i,:)).^2, 2);
                theta_init = fminbnd(R_real, 0, pi);
                t_init = 2*theta_init/pi - 1;
                
                % Newton solve with complex perturbation
                t0 = obj.newtonSolve(gamma, dgamma, targets(i,:), t_init + obj.config.newtonInitPerturbation);
                
                theta0(i) = imap(t0);
            end
        end
        
        function disp(obj)
            % DISP Display peanut parameters
            fprintf('Peanut: amplitude = %.4f\n', obj.amplitude);
        end
    end
    
    methods (Access = private)
        function t0 = newtonSolve(obj, gamma, dgamma, target, t_init)
            % NEWTONSOLVE Newton iteration for complex root finding
            %
            % Solves |gamma(t) - target|^2 = 0 for complex t
            
            % Functional and its derivative
            R2 = @(t) sum((gamma(t) - target).^2, 2);
            R2_t = @(t) 2 * sum((gamma(t) - target) .* dgamma(t), 2);
            
            % First attempt: standard Newton
            t0 = t_init;
            tol = obj.config.newtonTol;
            maxIter = obj.config.newtonMaxIter;
            
            for k = 1:maxIter
                step = R2(t0) / R2_t(t0);
                t0 = t0 - step;
                if abs(step) < tol
                    return;
                end
            end
            
            % If failed, try again with smaller step size
            if abs(step) >= tol || isnan(real(t0)) || isnan(imag(t0))
                quadest.util.Diagnostics.warn('Newton solver: retrying with smaller step size');
                
                t0 = t_init;
                maxIter = obj.config.newtonMaxIterFallback;
                stepScale = obj.config.newtonStepSizeFallback;
                
                for k = 1:maxIter
                    step = R2(t0) / R2_t(t0);
                    t0 = t0 - stepScale * step;
                    if abs(step) < tol
                        return;
                    end
                end
                
                if abs(step) >= tol
                    quadest.util.Diagnostics.warn('Newton solver did not converge');
                end
            end
        end
    end
end
