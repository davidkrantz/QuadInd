classdef (Abstract) AxsymGeometry
    % AXSYMGEOMETRY Abstract base class for axisymmetric geometries
    %
    % An axisymmetric geometry is parameterized by functions a(theta) and c(theta)
    % that define the surface in cylindrical coordinates:
    %
    %   x = a(theta) * cos(phi)
    %   y = a(theta) * sin(phi)
    %   z = c(theta)
    %
    % where theta in [0, pi] and phi in [0, 2*pi).
    %
    % Subclasses must implement the abstract methods to define specific shapes.
    %
    % See also: quadest.geometry.Spheroid, quadest.geometry.Peanut
    
    methods (Abstract)
        % Return radial distance in xy-plane: a(theta) * sin(theta) for spheroid
        val = at(obj, theta)
        
        % Return z-coordinate: c(theta) * cos(theta) for spheroid
        val = ct(obj, theta)
        
        % Return derivative d(at)/d(theta)
        val = dadt(obj, theta)
        
        % Return derivative d(ct)/d(theta)
        val = dcdt(obj, theta)
        
        % Return maximum radial extent
        val = maxRadius(obj)
        
        % Return maximum z extent
        val = maxHeight(obj)
        
        % Test if points are exterior to the geometry
        mask = isExterior(obj, points)
    end
    
    methods
        function theta0 = findThetaRoot(obj, targets, phi)
            % FINDTHETAROOT Find complex theta root using Newton iteration
            %
            %   theta0 = findThetaRoot(obj, targets, phi) computes the complex
            %   theta values where the kernel becomes singular.
            %
            % Default implementation using closure-free Newton iteration via
            % RootFinder.newtonSolve. Subclasses with analytic formulas
            % (e.g. Spheroid) should override this method.
            %
            % Inputs:
            %   targets - [M×3] target points
            %   phi     - [M×1] phi values (typically from closest grid point)
            %
            % Output:
            %   theta0 - [M×1] complex theta roots
            
            M = size(targets, 1);
            theta0 = zeros(M, 1);
            
            dfac = pi/2;
            config = quadest.util.Config();
            
            for i = 1:M
                % Find closest real theta as initial guess
                R_real = @(t) sum((obj.evaluate(t, phi(i)) - targets(i,:)).^2, 2);
                theta_init = fminbnd(R_real, 0, pi);
                t_init = 2*theta_init/pi - 1;
                
                % Closure-free Newton solve
                [t0, converged] = quadest.errorest.RootFinder.newtonSolve(...
                    obj, targets(i,:), phi(i), t_init, dfac, config);
                if ~converged
                    error('quadest:AxsymGeometry:thetaRootFailure', ...
                        'Theta-root Newton iteration failed for target [%g %g %g].', ...
                        targets(i,:));
                end
                
                theta0(i) = dfac * (t0 + 1);
            end
        end
        
        function pts = evaluate(obj, theta, phi)
            % EVALUATE Compute surface points at given (theta, phi)
            %
            %   pts = evaluate(obj, theta, phi) returns [N×3] array of
            %   surface points for column vectors theta, phi of same size.
            %
            % This is the parameterization gamma(theta, phi):
            %   [a(theta)*cos(phi), a(theta)*sin(phi), c(theta)]
            
            at_val = obj.at(theta);
            ct_val = obj.ct(theta);
            
            pts = [at_val .* cos(phi), at_val .* sin(phi), ct_val];
        end
        
        function grad = drdtheta(obj, theta, phi)
            % DRDTHETA Derivative of parameterization w.r.t. theta
            %
            %   grad = drdtheta(obj, theta, phi) returns [N×3] array
            
            dadt_val = obj.dadt(theta);
            dcdt_val = obj.dcdt(theta);
            
            grad = [dadt_val .* cos(phi), dadt_val .* sin(phi), dcdt_val];
        end
        
        function grad = drdphi(obj, theta, phi)
            % DRDPHI Derivative of parameterization w.r.t. phi
            %
            %   grad = drdphi(obj, theta, phi) returns [N×3] array
            
            at_val = obj.at(theta);
            
            grad = [-at_val .* sin(phi), at_val .* cos(phi), zeros(size(theta))];
        end
        
        function [J, normals] = jacobian(obj, theta, phi)
            % JACOBIAN Compute Jacobian determinant and unit normals
            %
            %   [J, normals] = jacobian(obj, theta, phi) computes:
            %   - J: Jacobian determinant (surface area element)
            %   - normals: [N×3] unit normal vectors
            %
            % The normal is computed as cross(drdtheta, drdphi) normalized.
            
            at_val = obj.at(theta);
            dadt_val = obj.dadt(theta);
            dcdt_val = obj.dcdt(theta);
            
            % Components of the Jacobian vector (cross product)
            J1 = -dcdt_val .* at_val .* cos(phi);
            J2 = -dcdt_val .* at_val .* sin(phi);
            J3 = dadt_val .* at_val;  % Independent of phi
            
            J = sqrt(J1.^2 + J2.^2 + J3.^2);
            
            % Unit normals
            NX = J1 ./ J;
            NY = J2 ./ J;
            NZ = J3 ./ J;
            
            % Handle poles where Jacobian is zero (NaN normals)
            idx_nan = isnan(NX) | isnan(NY) | isnan(NZ);
            if any(idx_nan(:))
                % North pole: theta ≈ 0
                idxN = idx_nan & abs(theta) < 1e-10;
                NX(idxN) = 0;
                NY(idxN) = 0;
                NZ(idxN) = 1;
                
                % South pole: theta ≈ pi
                idxS = idx_nan & abs(theta - pi) < 1e-10;
                NX(idxS) = 0;
                NY(idxS) = 0;
                NZ(idxS) = -1;
            end
            
            normals = [NX(:), NY(:), NZ(:)];
            
            if any(isnan(normals(:)))
                error('quadest:AxsymGeometry:nanNormals', ...
                    'NaN values in computed normals');
            end
        end
    end
end
