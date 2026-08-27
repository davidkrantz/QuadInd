classdef Peanut < quadind.geometry.AxsymGeometry
    % PEANUT Peanut-shaped axisymmetric geometry
    %
    % A peanut shape with parameterization:
    %   a(theta) = c(theta) = amplitude + cos(2*theta)
    %
    % where amplitude controls the size (default 2.2).
    %
    % Usage:
    %   geom = quadind.geometry.Peanut();           % default amplitude=2.2
    %   geom = quadind.geometry.Peanut('amplitude', 3.0);
    %
    % Unlike the spheroid, the peanut geometry requires Newton iteration
    % to find complex theta roots.
    %
    % See also: quadind.geometry.AxsymGeometry, quadind.geometry.Spheroid
    
    properties (SetAccess = private)
        amplitude   % Base amplitude (default 2.2)
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
        end
        
        function val = baseRadius(obj, theta)
            % BASERADIUS Base radius function: amplitude + cos(2*theta)
            val = obj.amplitude + cos(2*theta);
        end
        
        function val = dbaseRadius(~, theta)
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
        
        function disp(obj)
            % DISP Display peanut parameters
            fprintf('Peanut: amplitude = %.4f\n', obj.amplitude);
        end
    end
end
