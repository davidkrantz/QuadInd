classdef CustomAxsym < quadest.geometry.AxsymGeometry
    % CUSTOMAXSYM User-defined axisymmetric geometry
    %
    % Allows defining an arbitrary axisymmetric geometry by providing
    % function handles for a(theta), c(theta), and their derivatives.
    %
    % Usage:
    %   % Define a custom shape
    %   a_func = @(theta) 1 + 0.3*cos(5*theta);
    %   c_func = @(theta) 1 + 0.3*cos(5*theta);
    %   da_func = @(theta) -0.3*5*sin(5*theta);
    %   dc_func = @(theta) -0.3*5*sin(5*theta);
    %   
    %   geom = quadest.geometry.CustomAxsym(a_func, c_func, da_func, dc_func);
    %
    % Note: CustomAxsym uses Newton iteration for theta root finding
    % since no analytic formula is available for arbitrary shapes.
    %
    % See also: quadest.geometry.AxsymGeometry, quadest.geometry.Spheroid
    
    properties (SetAccess = private)
        aFunc       % Function handle for a(theta) * sin(theta)
        cFunc       % Function handle for c(theta) * cos(theta)
        daFunc      % Function handle for d/dtheta [a(theta) * sin(theta)]
        dcFunc      % Function handle for d/dtheta [c(theta) * cos(theta)]
        name        % Optional name for display
    end
    
    properties (Access = private)
        maxRadiusCache
        maxHeightCache
    end
    
    methods
        function obj = CustomAxsym(aFunc, cFunc, daFunc, dcFunc, varargin)
            % CUSTOMAXSYM Construct a custom axisymmetric geometry
            %
            %   geom = CustomAxsym(aFunc, cFunc, daFunc, dcFunc)
            %   geom = CustomAxsym(..., 'name', 'MyShape')
            %
            % Inputs:
            %   aFunc  - Function handle: at(theta) = a(theta)*sin(theta)
            %   cFunc  - Function handle: ct(theta) = c(theta)*cos(theta)
            %   daFunc - Function handle: d(at)/d(theta)
            %   dcFunc - Function handle: d(ct)/d(theta)
            %   'name' - Optional name for display
            %
            % Note: The functions should accept vector inputs.
            
            arguments
                aFunc function_handle
                cFunc function_handle
                daFunc function_handle
                dcFunc function_handle
            end
            arguments (Repeating)
                varargin
            end
            
            obj.aFunc = aFunc;
            obj.cFunc = cFunc;
            obj.daFunc = daFunc;
            obj.dcFunc = dcFunc;
            
            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'name', 'CustomAxsym', @ischar);
            parse(p, varargin{:});
            obj.name = p.Results.name;
            
            % Precompute max extents by sampling
            theta_sample = linspace(0, pi, 1000);
            obj.maxRadiusCache = max(abs(obj.aFunc(theta_sample)));
            obj.maxHeightCache = max(abs(obj.cFunc(theta_sample)));
        end
        
        function val = at(obj, theta)
            % AT Radial distance in xy-plane
            val = obj.aFunc(theta);
        end
        
        function val = ct(obj, theta)
            % CT Z-coordinate
            val = obj.cFunc(theta);
        end
        
        function val = dadt(obj, theta)
            % DADT Derivative of at w.r.t. theta
            val = obj.daFunc(theta);
        end
        
        function val = dcdt(obj, theta)
            % DCDT Derivative of ct w.r.t. theta
            val = obj.dcFunc(theta);
        end
        
        function val = maxRadius(obj)
            % MAXRADIUS Maximum radial extent
            val = obj.maxRadiusCache;
        end
        
        function val = maxHeight(obj)
            % MAXHEIGHT Maximum z extent
            val = obj.maxHeightCache;
        end
        
        function mask = isExterior(obj, points)
            % ISEXTERIOR Test if points are outside the geometry
            %
            % Uses iterative search (same as Peanut).
            
            arguments
                obj
                points (:,3) {mustBeNumeric}
            end
            
            M = size(points, 1);
            distX = sum(points.^2, 2);
            
            phiVal = mod(atan2(points(:,2), points(:,1)), 2*pi);
            
            theta = zeros(M, 1);
            for i = 1:M
                R = @(t) sum((obj.evaluate(t, phiVal(i)) - points(i,:)).^2, 2);
                theta(i) = fminbnd(R, 0, pi);
            end
            
            surfPts = obj.evaluate(theta, phiVal);
            surfDist = sum(surfPts.^2, 2);
            
            mask = surfDist < distX;
        end
        
        function disp(obj)
            % DISP Display geometry info
            fprintf('%s: maxRadius = %.4f, maxHeight = %.4f\n', ...
                obj.name, obj.maxRadiusCache, obj.maxHeightCache);
        end
    end
end
