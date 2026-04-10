classdef Capsule < quadest.geometry.AxsymGeometry
    % CAPSULE Smooth capsule/rod-like axisymmetric geometry
    %
    % A capsule shape with a nearly cylindrical midsection and smoothly
    % rounded endcaps, parameterized as:
    %   at(theta) = R * tanh(kappa * sin(theta)) / tanh(kappa)
    %   ct(theta) = (L/2) * cos(theta)
    %
    % where:
    %   R > 0 is the radius of the cylindrical middle
    %   L > 0 is the total length
    %   kappa > 0 controls how rod-like the shape is (larger = more capsule-like)
    %
    % Usage:
    %   geom = quadest.geometry.Capsule();                    % defaults: R=1, L=6, kappa=4
    %   geom = quadest.geometry.Capsule('R', 0.5, 'L', 4);   % custom parameters
    %
    % The shape is smooth everywhere. At the poles (theta = 0, pi):
    %   at(0) = at(pi) = 0 (surface closes to a point on z-axis)
    %   dadt(0) = R*kappa/tanh(kappa), dadt(pi) = -R*kappa/tanh(kappa)
    % There are no singularities in the parameterization or its derivatives.
    %
    % Unlike the spheroid, the capsule geometry requires Newton iteration
    % to find complex theta roots.
    %
    % See also: quadest.geometry.AxsymGeometry, quadest.geometry.Spheroid, quadest.geometry.Peanut

    properties (SetAccess = private)
        R       % Radius of cylindrical middle
        L       % Total length
        kappa   % Shape parameter (larger = more capsule-like)
    end

    properties (Access = private)
        tanhKappa   % Precomputed tanh(kappa)
    end

    methods
        function obj = Capsule(varargin)
            % CAPSULE Construct a capsule geometry
            %
            %   geom = Capsule() creates capsule with default R=1, L=6, kappa=4
            %   geom = Capsule('R', val, 'L', val, 'kappa', val)

            p = inputParser;
            addParameter(p, 'R', 1, @(x) isscalar(x) && x > 0);
            addParameter(p, 'L', 6, @(x) isscalar(x) && x > 0);
            addParameter(p, 'kappa', 4, @(x) isscalar(x) && x > 0);
            parse(p, varargin{:});

            obj.R = p.Results.R;
            obj.L = p.Results.L;
            obj.kappa = p.Results.kappa;
            obj.tanhKappa = tanh(obj.kappa);
        end

        function val = at(obj, theta)
            % AT Radial distance in xy-plane: R * tanh(kappa*sin(theta)) / tanh(kappa)
            %
            % Smooth and well-defined everywhere:
            %   at(0) = at(pi) = 0 (poles)
            %   at(pi/2) = R (equator)
            val = obj.R * tanh(obj.kappa * sin(theta)) / obj.tanhKappa;
        end

        function val = ct(obj, theta)
            % CT Z-coordinate: (L/2) * cos(theta)
            val = (obj.L / 2) * cos(theta);
        end

        function val = dadt(obj, theta)
            % DADT Derivative of at w.r.t. theta
            % R * kappa * cos(theta) * sech^2(kappa*sin(theta)) / tanh(kappa)
            %
            % Well-defined at poles: dadt(0) = R*kappa/tanh(kappa)
            val = obj.R * obj.kappa * cos(theta) ...
                ./ cosh(obj.kappa * sin(theta)).^2 / obj.tanhKappa;
        end

        function val = dcdt(obj, theta)
            % DCDT Derivative of ct w.r.t. theta: -(L/2) * sin(theta)
            val = -(obj.L / 2) * sin(theta);
        end

        function val = maxRadius(obj)
            % MAXRADIUS Maximum radial extent: R (achieved at theta=pi/2)
            val = obj.R;
        end

        function val = maxHeight(obj)
            % MAXHEIGHT Maximum z extent: L/2 (achieved at theta=0)
            val = obj.L / 2;
        end

        function mask = isExterior(obj, points)
            % ISEXTERIOR Test if points are outside the capsule
            %
            % Uses iterative search to find closest surface point.

            arguments
                obj
                points (:,3) {mustBeNumeric}
            end

            M = size(points, 1);
            distX = sum(points.^2, 2);

            phiVal = mod(atan2(points(:,2), points(:,1)), 2*pi);

            theta = zeros(M, 1);
            for i = 1:M
                Rfun = @(t) sum((obj.evaluate(t, phiVal(i)) - points(i,:)).^2, 2);
                theta(i) = fminbnd(Rfun, 0, pi);
            end

            surfPts = obj.evaluate(theta, phiVal);
            surfDist = sum(surfPts.^2, 2);

            mask = surfDist < distX;
        end

        function disp(obj)
            % DISP Display capsule parameters
            fprintf('Capsule: R = %.4f, L = %.4f, kappa = %.4f\n', ...
                obj.R, obj.L, obj.kappa);
        end
    end
end
