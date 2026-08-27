classdef AxsymGrid < handle
    % AXSYMGRID Gauss-Legendre × Trapezoidal grid on an axisymmetric geometry
    %
    % Creates a discretization of an axisymmetric surface using:
    %   - Gauss-Legendre nodes in theta direction [0, pi]
    %   - Uniform trapezoidal nodes in phi direction [0, 2*pi)
    %
    % Usage:
    %   geom = quadind.geometry.Spheroid('a', 0.05, 'c', 0.1);
    %   grid = quadind.grid.AxsymGrid(geom, 'nth', 40, 'nph', 60);
    %
    %   % Access grid data
    %   points = grid.x;       % [N×3] surface points
    %   normals = grid.n;      % [N×3] unit normals
    %   weights = grid.w;      % [N×1] quadrature weights
    %
    % See also: quadind.geometry.AxsymGeometry, quadind.grid.Upsampler
    
    properties (SetAccess = private)
        geometry    % AxsymGeometry object
        nth         % Number of theta (GL) points
        nph         % Number of phi (trapezoidal) points
        
        x           % [N×3] surface points (N = nth*nph)
        n           % [N×3] unit normal vectors
        w           % [N×1] quadrature weights (including Jacobian)
        
        theta       % [nth×1] Gauss-Legendre nodes in [0, pi]
        phi         % [1×nph] uniform nodes in [0, 2*pi)
        wtheta      % [nth×1] GL weights
        wphi        % [1×nph] trapezoidal weights
    end
    
    methods
        function obj = AxsymGrid(geometry, varargin)
            % AXSYMGRID Construct a discretized grid on the geometry
            %
            %   grid = AxsymGrid(geometry) uses default nth=40, nph=60
            %   grid = AxsymGrid(geometry, 'nth', N1, 'nph', N2) specifies sizes
            
            arguments
                geometry (1,1) quadind.geometry.AxsymGeometry
            end
            arguments (Repeating)
                varargin
            end
            
            % Parse options
            p = inputParser;
            cfg = quadind.util.Config();
            addParameter(p, 'nth', cfg.nth, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2 && mod(x,1) == 0);
            addParameter(p, 'nph', cfg.nph, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 2 && mod(x,1) == 0);
            parse(p, varargin{:});
            
            obj.geometry = geometry;
            obj.nth = p.Results.nth;
            obj.nph = p.Results.nph;
            
            % Build the grid
            obj.buildGrid();
        end
        
        function N = numPoints(obj)
            % NUMPOINTS Total number of grid points
            N = obj.nth * obj.nph;
        end
        
        function [theta_mat, phi_mat] = meshgrid(obj)
            % MESHGRID Return theta and phi as matrices
            %
            %   [theta_mat, phi_mat] = meshgrid(obj) returns nth×nph matrices
            
            [phi_mat, theta_mat] = meshgrid(obj.phi, obj.theta);
        end
        
        function idx = sub2ind(obj, itheta, iphi)
            % SUB2IND Convert (itheta, iphi) subscripts to linear index
            %
            % Grid is stored with theta varying fastest (column-major).
            
            idx = (iphi - 1) * obj.nth + itheta;
        end
        
        function [itheta, iphi] = ind2sub(obj, idx)
            % IND2SUB Convert linear index to (itheta, iphi) subscripts
            
            itheta = mod(idx - 1, obj.nth) + 1;
            iphi = floor((idx - 1) / obj.nth) + 1;
        end
        
        function newGrid = upsample(obj, factor)
            % UPSAMPLE Create a new grid with higher resolution
            %
            %   newGrid = upsample(obj, factor) creates grid with
            %   factor*nth × factor*nph points
            
            validateattributes(factor, {'numeric'}, {'scalar','integer','positive','finite'});
            newGrid = quadind.grid.AxsymGrid(obj.geometry, ...
                'nth', factor * obj.nth, ...
                'nph', factor * obj.nph);
        end
        
        function interpolator = getInterpolator(obj, targetGrid)
            % GETINTERPOLATOR Get interpolation operator to target grid
            %
            %   interp = getInterpolator(obj, targetGrid)
            %   density_fine = interp.apply(density_coarse)
            %
            % See also: quadind.grid.Upsampler
            
            interpolator = quadind.grid.Upsampler(obj, targetGrid);
        end
        
        function disp(obj)
            % DISP Display grid summary
            fprintf('AxsymGrid: %d×%d = %d points\n', obj.nth, obj.nph, obj.numPoints());
            fprintf('  Geometry: ');
            disp(obj.geometry);
        end
    end
    
    methods (Access = private)
        function buildGrid(obj)
            % BUILDGRID Construct grid points, normals, and weights
            
            % Gauss-Legendre nodes and weights on [0, pi]
            [obj.theta, obj.wtheta] = quadind.util.GaussLegendre(obj.nth, 0, pi);
            
            % Uniform trapezoidal nodes on [0, 2*pi)
            obj.phi = linspace(0, 2*pi, obj.nph + 1);
            obj.phi = obj.phi(1:obj.nph);  % Remove endpoint
            dphi = obj.phi(2) - obj.phi(1);
            obj.wphi = dphi * ones(1, obj.nph);
            
            % Create meshgrid (theta varies fastest)
            [phi_mat, theta_mat] = meshgrid(obj.phi, obj.theta);
            
            % Evaluate surface points
            obj.x = obj.geometry.evaluate(theta_mat(:), phi_mat(:));
            
            % Compute Jacobian and normals
            [J, obj.n] = obj.geometry.jacobian(theta_mat(:), phi_mat(:));
            
            % Quadrature weights: GL weight × trap weight × Jacobian
            % J is column vector, need to reshape to nth × nph
            J_mat = reshape(J, obj.nth, obj.nph);
            W_mat = J_mat .* (obj.wtheta * obj.wphi);
            obj.w = W_mat(:);
            
            % Check that the geometry is adequately resolved by this grid
            quadind.util.Diagnostics.checkGeometryResolution(obj.geometry, obj.nth);
        end
    end
end
