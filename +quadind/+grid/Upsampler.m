classdef Upsampler
    % UPSAMPLER Interpolation operator between two AxsymGrid instances
    %
    % Performs interpolation from a coarse grid to a fine grid using:
    %   - Trigonometric interpolation in phi direction (periodic)
    %   - Barycentric Lagrange interpolation in theta direction (GL nodes)
    %
    % Usage:
    %   coarseGrid = quadind.grid.AxsymGrid(geom, 'nth', 40, 'nph', 60);
    %   fineGrid = coarseGrid.upsample(2);  % 80×120 grid
    %   interp = quadind.grid.Upsampler(coarseGrid, fineGrid);
    %   
    %   density_fine = interp.apply(density_coarse);
    %
    % See also: quadind.grid.AxsymGrid
    
    properties (SetAccess = private)
        sourceGrid      % Source (coarse) grid
        targetGrid      % Target (fine) grid
        T               % Trigonometric interpolation matrix (phi)
        B               % Barycentric Lagrange matrix (theta)
    end
    
    methods
        function obj = Upsampler(sourceGrid, targetGrid)
            % UPSAMPLER Construct interpolation operator
            %
            %   interp = Upsampler(sourceGrid, targetGrid)
            
            arguments
                sourceGrid (1,1) quadind.grid.AxsymGrid
                targetGrid (1,1) quadind.grid.AxsymGrid
            end
            
            obj.sourceGrid = sourceGrid;
            obj.targetGrid = targetGrid;
            
            % Build interpolation matrices
            obj.T = obj.trigInterpMatrix(sourceGrid.phi(:), targetGrid.phi(:));
            
            bclag_weights = obj.bclagWeights(sourceGrid.theta(:));
            obj.B = obj.bclagMatrix(sourceGrid.theta(:), bclag_weights, targetGrid.theta(:));
        end
        
        function density_fine = apply(obj, density_coarse)
            % APPLY Interpolate density from coarse to fine grid
            %
            %   density_fine = apply(obj, density_coarse)
            %
            % Input:
            %   density_coarse - [N_coarse × ncomp] density on source grid
            %
            % Output:
            %   density_fine - [N_fine × ncomp] density on target grid
            
            arguments
                obj
                density_coarse (:,:) {mustBeNumeric}
            end
            
            ncomp = size(density_coarse, 2);
            if size(density_coarse, 1) ~= obj.sourceGrid.numPoints()
                error('quadind:Upsampler:invalidDensity', ...
                    'Density must have %d rows (got %d).', ...
                    obj.sourceGrid.numPoints(), size(density_coarse, 1));
            end
            n_source = obj.sourceGrid.nth;
            m_source = obj.sourceGrid.nph;
            
            % Interpolation: B * reshape(density, nth, nph) * T'
            T_t = obj.T.';
            
            density_fine = zeros(obj.targetGrid.numPoints(), ncomp);
            for c = 1:ncomp
                q_mat = reshape(density_coarse(:, c), n_source, m_source);
                q_fine = obj.B * q_mat * T_t;
                density_fine(:, c) = q_fine(:);
            end
        end
        
        function disp(obj)
            % DISP Display upsampler info
            fprintf('Upsampler: %d×%d → %d×%d\n', ...
                obj.sourceGrid.nth, obj.sourceGrid.nph, ...
                obj.targetGrid.nth, obj.targetGrid.nph);
        end
    end
    
    methods (Static)
        function T = trigInterpMatrix(x, xi)
            % TRIGINTERPMATRIX Build trigonometric interpolation matrix
            %
            % For interpolating from equispaced points x to arbitrary points xi
            % in a periodic domain.
            %
            %   T = trigInterpMatrix(x, xi)
            %   f_interp = T * f
            
            assert(size(x, 2) == 1, 'x must be column vector');
            assert(size(xi, 2) == 1, 'xi must be column vector');
            
            n = numel(x);
            N = numel(xi);
            
            % DFT matrix
            D = fft(eye(n));
            p = ceil(n/2);
            D = [D(p+1:end, :); D(1:p, :)];  % fftshift
            
            % Non-uniform inverse DFT matrix
            I = complex(zeros(N, n));
            q = floor(n/2);
            z = exp(1i * xi);
            
            % Build matrix efficiently
            w = z.^(-q);
            I(:, 1) = w;
            for l = 2:n
                w = w .* z;
                I(:, l) = w;
            end
            I = I / n;
            
            % Combine: real(I * D)
            T = real(I) * real(D) - imag(I) * imag(D);
        end
        
        function w = bclagWeights(x)
            % BCLAGWEIGHTS Barycentric Lagrange interpolation weights
            %
            % Reference: Berrut & Trefethen, SIAM Review 46(3), 2004.
            
            assert(size(x, 2) == 1, 'x must be column vector');
            n = numel(x);
            
            w = zeros(size(x));
            for j = 1:n
                w(j) = 1 / prod(x(j) - x(1:n ~= j));
            end
        end
        
        function B = bclagMatrix(x, w, xi)
            % BCLAGMATRIX Barycentric Lagrange interpolation matrix
            %
            % For interpolating from nodes x (with weights w) to points xi.
            %
            % Reference: Berrut & Trefethen, SIAM Review 46(3), 2004.
            
            assert(size(x, 2) == 1, 'x must be column vector');
            assert(size(xi, 2) == 1, 'xi must be column vector');
            assert(all(size(x) == size(w)), 'x and w must have same size');
            
            n = numel(x);
            N = numel(xi);
            
            B = zeros(N, n);
            denom = zeros(size(xi));
            exact = zeros(size(xi));
            
            for j = 1:n
                xdiff = xi - x(j);
                temp = w(j) ./ xdiff;
                B(:, j) = temp;
                denom = denom + temp;
                exact(xdiff == 0) = j;
            end
            
            B = bsxfun(@rdivide, B, denom);
            
            % Handle exact matches (where xi == x(j))
            jj = find(exact);
            B(jj, :) = 0;
            B(jj + N * (exact(jj) - 1)) = 1;
        end
    end
end
