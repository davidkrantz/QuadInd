classdef ErrorEstimator < handle
    % ERRORESTIMATOR Main class for adaptive quadrature error estimation
    %
    % Precomputes and evaluates error estimates for layer potentials on
    % axisymmetric geometries. The typical workflow is:
    %
    %   1. Create geometry and kernel
    %   2. Construct ErrorEstimator (performs precomputation)
    %   3. Evaluate error estimates at target points with given density
    %
    % Usage:
    %   geom = quadest.geometry.Spheroid('a', 0.05, 'c', 0.1);
    %   kernel = quadest.kernel.StokesStresslet();
    %   estimator = quadest.errorest.ErrorEstimator(geom, kernel, ...
    %       'nth', 40, 'nph', 60, 'upsampFactors', 1:6);
    %
    %   density = ones(40*60, 3);  % stresslet identity
    %   targets = [0.1, 0, 0; 0.05, 0, 0.05];
    %
    %   % Get raw estimates
    %   estimates = estimator.evaluate(targets, density);
    %
    %   % Get estimates with classification
    %   [estimates, classification] = estimator.evaluate(targets, density, 'tol', 1e-6);
    %
    % See also: quadest.geometry.AxsymGeometry, quadest.kernel.Kernel
    
    properties (SetAccess = private)
        geometry            % AxsymGeometry object
        kernel              % Kernel object
        grid                % AxsymGrid object (base discretization)
        
        interpolants        % Cell array of griddedInterpolant for error estimates
        rootInterpolants    % Struct with theta root interpolants (optional)
        upsampFactors       % Array of upsampling factors
        
        config              % Configuration object
        interpolateRoots    % Whether theta roots are interpolated
    end
    
    methods
        function obj = ErrorEstimator(geometry, kernel, varargin)
            % ERRORESTIMATOR Construct error estimator and perform precomputation
            %
            %   estimator = ErrorEstimator(geometry, kernel)
            %   estimator = ErrorEstimator(geometry, kernel, 'nth', 40, 'nph', 60, ...)
            %
            % Options:
            %   'nth'              - Number of theta (GL) points (default: 40)
            %   'nph'              - Number of phi (trapezoidal) points (default: 60)
            %   'upsampFactors'    - Upsampling factors to precompute (default: 1:6)
            %   'interpolateRoots' - Cache theta roots for fast evaluation (default: true)
            %   'config'           - Config object with additional settings
            
            arguments
                geometry (1,1) quadest.geometry.AxsymGeometry
                kernel (1,1) quadest.kernel.Kernel
            end
            arguments (Repeating)
                varargin
            end
            
            % Parse options
            p = inputParser;
            defaultConfig = quadest.util.Config();
            addParameter(p, 'nth', defaultConfig.nth);
            addParameter(p, 'nph', defaultConfig.nph);
            addParameter(p, 'upsampFactors', defaultConfig.upsampFactors);
            addParameter(p, 'interpolateRoots', defaultConfig.interpolateRoots);
            addParameter(p, 'config', defaultConfig);
            parse(p, varargin{:});
            
            obj.geometry = geometry;
            obj.kernel = kernel;
            if ~isa(kernel, 'quadest.kernel.StokesStresslet')
                error('quadest:ErrorEstimator:unsupportedKernel', ...
                    'ErrorEstimator currently supports only StokesStresslet.');
            end
            obj.config = p.Results.config;
            obj.config.nth = p.Results.nth;
            obj.config.nph = p.Results.nph;
            if ~isscalar(obj.config.tabulationSafetyFactor) || ...
                    ~isfinite(obj.config.tabulationSafetyFactor) || ...
                    obj.config.tabulationSafetyFactor < 1
                error('quadest:ErrorEstimator:invalidTabulationSafetyFactor', ...
                    'tabulationSafetyFactor must be a finite scalar >= 1.');
            end
            obj.upsampFactors = unique([1, p.Results.upsampFactors]);
            obj.interpolateRoots = p.Results.interpolateRoots;
            
            % Create quadrature grid
            obj.grid = quadest.grid.AxsymGrid(geometry, ...
                'nth', obj.config.nth, 'nph', obj.config.nph);
            
            % Precompute error estimate interpolants
            [obj.interpolants, obj.rootInterpolants, ~] = ...
                quadest.errorest.UniformEstimateBuilder.build(...
                    geometry, obj.grid, kernel, obj.upsampFactors, ...
                    obj.interpolateRoots, obj.config);
        end
        
        function [result, varargout] = evaluate(obj, targets, density, varargin)
            % EVALUATE Compute error estimates at target points
            %
            %   estimates = evaluate(obj, targets, density)
            %   [estimates, classification] = evaluate(obj, targets, density, 'tol', 1e-6)
            %
            % Inputs:
            %   targets - [M×3] target points
            %   density - [N×3] density on grid (N = nth*nph)
            %
            % Options:
            %   'tol' - Error tolerance for classification (optional)
            %
            % Outputs:
            %   estimates - [M×1] componentwise infinity-norm error indicator
            %
            % If 'tol' is provided, also returns classification struct:
            %   classification.upsampfac - [M×1] recommended upsampling factor
            %   classification.isSQ      - [M×1] logical, true if needs special quadrature
            %   classification.masks     - Cell array of logical masks per upsampling level
            
            arguments
                obj
                targets (:,3) {mustBeNumeric}
                density (:,:) {mustBeNumeric}
            end
            arguments (Repeating)
                varargin
            end
            
            % Parse options
            p = inputParser;
            addParameter(p, 'tol', []);
            parse(p, varargin{:});
            tol = p.Results.tol;
            
            % Validate density size
            expectedN = obj.grid.numPoints();
            if size(density, 1) ~= expectedN
                error('quadest:ErrorEstimator:invalidDensity', ...
                    'Density must have %d rows (got %d)', expectedN, size(density, 1));
            end
            if size(density, 2) ~= obj.kernel.numComponents()
                error('quadest:ErrorEstimator:invalidDensity', ...
                    'Density must have %d columns (got %d)', ...
                    obj.kernel.numComponents(), size(density, 2));
            end
            
            M = size(targets, 1);
            ncomp = obj.kernel.numComponents();
            
            % Map targets to interpolation coordinates
            rxy = sqrt(targets(:,1).^2 + targets(:,2).^2);
            zabs = abs(targets(:,3));
            maxExtent = 2 * (obj.geometry.maxRadius() + obj.geometry.maxHeight());
            insideTabulation = rxy <= maxExtent & zabs <= maxExtent;
            % Holding the boundary value is conservative for the far field
            % and avoids uncontrolled linear extrapolation of log indicators.
            rquery = min(rxy, maxExtent);
            zquery = min(zabs, maxExtent);
            
            % Interpolate uniform estimates (for direct quadrature, upfac=1)
            estval = zeros(M, ncomp);
            for c = 1:ncomp
                estval(:, c) = 10.^(obj.interpolants{1, c}(rquery, zquery));
            end
            estval = obj.config.tabulationSafetyFactor * estval;
            
            % Compute density modification
            if obj.interpolateRoots && ~isempty(obj.rootInterpolants)
                % Use interpolated theta roots
                theta0_interp = obj.rootInterpolants.Re{1}(rquery, zquery) + ...
                                1i * obj.rootInterpolants.Im{1}(rquery, zquery);
                % Root tables use the z >= 0 half-plane. Map the cached root
                % back to the physical hemisphere for reflected targets.
                negz = targets(:,3) < 0;
                theta0_interp(negz) = pi - theta0_interp(negz);

                % % Extra check: bad root
                % [~, iphi] = quadest.errorest.UniformEstimateBuilder.findNearestNodes( ...
                %     obj.grid, targets);
                % phi_star = obj.grid.phi(iphi).';
                % rootResidual = Inf(M,1);
                % if any(insideTabulation)
                %     displacement = obj.geometry.evaluate( ...
                %         theta0_interp(insideTabulation), phi_star(insideTabulation)) ...
                %         - targets(insideTabulation,:);
                %     rootResidual(insideTabulation) = abs(sum(displacement.^2,2)) ./ ...
                %         max(sum(abs(displacement).^2,2), realmin);
                % end
                % invalidRoot = ~insideTabulation | ~isfinite(theta0_interp) | ...
                %     rootResidual > obj.config.interpolatedRootResidualTolerance;
                % if any(invalidRoot)
                %     theta0_interp(invalidRoot) = obj.geometry.findThetaRoot( ...
                %         targets(invalidRoot,:), phi_star(invalidRoot));
                % end

                [qphi, qtheta, ~, ~] = quadest.errorest.DensityModifier.computeDensityAtRoots(...
                    obj.grid, targets, density, theta0_interp);
            else
                % Compute theta roots on-the-fly
                [qphi, qtheta, ~, ~] = quadest.errorest.DensityModifier.computeDensityAtRoots(...
                    obj.grid, targets, density);
            end
            
            % Modification factor: max over phi and theta directions
            modifier = quadest.errorest.DensityModifier.computeModifier(qphi, qtheta, targets);

            % Modified estimate: conservatively sum component indicators.
            result = sum(estval .* modifier, 2);
            % Legacy/direct-mask parity version: max over components.
            %result = max(estval .* modifier, [], 2);
            
            % Classification (only if tolerance provided)
            if ~isempty(tol) && nargout > 1
                classification = obj.classifyTargets(targets, density, modifier, tol, rquery, zquery);
            elseif nargout > 1
                classification = [];
            end
            
            if nargout > 1
                % Return both estimates and classification
                varargout{1} = classification;
            end
        end
        
        function grid = getGrid(obj)
            % GETGRID Return the underlying quadrature grid
            %
            %   grid = estimator.getGrid()
            %
            % Useful for kernel evaluation with the same discretization.
            
            grid = obj.grid;
        end
        
        function disp(obj)
            % DISP Display estimator summary
            fprintf('ErrorEstimator:\n');
            fprintf('  Geometry: ');
            disp(obj.geometry);
            fprintf('  Kernel: %s\n', obj.kernel.kernelName());
            fprintf('  Grid: %d×%d = %d points\n', ...
                obj.config.nth, obj.config.nph, obj.grid.numPoints());
            fprintf('  Upsampling factors: [%s]\n', num2str(obj.upsampFactors));
            fprintf('  Interpolate roots: %s\n', mat2str(obj.interpolateRoots));
        end
        
        function result = evaluateDirect(obj, targets, density, varargin)
            % EVALUATEDIRECT Compute error estimates directly (without tabulation)
            %
            %   estimates = evaluateDirect(obj, targets, density)
            %
            % Computes the same per-component uniform indicators used to build
            % the tabulated interpolants, but evaluates them directly at each
            % target point instead of interpolating. The supplied grid density is
            % then applied using the same root-density modifier as evaluate().
            %
            % Unlike evaluate(), this method does not support classification
            % or upsampling factors — it only computes direct (upfac=1) estimates.
            %
            % Inputs:
            %   targets - [M×3] target points
            %   density - [N×ncomp] density on grid (N = nth*nph)
            %
            % Outputs:
            %   estimates - [M×1] componentwise infinity-norm error indicator
            %
            % See also: evaluate, quadest.errorest.UniformEstimateBuilder.computeDirectEstimates
            
            arguments
                obj
                targets (:,3) {mustBeNumeric}
                density (:,:) {mustBeNumeric}
            end
            arguments (Repeating)
                varargin
            end

            if ~isempty(varargin)
                error('quadest:ErrorEstimator:unsupportedDirectOption', ...
                    'evaluateDirect uses grid density only and does not support additional options.');
            end

            expectedN = obj.grid.numPoints();
            if size(density, 1) ~= expectedN
                error('quadest:ErrorEstimator:invalidDensity', ...
                    'Density must have %d rows (got %d)', expectedN, size(density, 1));
            end
            if size(density, 2) ~= obj.kernel.numComponents()
                error('quadest:ErrorEstimator:invalidDensity', ...
                    'Density must have %d columns (got %d)', ...
                    obj.kernel.numComponents(), size(density, 2));
            end
            
            result = quadest.errorest.UniformEstimateBuilder.computeDirectEstimates(...
                obj.geometry, obj.grid, obj.kernel, targets, density);
        end
    end
    
    methods (Access = private)
        function classification = classifyTargets(obj, targets, density, modifier, tol, rxy, zabs)
            % CLASSIFYTARGETS Assign targets to quadrature regions based on estimates
            
            M = size(targets, 1);
            ncomp = obj.kernel.numComponents();
            nfac = length(obj.upsampFactors);
            
            % Initialize classification
            classification.upsampfac = ones(M, 1);
            classification.isSQ = false(M, 1);
            classification.masks = cell(1, nfac);
            classification.estimates = cell(1, nfac);
            classification.tol = tol;
            
            % Track which targets are still available (not yet classified)
            available = true(M, 1);
            
            for k = 1:nfac
                upfac = obj.upsampFactors(k);
                
                % Interpolate estimate for this upsampling factor
                estval = zeros(M, ncomp);
                for c = 1:ncomp
                    estval(:, c) = 10.^(obj.interpolants{upfac, c}(rxy, zabs));
                end
                estval = obj.config.tabulationSafetyFactor * estval;
                
                % Modified estimate
                estimate_up = sum(estval .* modifier, 2);
                % Legacy/direct-mask parity version: max over components.
                %estimate_up = max(estval .* modifier, [], 2);
                
                % Store per-factor estimate for all targets
                classification.estimates{k} = estimate_up;
                
                % Classify: targets where estimate < tol
                mask = available & (estimate_up < tol);
                classification.masks{k} = mask;
                classification.upsampfac(mask) = upfac;
                
                % Update available targets
                available = available & ~mask;
            end
            
            % Remaining targets require special quadrature
            classification.isSQ = available;
            classification.upsampfac(available) = max(obj.upsampFactors) + 1;  % Flag as SQ
        end
    end
end
