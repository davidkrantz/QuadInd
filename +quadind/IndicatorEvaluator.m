classdef IndicatorEvaluator < handle
    % INDICATOREVALUATOR Tabulate and evaluate Stokes stresslet error indicators
    %
    % Precomputes uniform-density quadrature error indicators on a meridional
    % half-plane and evaluates density-dependent indicators at target points.
    % A typical workflow is:
    %
    %   1. Create an axisymmetric geometry and a Stokes stresslet kernel.
    %   2. Construct an IndicatorEvaluator, which performs the tabulation.
    %   3. Obtain its quadrature grid and define a density on that grid.
    %   4. Evaluate indicators, optionally classifying targets for a tolerance.
    %
    % Example:
    %   geometry = quadind.geometry.Spheroid('a', 0.5, 'c', 1);
    %   kernel = quadind.kernel.StokesStresslet();
    %   evaluator = quadind.IndicatorEvaluator(geometry, kernel, ...
    %       'nth', 20, 'nph', 30, 'upsampFactors', 1:3);
    %   grid = evaluator.getGrid();
    %   density = ones(grid.numPoints(), kernel.numComponents());
    %   targets = [0.6, 0, 0; 0.5, 0, 0.5];
    %   [indicators, classification] = evaluator.evaluate( ...
    %       targets, density, 'tol', 1e-6);
    %
    % Indicator construction currently supports only the Stokes stresslet.
    % Because the tabulation is performed for z >= 0 and reused through
    % reflection, the geometry must be symmetric about z = 0.
    %
    % See also: quadind.geometry.AxsymGeometry,
    %   quadind.kernel.StokesStresslet, quadind.grid.AxsymGrid

    properties (SetAccess = private)
        geometry           % Axisymmetric source geometry
        kernel             % Kernel for which indicators are constructed
        grid               % Base surface quadrature grid
        interpolants       % Tabulated uniform-density indicator interpolants
        rootInterpolants   % Cached complex theta-root interpolants, if enabled
        upsampFactors      % Available integer upsampling factors
        config             % Tabulation and numerical configuration
        interpolateRoots   % Whether cached theta roots are used at evaluation
    end

    methods
        function obj = IndicatorEvaluator(geometry, kernel, varargin)
            % INDICATOREVALUATOR Construct an evaluator and tabulate indicators
            %
            %   evaluator = quadind.IndicatorEvaluator(geometry, kernel)
            %   evaluator = quadind.IndicatorEvaluator(geometry, kernel, ...
            %       Name=Value)
            %
            % Name-value options:
            %   config           - quadind.util.Config object
            %   nth              - Number of Gauss-Legendre theta nodes
            %   nph              - Number of equispaced azimuthal nodes
            %   upsampFactors    - Integer factors tabulated for classification
            %   interpolateRoots - Use cached theta-root interpolants
            %
            % Construction performs the potentially expensive precomputation.
            % The supplied geometry is checked for the equatorial symmetry
            % required by the half-plane tabulation.
            arguments
                geometry (1,1) quadind.geometry.AxsymGeometry
                kernel (1,1) quadind.kernel.Kernel
            end
            arguments (Repeating)
                varargin
            end
            if ~isa(kernel, 'quadind.kernel.StokesStresslet')
                error('quadind:IndicatorEvaluator:unsupportedKernel', ...
                    'IndicatorEvaluator currently supports only quadind.kernel.StokesStresslet.');
            end
            quadind.IndicatorEvaluator.validateSymmetry(geometry);

            p = inputParser;
            p.FunctionName = 'quadind.IndicatorEvaluator';
            defaultConfig = quadind.util.Config();
            addParameter(p, 'config', defaultConfig, ...
                @(x) isa(x, 'quadind.util.Config') && isscalar(x));
            addParameter(p, 'nth', defaultConfig.nth, ...
                @(x) quadind.IndicatorEvaluator.isIntegerAtLeast(x, 2));
            addParameter(p, 'nph', defaultConfig.nph, ...
                @(x) quadind.IndicatorEvaluator.isIntegerAtLeast(x, 2));
            addParameter(p, 'upsampFactors', defaultConfig.upsampFactors, ...
                @(x) isnumeric(x) && isvector(x) && all(isfinite(x)) && ...
                all(x >= 1) && all(mod(x, 1) == 0));
            addParameter(p, 'interpolateRoots', defaultConfig.interpolateRoots, ...
                @(x) islogical(x) && isscalar(x));
            parse(p, varargin{:});

            cfg = p.Results.config;
            cfg.nth = p.Results.nth;
            cfg.nph = p.Results.nph;
            cfg.validate();

            obj.geometry = geometry;
            obj.kernel = kernel;
            obj.config = cfg;
            obj.upsampFactors = unique([1, p.Results.upsampFactors], 'sorted');
            obj.interpolateRoots = p.Results.interpolateRoots;
            obj.grid = quadind.grid.AxsymGrid(geometry, 'nth', cfg.nth, 'nph', cfg.nph);
            [obj.interpolants, obj.rootInterpolants] = ...
                quadind.indicator.UniformIndicatorBuilder.build(geometry, obj.grid, ...
                kernel, obj.upsampFactors, obj.interpolateRoots, cfg);
        end

        function [indicators, classification] = evaluate(obj, targets, density, varargin)
            % EVALUATE Compute density-dependent indicators at target points
            %
            %   indicators = evaluator.evaluate(targets, density)
            %   [indicators, classification] = evaluator.evaluate( ...
            %       targets, density, 'tol', tolerance)
            %
            % Inputs:
            %   targets - M-by-3 Cartesian target coordinates
            %   density - N-by-C density values on evaluator.getGrid(), where
            %             N is grid.numPoints() and C is kernel.numComponents()
            %
            % Output:
            %   indicators - M-by-1 quadrature error indicators for the base grid
            %
            % When a tolerance is supplied, classification contains:
            %   upsamplingFactor         - Recommended factor for each target
            %   requiresSpecialQuadrature- True if no tabulated factor suffices
            %   masks                    - Newly accepted targets at each factor
            %   indicators               - Indicator values at each factor
            %   tol                      - The requested tolerance
            arguments
                obj
                targets (:,3) {mustBeNumeric, mustBeFinite}
                density (:,:) {mustBeNumeric, mustBeFinite}
            end
            arguments (Repeating)
                varargin
            end
            obj.validateDensity(density);
            p = inputParser;
            p.FunctionName = 'quadind.IndicatorEvaluator.evaluate';
            addParameter(p, 'tol', [], @(x) isempty(x) || ...
                (isnumeric(x) && isscalar(x) && isfinite(x) && x > 0));
            parse(p, varargin{:});
            tol = p.Results.tol;

            [rquery, zquery] = obj.mapTargets(targets);
            uniformIndicators = obj.interpolateIndicators(1, rquery, zquery);
            if obj.interpolateRoots && ~isempty(obj.rootInterpolants)
                theta0 = obj.rootInterpolants.Re{1}(rquery, zquery) + ...
                    1i * obj.rootInterpolants.Im{1}(rquery, zquery);
                negativeZ = targets(:,3) < 0;
                theta0(negativeZ) = pi - theta0(negativeZ);
                [qphi, qtheta] = quadind.indicator.DensityModifier.computeDensityAtRoots(...
                    obj.grid, targets, density, theta0);
            else
                [qphi, qtheta] = quadind.indicator.DensityModifier.computeDensityAtRoots(...
                    obj.grid, targets, density);
            end
            modifier = quadind.indicator.DensityModifier.computeModifier(qphi, qtheta, targets);
            indicators = sum(uniformIndicators .* modifier, 2);

            if nargout > 1
                if isempty(tol)
                    classification = [];
                else
                    classification = obj.classifyTargets(modifier, tol, rquery, zquery);
                end
            end
        end

        function indicators = evaluateDirect(obj, targets, density, varargin)
            % EVALUATEDIRECT Compute indicators without tabulation interpolation
            %
            %   indicators = evaluator.evaluateDirect(targets, density)
            %
            % Evaluates the same base-grid indicator construction directly at
            % every target. This is useful for validating the interpolated
            % result, but is substantially more expensive than evaluate(). It
            % does not classify targets or accept additional options.
            arguments
                obj
                targets (:,3) {mustBeNumeric, mustBeFinite}
                density (:,:) {mustBeNumeric, mustBeFinite}
            end
            arguments (Repeating)
                varargin
            end
            if ~isempty(varargin)
                error('quadind:IndicatorEvaluator:unsupportedDirectOption', ...
                    'evaluateDirect uses grid density only and does not support additional options.');
            end
            obj.validateDensity(density);
            indicators = quadind.indicator.UniformIndicatorBuilder.computeDirectIndicators(...
                obj.geometry, obj.grid, obj.kernel, targets, density);
        end

        function grid = getGrid(obj)
            % GETGRID Return the base quadrature grid used by the evaluator
            grid = obj.grid;
        end

        function disp(obj)
            % DISP Display a summary of the evaluator and its tabulation
            fprintf('IndicatorEvaluator:\n');
            fprintf('  Geometry: '); disp(obj.geometry);
            fprintf('  Kernel: %s\n', obj.kernel.kernelName());
            fprintf('  Grid: %d x %d = %d points\n', obj.grid.nth, obj.grid.nph, obj.grid.numPoints());
            fprintf('  Upsampling factors: [%s]\n', num2str(obj.upsampFactors));
            fprintf('  Interpolate roots: %s\n', mat2str(obj.interpolateRoots));
        end
    end

    methods (Access = private)
        function validateDensity(obj, density)
            % Validate that density matches the base grid and kernel layout.
            if size(density, 1) ~= obj.grid.numPoints() || ...
                    size(density, 2) ~= obj.kernel.numComponents()
                error('quadind:IndicatorEvaluator:invalidDensity', ...
                    'Density must have size %d-by-%d (got %d-by-%d).', ...
                    obj.grid.numPoints(), obj.kernel.numComponents(), size(density,1), size(density,2));
            end
        end

        function [rquery, zquery] = mapTargets(obj, targets)
            % Map Cartesian targets to the symmetric tabulation half-plane.
            rxy = hypot(targets(:,1), targets(:,2));
            zabs = abs(targets(:,3));
            extent = 2 * (obj.geometry.maxRadius() + obj.geometry.maxHeight());
            rquery = min(rxy, extent);
            zquery = min(zabs, extent);
        end

        function values = interpolateIndicators(obj, upsamplingFactor, rquery, zquery)
            % Recover per-component indicators from their log10 interpolants.
            values = zeros(numel(rquery), obj.kernel.numComponents());
            for component = 1:size(values, 2)
                values(:,component) = 10.^obj.interpolants{upsamplingFactor,component}(rquery, zquery);
            end
            values = obj.config.tabulationSafetyFactor * values;
        end

        function classification = classifyTargets(obj, modifier, tol, rquery, zquery)
            % Assign each target to the first tabulated factor below tolerance.
            targetCount = numel(rquery);
            factorCount = numel(obj.upsampFactors);
            classification.upsamplingFactor = ones(targetCount, 1);
            classification.requiresSpecialQuadrature = false(targetCount, 1);
            classification.masks = cell(1, factorCount);
            classification.indicators = cell(1, factorCount);
            classification.tol = tol;
            available = true(targetCount, 1);
            for factorIndex = 1:factorCount
                upsamplingFactor = obj.upsampFactors(factorIndex);
                values = obj.interpolateIndicators(upsamplingFactor, rquery, zquery);
                indicator = sum(values .* modifier, 2);
                mask = available & indicator < tol;
                classification.indicators{factorIndex} = indicator;
                classification.masks{factorIndex} = mask;
                classification.upsamplingFactor(mask) = obj.upsampFactors(factorIndex);
                available(mask) = false;
            end
            classification.requiresSpecialQuadrature = available;
            classification.upsamplingFactor(available) = max(obj.upsampFactors) + 1;
        end
    end

    methods (Static, Access = private)
        function tf = isIntegerAtLeast(value, lowerBound)
            tf = isnumeric(value) && isscalar(value) && isfinite(value) && ...
                value >= lowerBound && mod(value, 1) == 0;
        end

        function validateSymmetry(geometry)
            theta = linspace(0, pi, 65).';
            reflected = pi - theta;
            differences = {geometry.at(theta) - geometry.at(reflected), ...
                geometry.ct(theta) + geometry.ct(reflected), ...
                geometry.dadt(theta) + geometry.dadt(reflected), ...
                geometry.dcdt(theta) - geometry.dcdt(reflected)};
            scales = {geometry.at(theta), geometry.ct(theta), ...
                geometry.dadt(theta), geometry.dcdt(theta)};
            for k = 1:numel(differences)
                scale = max(1, max(abs(scales{k})));
                if any(~isfinite(differences{k})) || max(abs(differences{k})) > 1e-9 * scale
                    error('quadind:IndicatorEvaluator:asymmetricGeometry', ...
                        ['The half-plane indicator tabulation requires equatorial symmetry: ', ...
                         'a(pi-theta)=a(theta) and c(pi-theta)=-c(theta), with matching derivative parity.']);
                end
            end
        end
    end
end
