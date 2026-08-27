classdef Config
    % CONFIG Validated parameters for the QuadInd package.
    %
    % This class provides default values for all configurable parameters
    % used throughout the error-indicator evaluation workflow.
    %
    % Usage:
    %   cfg = quadind.util.Config();
    %   cfg.tol = 1e-8;  % Override default tolerance
    %
    % See also: quadind.IndicatorEvaluator
    
    properties
        % Quadrature grid defaults
        nth = 40              % Number of theta (GL) points
        nph = 60              % Number of phi (trapezoidal) points
        
        % Error-indicator evaluation defaults
        tol = 1e-6            % Default error tolerance for classification
        upsampFactors = 1:6   % Upsampling factors to precompute
        interpolateRoots = true  % Whether to interpolate theta roots
        interpolatedRootResidualTolerance = NaN  % Root interpolation/fallback check
        tabulationSafetyFactor = 1  % Guard against log-interpolation underestimation
        
        % Tabulation grid for precomputation
        ntab = 100            % Tabulation grid size (xy direction)
        nztab = 100           % Tabulation grid size (z direction)
        
        % Newton solver parameters
        newtonTol = 1e-12     % Convergence tolerance
        newtonMaxIter = 100   % Maximum iterations (first attempt)
        newtonMaxIterFallback = 2000  % Maximum iterations (fallback)
        newtonStepSizeFallback = 0.1  % Step size for fallback
        newtonInitPerturbation = 0.1i  % Initial complex perturbation
        
        % Numerical safety
        interiorIndicator = 1/eps^2  % Indicator assigned to interior table points
        indicatorLowerBound = eps    % Positive floor before log10 tabulation
    end
    
    properties (Constant)
        % Stresslet singularity order (fixed by kernel)
        stressletSingularityOrder = 5/2

        % Density components (fixed for vector problems)
        numComponents = 3
        
        % Package version
        version = '0.1.0'
    end
    
    methods
        function obj = Config(varargin)
            % CONFIG Construct configuration object
            %
            %   cfg = Config() creates config with defaults
            %   cfg = Config('nth', 80, 'tol', 1e-8) overrides specific values
            
            if mod(nargin, 2) ~= 0
                error('quadind:Config:invalidArguments', ...
                    'Configuration values must be supplied as name-value pairs.');
            end
            for i = 1:2:length(varargin)
                name = varargin{i};
                if ~(ischar(name) || (isstring(name) && isscalar(name)))
                    error('quadind:Config:invalidName', 'Configuration names must be text scalars.');
                end
                name = char(name);
                constantNames = {'stressletSingularityOrder','numComponents','version'};
                if isprop(obj, name) && ~ismember(name, constantNames)
                    obj.(name) = varargin{i+1};
                else
                    error('quadind:Config:unknownProperty', 'Unknown configuration property: %s.', name);
                end
            end
            obj.validate();
        end

        function validate(obj)
            validateattributes(obj.nth, {'numeric'}, {'scalar','integer','>=',2,'finite'}, 'Config', 'nth');
            validateattributes(obj.nph, {'numeric'}, {'scalar','integer','>=',2,'finite'}, 'Config', 'nph');
            validateattributes(obj.ntab, {'numeric'}, {'scalar','integer','>=',2,'finite'}, 'Config', 'ntab');
            validateattributes(obj.nztab, {'numeric'}, {'scalar','integer','>=',2,'finite'}, 'Config', 'nztab');
            validateattributes(obj.tol, {'numeric'}, {'scalar','positive','finite'}, 'Config', 'tol');
            validateattributes(obj.upsampFactors, {'numeric'}, {'vector','integer','positive','finite'}, 'Config', 'upsampFactors');
            validateattributes(obj.interpolateRoots, {'logical'}, {'scalar'}, 'Config', 'interpolateRoots');
            validateattributes(obj.tabulationSafetyFactor, {'numeric'}, {'scalar','finite','>=',1}, 'Config', 'tabulationSafetyFactor');
            if ~(isnumeric(obj.interpolatedRootResidualTolerance) && ...
                    isscalar(obj.interpolatedRootResidualTolerance) && ...
                    (isnan(obj.interpolatedRootResidualTolerance) || ...
                    (isfinite(obj.interpolatedRootResidualTolerance) && obj.interpolatedRootResidualTolerance > 0)))
                error('quadind:Config:invalidRootTolerance', ...
                    'interpolatedRootResidualTolerance must be NaN or a positive finite scalar.');
            end
        end
        
        function s = toStruct(obj)
            % TOSTRUCT Convert configuration to struct for saving
            props = properties(obj);
            s = struct();
            for i = 1:length(props)
                s.(props{i}) = obj.(props{i});
            end
        end
        
        function disp(obj)
            % DISP Display configuration summary
            fprintf('quadind.util.Config (v%s)\n', obj.version);
            fprintf('  Grid: nth=%d, nph=%d\n', obj.nth, obj.nph);
            fprintf('  Tolerance: %.2e\n', obj.tol);
            fprintf('  Upsampling factors: [%s]\n', num2str(obj.upsampFactors));
            fprintf('  Interpolate roots: %s\n', mat2str(obj.interpolateRoots));
        end
    end
end
