classdef Config
    % CONFIG Default parameters and configuration for quadest package
    %
    % This class provides default values for all configurable parameters
    % used throughout the error estimation workflow.
    %
    % Usage:
    %   cfg = quadest.util.Config();
    %   cfg.tol = 1e-8;  % Override default tolerance
    %
    % See also: quadest.errorest.ErrorEstimator
    
    properties
        % Quadrature grid defaults
        nth = 40              % Number of theta (GL) points
        nph = 60              % Number of phi (trapezoidal) points
        
        % Error estimation defaults
        tol = 1e-6            % Default error tolerance for classification
        upsampFactors = 1:6   % Upsampling factors to precompute
        interpolateRoots = true  % Whether to interpolate theta roots
        interpolatedRootResidualTolerance = NaN  % Root interpolation diagnostic
        
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
        interiorEstimate = 1/eps^2  % Error estimate for interior points
        estimateLowerBound = eps    % Lower bound on error estimates
        
        % Gauss-Laguerre quadrature order
        laguerreOrder = 8     % 8-point Gauss-Laguerre (hardcoded)
        
        % Density resolution check
        densityDecayThreshold = 1e-12  % FFT decay threshold
        densityDecayModes = 5          % Number of high-frequency modes to check
        
        % Reference solution upsampling
        sqUpsampFactor = 20   % Upsampling factor for SQ region
        refUpsampFactor = 10  % Upsampling factor for reference solutions
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
            
            % Parse name-value pairs
            for i = 1:2:length(varargin)
                if isprop(obj, varargin{i})
                    obj.(varargin{i}) = varargin{i+1};
                else
                    warning('quadest:Config:unknownProperty', ...
                        'Unknown property: %s', varargin{i});
                end
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
            fprintf('quadest.util.Config (v%s)\n', obj.version);
            fprintf('  Grid: nth=%d, nph=%d\n', obj.nth, obj.nph);
            fprintf('  Tolerance: %.2e\n', obj.tol);
            fprintf('  Upsampling factors: [%s]\n', num2str(obj.upsampFactors));
            fprintf('  Interpolate roots: %s\n', mat2str(obj.interpolateRoots));
        end
    end
end
