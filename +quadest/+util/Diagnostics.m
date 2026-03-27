classdef Diagnostics
    % DIAGNOSTICS Logging, warnings, and stability monitoring
    %
    % This class provides diagnostic utilities for monitoring numerical
    % stability and reporting issues during error estimation.
    %
    % Usage:
    %   quadest.util.Diagnostics.warn('Newton solver did not converge');
    %   quadest.util.Diagnostics.checkDensityResolution(density, nph);
    %
    % See also: quadest.util.Config
    
    methods (Static)
        function warn(msg, varargin)
            % WARN Issue a warning with quadest prefix
            %
            %   Diagnostics.warn('message') issues warning
            %   Diagnostics.warn('message with %d', value) with formatting
            
            if nargin > 1
                msg = sprintf(msg, varargin{:});
            end
            warning('quadest:Diagnostics', '%s', msg);
        end
        
        function info(msg, varargin)
            % INFO Print informational message
            %
            %   Diagnostics.info('message') prints to console
            
            if nargin > 1
                msg = sprintf(msg, varargin{:});
            end
            fprintf('[quadest] %s\n', msg);
        end
        
        function resolved = checkDensityResolution(density, nph, threshold)
            % CHECKDENSITYRESOLUTION Check if density is well-resolved via FFT decay
            %
            %   resolved = checkDensityResolution(density, nph) checks if
            %   the last 5 Fourier modes in phi direction decay below 1e-12.
            %
            %   resolved = checkDensityResolution(density, nph, threshold)
            %   uses custom threshold.
            %
            % Inputs:
            %   density   - [nth*nph × ncomp] density array
            %   nph       - Number of phi points
            %   threshold - (optional) Decay threshold, default 1e-12
            %
            % Outputs:
            %   resolved - true if density is well-resolved
            
            arguments
                density (:,:) {mustBeNumeric}
                nph (1,1) {mustBePositive, mustBeInteger}
                threshold (1,1) {mustBePositive} = 1e-12
            end
            
            nth = size(density, 1) / nph;
            ncomp = size(density, 2);
            ncheck = 5;  % Number of high-frequency modes to check
            
            maxErr = 0;
            for comp = 1:ncomp
                qmat = reshape(density(:, comp), nth, nph);
                for ii = 1:nth
                    qhat = fft(qmat(ii, :).') / (nph/2);
                    qhat = qhat(1:nph/2);  % Assume density is real-valued
                    qhatMax = max(abs(qhat));
                    if qhatMax > 0
                        err = max(abs(qhat(end-ncheck+1:end) / qhatMax));
                        maxErr = max(maxErr, err);
                    end
                end
            end
            
            resolved = maxErr <= threshold;
            
            if ~resolved
                quadest.util.Diagnostics.warn( ...
                    'Density may not be well-resolved. Max FFT decay error: %.2e (threshold: %.2e)', ...
                    maxErr, threshold);
            end
        end
        
        function ok = checkNewtonConvergence(iterations, maxIter, residual, tol)
            % CHECKNEWTONCONVERGENCE Check if Newton solver converged
            %
            %   ok = checkNewtonConvergence(iterations, maxIter, residual, tol)
            %
            % Returns true if converged, issues warning if not.
            
            ok = (iterations < maxIter) && (residual < tol);
            
            if ~ok
                quadest.util.Diagnostics.warn( ...
                    'Newton solver: %d iterations, residual %.2e (tol %.2e)', ...
                    iterations, residual, tol);
            end
        end
        
        function checkRootQuality(root, gamma, target, tol)
            % CHECKROOTQUALITY Verify that a computed root satisfies |gamma(root) - target|^2 ≈ 0
            %
            %   checkRootQuality(root, gamma, target, tol)
            %
            % Inputs:
            %   root   - Complex root value
            %   gamma  - Parameterization function handle
            %   target - Target point [1×3]
            %   tol    - Tolerance for residual check
            
            residual = sum(abs(gamma(root) - target).^2);
            if residual > tol
                quadest.util.Diagnostics.warn( ...
                    'Root quality check failed: |gamma(root) - target|^2 = %.2e', ...
                    residual);
            end
        end
        
        function nearAxis = checkAxisSingularity(targets, threshold)
            % CHECKAXISSINGULARITY Check for targets near the z-axis
            %
            %   nearAxis = checkAxisSingularity(targets, threshold)
            %
            % Targets with sqrt(x^2 + y^2) < threshold are flagged.
            
            arguments
                targets (:,3) {mustBeNumeric}
                threshold (1,1) {mustBePositive} = 1e-10
            end
            
            rxy = sqrt(targets(:,1).^2 + targets(:,2).^2);
            nearAxis = rxy < threshold;
            
            if any(nearAxis)
                quadest.util.Diagnostics.warn( ...
                    '%d target(s) near z-axis (r_xy < %.2e). Special handling may be needed.', ...
                    sum(nearAxis), threshold);
            end
        end
        
        function meta = getMetadata()
            % GETMETADATA Return package metadata for reproducibility
            %
            %   meta = getMetadata() returns struct with version, date, MATLAB info
            
            meta.version = quadest.util.Config.version;
            meta.timestamp = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
            meta.matlabVersion = version;
            meta.computer = computer;
        end
    end
end
