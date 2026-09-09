classdef Diagnostics
    % DIAGNOSTICS Logging, warnings, and stability monitoring
    %
    % This class provides diagnostic utilities for monitoring numerical
    % stability and reporting issues during error-indicator evaluation.
    %
    % Usage:
    %   quadind.util.Diagnostics.warn('Newton solver did not converge');
    %   quadind.util.Diagnostics.checkDensityResolution(density, nph);
    %   quadind.util.Diagnostics.checkGeometryResolution(geometry, nth);
    %
    % See also: quadind.util.Config
    
    methods (Static)
        function warn(msg, varargin)
            % WARN Issue a warning with quadind prefix
            %
            %   Diagnostics.warn('message') issues warning
            %   Diagnostics.warn('message with %d', value) with formatting
            
            if nargin > 1
                msg = sprintf(msg, varargin{:});
            end
            warning('quadind:Diagnostics', '%s', msg);
        end
        
        function info(msg, varargin)
            % INFO Print informational message
            %
            %   Diagnostics.info('message') prints to console
            
            if nargin > 1
                msg = sprintf(msg, varargin{:});
            end
            fprintf('[quadind] %s\n', msg);
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
                quadind.util.Diagnostics.warn( ...
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
                quadind.util.Diagnostics.warn( ...
                    'Newton solver: %d iterations, residual %.2e (tol %.2e)', ...
                    iterations, residual, tol);
            end
        end
        
        function [valid, residual] = checkInterpolatedThetaRoots(geometry, roots, phi, targets, tol, warnOnFailure)
            % CHECKINTERPOLATEDTHETAROOTS Check the analytic complex root equation
            %
            %   [valid, residual] = checkInterpolatedThetaRoots(geometry, roots, phi, targets, tol)
            %   [...] = checkInterpolatedThetaRoots(..., warnOnFailure)
            %
            % A complex singularity root satisfies the analytic equation
            % sum((gamma(theta,phi)-target).^2) = 0. The numerator must not
            % conjugate the complex displacement.

            arguments
                geometry (1,1) quadind.geometry.AxsymGeometry
                roots (:,1) {mustBeNumeric}
                phi (:,1) {mustBeNumeric}
                targets (:,3) {mustBeNumeric}
                tol (1,1) {mustBePositive}
                warnOnFailure (1,1) logical = true
            end

            if size(roots,1) ~= size(targets,1) || size(phi,1) ~= size(targets,1)
                error('quadind:Diagnostics:targetCountMismatch', ...
                    'roots, phi, and targets must have the same number of rows.');
            end

            displacement = geometry.evaluate(roots, phi) - targets;
            residualAbs = abs(sum(displacement.^2, 2));
            residualScale = max(sum(abs(displacement).^2, 2), realmin);
            residual = residualAbs ./ residualScale;
            valid = isfinite(roots) & residual <= tol;

            if warnOnFailure && any(~valid)
                warning('quadind:Diagnostics:interpolatedRootResidual', ...
                    ['%d interpolated theta root(s) failed the analytic root check ' ...
                     '(maximum normalized residual %.2e, tolerance %.2e).'], ...
                    sum(~valid), max(residual(~valid)), tol);
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
                quadind.util.Diagnostics.warn( ...
                    '%d target(s) near z-axis (r_xy < %.2e). Special handling may be needed.', ...
                    sum(nearAxis), threshold);
            end
        end
        
        function meta = getMetadata()
            % GETMETADATA Return package metadata for reproducibility
            %
            %   meta = getMetadata() returns struct with version, date, MATLAB info
            
            meta.version = quadind.util.Config.version;
            meta.timestamp = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
            meta.matlabVersion = version;
            meta.computer = computer;
        end
        
        function result = checkGeometryResolution(geometry, nth, threshold)
            % CHECKGEOMETRYRESOLUTION Check if geometry is well-resolved by GL grid
            %
            %   result = checkGeometryResolution(geometry, nth) checks whether
            %   nth Gauss-Legendre nodes adequately resolve the geometry.
            %
            %   result = checkGeometryResolution(geometry, nth, threshold)
            %   uses a custom relative error threshold (default: 1e-10).
            %
            % Two complementary checks are performed:
            %   1. Surface area convergence: compares sum(w) at nth vs 2*nth
            %      GL nodes. This directly tests whether the Jacobian (which
            %      enters quadrature weights) is adequately resolved.
            %   2. Pointwise interpolation of at(theta) and ct(theta) from
            %      nth to 2*nth GL nodes via barycentric Lagrange. This
            %      ensures individual surface points are accurately placed.
            %
            % Inputs:
            %   geometry  - AxsymGeometry object
            %   nth       - Number of theta (GL) points
            %   threshold - (optional) Relative error threshold, default 1e-10
            %
            % Outputs:
            %   result - struct with fields:
            %     resolved  - true if all checks pass
            %     maxError  - maximum error across all checks
            %     areaError - relative surface area convergence error
            %     atError   - relative interpolation error for at(theta)
            %     ctError   - relative interpolation error for ct(theta)
            
            arguments
                geometry (1,1) quadind.geometry.AxsymGeometry
                nth (1,1) {mustBePositive, mustBeInteger}
                threshold (1,1) {mustBePositive} = 1e-10
            end
            
            % GL nodes and weights on coarse and fine grids
            [theta_c, wtheta_c] = quadind.util.GaussLegendre(nth, 0, pi);
            [theta_f, wtheta_f] = quadind.util.GaussLegendre(2*nth, 0, pi);
            
            % --- Surface area convergence check ---
            % Compute Jacobian J(theta) = |at| * sqrt(dadt^2 + dcdt^2)
            at_c = geometry.at(theta_c);
            dadt_c = geometry.dadt(theta_c);
            dcdt_c = geometry.dcdt(theta_c);
            J_c = abs(at_c) .* sqrt(dadt_c.^2 + dcdt_c.^2);
            
            at_f = geometry.at(theta_f);
            dadt_f = geometry.dadt(theta_f);
            dcdt_f = geometry.dcdt(theta_f);
            J_f = abs(at_f) .* sqrt(dadt_f.^2 + dcdt_f.^2);
            
            % Surface area = 2*pi * integral of J over theta
            area_c = 2 * pi * sum(wtheta_c .* J_c);
            area_f = 2 * pi * sum(wtheta_f .* J_f);
            
            if area_f > 0
                areaErr = abs(area_c - area_f) / area_f;
            else
                areaErr = abs(area_c - area_f);
            end
            
            % --- Pointwise interpolation check for at and ct ---
            ct_c = geometry.ct(theta_c);
            ct_f = geometry.ct(theta_f);
            
            at_interp = quadind.util.Diagnostics.bclagInterp(theta_c, at_c, theta_f);
            ct_interp = quadind.util.Diagnostics.bclagInterp(theta_c, ct_c, theta_f);
            
            atErr = quadind.util.Diagnostics.relativeError(at_interp, at_f);
            ctErr = quadind.util.Diagnostics.relativeError(ct_interp, ct_f);
            
            % --- Aggregate results ---
            maxErr = max([areaErr, atErr, ctErr]);
            resolved = maxErr <= threshold;
            
            result.resolved = resolved;
            result.maxError = maxErr;
            result.areaError = areaErr;
            result.atError = atErr;
            result.ctError = ctErr;
            
            if ~resolved
                names = {'surface area', 'at(theta)', 'ct(theta)'};
                errs = [areaErr, atErr, ctErr];
                [~, idx] = max(errs);
                quadind.util.Diagnostics.warn( ...
                    ['Geometry may not be well-resolved with nth = %d. ' ...
                     'Max error in %s: %.2e (threshold: %.2e). ' ...
                     'Consider increasing nth.'], ...
                    nth, names{idx}, maxErr, threshold);
            end
        end
    end
    
    methods (Static, Access = private)
        function f_tgt = bclagInterp(x_src, f_src, x_tgt)
            % BCLAGINTERP Barycentric Lagrange interpolation
            %
            % Interpolates data f_src at nodes x_src to target points x_tgt.
            % Uses the Berrut-Trefethen second-form barycentric formula.
            %
            % Reference: Berrut & Trefethen, SIAM Review 46(3), 2004.
            
            n = numel(x_src);
            N = numel(x_tgt);
            
            % Compute barycentric weights
            w = zeros(n, 1);
            for j = 1:n
                w(j) = 1 / prod(x_src(j) - x_src([1:j-1, j+1:n]));
            end
            
            % Evaluate interpolant at target points
            numer = zeros(N, 1);
            denom = zeros(N, 1);
            exact = zeros(N, 1);
            
            for j = 1:n
                xdiff = x_tgt - x_src(j);
                temp = w(j) ./ xdiff;
                numer = numer + temp .* f_src(j);
                denom = denom + temp;
                exact(xdiff == 0) = j;
            end
            
            f_tgt = numer ./ denom;
            
            % Handle exact matches
            jj = find(exact);
            f_tgt(jj) = f_src(exact(jj));
        end
        
        function err = relativeError(approx, exact)
            % RELATIVEERROR Max relative error with safe denominator
            %
            % Uses max(|exact|) as the normalization to avoid dividing
            % by zero at individual nodes (e.g., at poles where at=0).
            
            scale = max(abs(exact));
            if scale == 0
                err = max(abs(approx - exact));
            else
                err = max(abs(approx - exact)) / scale;
            end
        end
    end
end
