classdef UniformEstimateBuilder
    % UNIFORMESTIMATEBUILDER Precompute uniform error estimates on tabulation grid
    %
    % Builds gridded interpolants for error estimates assuming unit density.
    % The actual error estimate is obtained by multiplying by the density
    % magnitude at complex roots (see DensityModifier).
    %
    % This is the most computationally expensive part of the precomputation
    % and should only be done once per geometry/discretization.
    %
    % See also: quadest.errorest.ErrorEstimator, quadest.errorest.DensityModifier
    
    properties (Constant, Access = private)
        % Singularity order for stresslet kernel
        STRESSLET_P = 5/2
    end
    
    methods (Static)
        function [interpolants, rootInterp, tabGrid] = build(geometry, grid, kernel, upsampFactors, interpolateRoots, config)
            % BUILD Create gridded interpolants for error estimates
            %
            %   [interpolants, rootInterp, tabGrid] = build(geometry, grid, kernel, upsampFactors, interpolateRoots, config)
            %
            % Inputs:
            %   geometry         - AxsymGeometry object
            %   grid             - AxsymGrid object (base discretization)
            %   kernel           - Kernel object
            %   upsampFactors    - Array of upsampling factors [1, 2, 3, ...]
            %   interpolateRoots - Whether to also build theta root interpolants
            %   config           - Config object
            %
            % Outputs:
            %   interpolants - Cell array {nfac × ncomp} of griddedInterpolant objects
            %                  for log10(error estimate)
            %   rootInterp   - Struct with theta root interpolants (if interpolateRoots)
            %   tabGrid      - Struct with tabulation grid info (for debugging)
            
            arguments
                geometry (1,1) quadest.geometry.AxsymGeometry
                grid (1,1) quadest.grid.AxsymGrid
                kernel (1,1) quadest.kernel.Kernel
                upsampFactors (1,:) {mustBePositive, mustBeInteger}
                interpolateRoots (1,1) logical = true
                config (1,1) quadest.util.Config = quadest.util.Config()
            end
            
            quadest.util.Diagnostics.info('Building uniform error estimate interpolants...');
            
            % Ensure upsampling factors include 1 (direct quadrature)
            upsampFactors = unique([1, upsampFactors]);
            nfac = length(upsampFactors);
            ncomp = kernel.numComponents();
            p = kernel.singularityOrder();
            
            % Build tabulation grid
            [tabGrid, evalPoints, evalMask] = quadest.errorest.UniformEstimateBuilder.buildTabGrid(geometry, config);
            
            % Initialize storage
            interpolants = cell(max(upsampFactors), ncomp);
            if interpolateRoots
                rootInterp.Re = cell(max(upsampFactors), 1);
                rootInterp.Im = cell(max(upsampFactors), 1);
            else
                rootInterp = [];
            end
            
            % Find axis points that need special handling
            axisPoints = quadest.errorest.UniformEstimateBuilder.findAxisPoints(evalPoints, evalMask, geometry);
            
            % Loop through upsampling factors
            for j = 1:nfac
                upfac = upsampFactors(j);
                quadest.util.Diagnostics.info('  Upsampling factor %d/%d (kappa=%d)', j, nfac, upfac);
                
                % Create upsampled grid
                if upfac == 1
                    refinedGrid = grid;
                else
                    refinedGrid = grid.upsample(upfac);
                end
                
                % Compute error estimates at exterior tabulation points
                [estimates, thetaRoots] = quadest.errorest.UniformEstimateBuilder.computeEstimates(...
                    geometry, refinedGrid, kernel, evalPoints(evalMask, :), p);
                
                % Build full estimate array (interior points get large values)
                fullEstimates = ones(size(evalPoints, 1), ncomp) * config.interiorEstimate;
                fullEstimates(evalMask, :) = estimates;
                
                % Handle NaN values
                fullEstimates(isnan(fullEstimates)) = config.interiorEstimate;
                
                % Safety: axis point closest to surface
                if ~isempty(axisPoints)
                    fullEstimates(axisPoints, :) = config.interiorEstimate;
                end
                
                % Add lower bound
                fullEstimates = fullEstimates + config.estimateLowerBound;
                
                % Build interpolants for each component
                for c = 1:ncomp
                    estMat = reshape(fullEstimates(:, c), config.ntab, config.nztab);
                    interpolants{upfac, c} = griddedInterpolant(...
                        tabGrid.xy_grid, tabGrid.z_grid, log10(estMat), 'linear');
                end
                
                % Store theta root interpolants if requested
                if interpolateRoots
                    fullRoots = zeros(size(evalPoints, 1), 1);
                    fullRoots(evalMask) = thetaRoots;
                    rootMat = reshape(fullRoots, config.ntab, config.nztab);
                    
                    rootInterp.Re{upfac} = griddedInterpolant(...
                        tabGrid.xy_grid, tabGrid.z_grid, real(rootMat), 'linear');
                    rootInterp.Im{upfac} = griddedInterpolant(...
                        tabGrid.xy_grid, tabGrid.z_grid, abs(imag(rootMat)), 'linear');
                end
            end
            
            quadest.util.Diagnostics.info('Done building interpolants.');
        end
        
        function [estimates, thetaRoots] = computeEstimates(geometry, grid, kernel, targets, p)
            % COMPUTEESTIMATES Compute uniform error estimates at target points
            %
            %   [estimates, thetaRoots] = computeEstimates(geometry, grid, kernel, targets, p)
            %
            % Computes per-component error estimates assuming unit density e_c
            % for each component c. This is the expensive per-target computation
            % that the tabulated interpolant replaces.
            %
            % Inputs:
            %   geometry - AxsymGeometry object
            %   grid     - AxsymGrid object
            %   kernel   - Kernel object
            %   targets  - [M×3] target points (must be exterior)
            %   p        - Singularity order (kernel.singularityOrder())
            %
            % Outputs:
            %   estimates  - [M×ncomp] per-component error estimates
            %   thetaRoots - [M×1] complex theta roots
            
            M = size(targets, 1);
            ncomp = kernel.numComponents();
            
            estimates = zeros(M, ncomp);
            thetaRoots = zeros(M, 1);
            
            % Find closest grid points
            idx_closest = knnsearch(grid.x, targets);
            
            % Get Gauss-Laguerre quadrature
            [xlag, wlag] = quadest.util.GaussLaguerre8();
            
            % Loop through target points
            for ii = 1:M
                target = targets(ii, :);
                
                [est, theta0] = quadest.errorest.UniformEstimateBuilder.computeSingleEstimate(...
                    geometry, grid, target, idx_closest(ii), p, xlag, wlag);
                
                estimates(ii, :) = est;
                thetaRoots(ii) = theta0;
            end
        end
        
        function estimates = computeDirectEstimates(geometry, grid, kernel, targets, density)
            % COMPUTEDIRECTESTIMATES Compute error estimates with given density
            %
            %   estimates = computeDirectEstimates(geometry, grid, kernel, targets, density)
            %
            % Computes error estimates directly for each target point using the
            % actual density, without precomputed tabulation. This is the expensive
            % per-target computation; use for benchmarking against the fast
            % tabulated method (ErrorEstimator.evaluate).
            %
            % Inputs:
            %   geometry - AxsymGeometry object
            %   grid     - AxsymGrid object
            %   kernel   - Kernel object
            %   targets  - [M×3] target points (must be exterior)
            %   density  - [N×ncomp] density on grid (N = nth*nph)
            %
            % Outputs:
            %   estimates - [M×1] scalar error estimates
            
            M = size(targets, 1);
            ncomp = kernel.numComponents();
            p = kernel.singularityOrder();
            
            % Reshape density for interpolation: [nth × nph × ncomp]
            dens_3d = reshape(density, grid.nth, grid.nph, ncomp);
            
            estimates = zeros(M, 1);
            
            % Find closest grid points
            idx_closest = knnsearch(grid.x, targets);
            
            % Get Gauss-Laguerre quadrature
            [xlag, wlag] = quadest.util.GaussLaguerre8();
            
            % Loop through target points
            for ii = 1:M
                target = targets(ii, :);
                
                estimates(ii) = quadest.errorest.UniformEstimateBuilder.computeSingleDirectEstimate(...
                    geometry, grid, target, idx_closest(ii), p, xlag, wlag, dens_3d);
            end
        end
    end
    
    methods (Static, Access = private)
        function [tabGrid, evalPoints, evalMask] = buildTabGrid(geometry, config)
            % BUILDTABGRID Create tabulation grid for interpolants
            
            % Grid extends to 2x the geometry extent
            maxExtent = 2 * (geometry.maxRadius() + geometry.maxHeight());
            
            % XY coordinates (radial distance from z-axis)
            xy = linspace(0, maxExtent, config.ntab - 1);
            xy = [-xy(2), xy];  % Add point for axis singularity handling
            
            % Z coordinates
            z = linspace(0, maxExtent, config.nztab - 1);
            z = [-z(2), z];
            
            [xy_grid, z_grid] = ndgrid(xy, z);
            
            % Create 3D evaluation points (y = 0 slice)
            [X, Y, Z] = ndgrid(xy, 0, z);
            evalPoints = [X(:), Y(:), Z(:)];
            
            % Find exterior points
            evalMask = geometry.isExterior(evalPoints);
            
            tabGrid.xy = xy;
            tabGrid.z = z;
            tabGrid.xy_grid = xy_grid;
            tabGrid.z_grid = z_grid;
        end
        
        function axisPoints = findAxisPoints(evalPoints, evalMask, geometry)
            % FINDAXISPOINTS Find axis points that need special handling
            
            % Points on z-axis (xy = 0) that are exterior
            onAxis = evalPoints(:, 1) == 0 & evalPoints(:, 2) == 0 & evalMask;
            
            if any(onAxis)
                axisIdx = find(onAxis);
                % Find the one closest to the surface
                [~, minIdx] = min(abs(evalPoints(axisIdx, 3)));
                axisPoints = axisIdx(minIdx);
            else
                axisPoints = [];
            end
        end
        
        function [est, theta0] = computeSingleEstimate(geometry, grid, target, idx_R, p, xlag, wlag)
            % COMPUTESINGLEESTIMATE Compute error estimate for single target
                        
            % Get closest grid point indices
            [itheta, iphi] = grid.ind2sub(idx_R);
            theta_star = grid.theta(itheta);
            phi_star = grid.phi(iphi);
            
            % Linear map for theta
            t_star = (2/pi) * theta_star - 1;
            imap = @(t) (pi/2) * (t + 1);
            dfac = pi/2;
            
            % Compute phi root
            [phi0, Gphi] = quadest.errorest.RootFinder.computePhiRoot(geometry, target, theta_star);
            
            % Compute theta root
            if target(3) ~= 0 && isa(geometry, 'quadest.geometry.Spheroid')
                % Use analytic formula for spheroid
                theta0 = geometry.findThetaRoot(target, phi_star);
            else
                % Use Newton solver
                gamma_t = @(t) geometry.evaluate(imap(t), phi_star);
                dgamma_t = @(t) dfac * geometry.drdtheta(imap(t), phi_star);
                [t0, ~] = quadest.errorest.RootFinder.newtonSolve(gamma_t, dgamma_t, target, t_star);
                theta0 = imap(t0);
            end
            t0 = (2/pi) * theta0 - 1;
            
            % Geometric factor for theta
            % Note: Use dot() to match legacy behavior with complex vectors
            gamma_t = @(t) geometry.evaluate(imap(t), phi_star);
            dgamma_t = @(t) dfac * geometry.drdtheta(imap(t), phi_star);
            r_theta = gamma_t(t0) - target;
            Gtheta = 2 * dot(r_theta, dgamma_t(t0));
            
            % Surface vectors at closest point
            srfvec = geometry.evaluate(theta_star, phi_star);
            rvec = srfvec - target;
            
            % Semi-analytical root functions
            drdphi = geometry.drdphi(theta_star, phi_star);
            drdtheta = geometry.drdtheta(theta_star, phi_star);
            dp0fun = quadest.errorest.RootFinder.semiRoot(rvec, drdphi, drdtheta);
            dt0fun = quadest.errorest.RootFinder.semiRoot(rvec, drdtheta, drdphi);
            
            % Scaling constants
            C_TZ = grid.nph * norm(drdtheta) / norm(drdphi);
            C_GL = max(grid.nth * norm(drdphi) / norm(drdtheta), grid.nth);
            
            % Compute estimates for phi direction (trapezoidal rule)
            estTZ = quadest.errorest.UniformEstimateBuilder.computePhiEstimate(...
                geometry, grid, target, theta_star, phi0, Gphi, p, ...
                dp0fun, C_TZ, xlag, wlag);
            
            % Compute estimates for theta direction (Gauss-Legendre)
            estGL = quadest.errorest.UniformEstimateBuilder.computeThetaEstimate(...
                geometry, grid, target, phi_star, t0, Gtheta, p, ...
                dt0fun, C_GL, xlag, wlag, imap, dfac);
            
            % Combined estimate (sum of both directions)
            est = estTZ + estGL;
        end
        
        function est = computeSingleDirectEstimate(geometry, grid, target, idx_R, p, xlag, wlag, dens_3d)
            % COMPUTESINGLEDIRECTESTIMATE Compute error estimate for single target with density
            %
            % Like computeSingleEstimate but uses the actual density (interpolated
            % at complex roots) instead of unit density vectors. Returns a scalar.
            
            % Get closest grid point indices
            [itheta, iphi] = grid.ind2sub(idx_R);
            theta_star = grid.theta(itheta);
            phi_star = grid.phi(iphi);
            
            % Linear map for theta
            t_star = (2/pi) * theta_star - 1;
            imap = @(t) (pi/2) * (t + 1);
            dfac = pi/2;
            
            % Compute phi root
            [phi0, Gphi] = quadest.errorest.RootFinder.computePhiRoot(geometry, target, theta_star);
            
            % Compute theta root
            if target(3) ~= 0 && isa(geometry, 'quadest.geometry.Spheroid')
                theta0 = geometry.findThetaRoot(target, phi_star);
            else
                gamma_t = @(t) geometry.evaluate(imap(t), phi_star);
                dgamma_t = @(t) dfac * geometry.drdtheta(imap(t), phi_star);
                [t0, ~] = quadest.errorest.RootFinder.newtonSolve(gamma_t, dgamma_t, target, t_star);
                theta0 = imap(t0);
            end
            t0 = (2/pi) * theta0 - 1;
            
            % Geometric factor for theta
            gamma_t = @(t) geometry.evaluate(imap(t), phi_star);
            dgamma_t = @(t) dfac * geometry.drdtheta(imap(t), phi_star);
            r_theta = gamma_t(t0) - target;
            Gtheta = 2 * dot(r_theta, dgamma_t(t0));
            
            % Surface vectors at closest point
            srfvec = geometry.evaluate(theta_star, phi_star);
            rvec = srfvec - target;
            
            % Semi-analytical root functions
            drdphi = geometry.drdphi(theta_star, phi_star);
            drdtheta = geometry.drdtheta(theta_star, phi_star);
            dp0fun = quadest.errorest.RootFinder.semiRoot(rvec, drdphi, drdtheta);
            dt0fun = quadest.errorest.RootFinder.semiRoot(rvec, drdtheta, drdphi);
            
            % Scaling constants
            C_TZ = grid.nph * norm(drdtheta) / norm(drdphi);
            C_GL = max(grid.nth * norm(drdphi) / norm(drdtheta), grid.nth);
            
            % Interpolate density at phi root (2-point linear in phi direction)
            q_phi = quadest.errorest.UniformEstimateBuilder.interpDensityPhiSingle(...
                grid, dens_3d, phi0, itheta);
            
            % Interpolate density at theta root (2-point linear in theta direction)
            q_theta = quadest.errorest.UniformEstimateBuilder.interpDensityThetaSingle(...
                grid, dens_3d, theta0, theta_star, iphi);
            
            % Compute estimates for phi direction with actual density
            estTZ = quadest.errorest.UniformEstimateBuilder.computePhiEstimateDirect(...
                geometry, grid, target, theta_star, phi0, Gphi, p, ...
                dp0fun, C_TZ, xlag, wlag, q_phi);
            
            % Compute estimates for theta direction with actual density
            estGL = quadest.errorest.UniformEstimateBuilder.computeThetaEstimateDirect(...
                geometry, grid, target, phi_star, t0, Gtheta, p, ...
                dt0fun, C_GL, xlag, wlag, imap, dfac, q_theta);
            
            % Combined estimate (sum of both directions)
            est = estTZ + estGL;
        end
        
        function q_interp = interpDensityPhiSingle(grid, dens_3d, phi0, itheta)
            % INTERPDENSITYPHISINGLE 2-point linear interpolation in phi for single target
            
            ncomp = size(dens_3d, 3);
            phi_grid = grid.phi(:);
            nph = grid.nph;
            phi_val = real(phi0);
            
            % Find two nearest phi points (accounting for periodicity)
            [~, sorted_idx] = sort(abs(phi_grid - phi_val));
            
            if sorted_idx(1) < 2
                tmpphi = phi_grid;
                tmpphi(end) = tmpphi(end) - 2*pi;
                [~, sorted_idx] = sort(abs(tmpphi - phi_val));
            elseif sorted_idx(1) > nph - 1
                tmpphi = phi_grid;
                tmpphi(1) = tmpphi(1) + 2*pi;
                [~, sorted_idx] = sort(abs(tmpphi - phi_val));
            end
            
            local_ip = sort(sorted_idx(1:2));
            p_pts = phi_grid(local_ip);
            
            % Handle wrapped values
            if abs(p_pts(2) - p_pts(1)) > pi
                if p_pts(1) < pi
                    p_pts(2) = p_pts(2) - 2*pi;
                else
                    p_pts(1) = p_pts(1) - 2*pi;
                end
            end
            
            q_interp = zeros(1, ncomp);
            for c = 1:ncomp
                q_vals = squeeze(dens_3d(itheta, local_ip, c));
                q_interp(c) = (q_vals(1) * (p_pts(2) - phi0) + ...
                               q_vals(2) * (phi0 - p_pts(1))) / ...
                              (p_pts(2) - p_pts(1));
            end
        end
        
        function q_interp = interpDensityThetaSingle(grid, dens_3d, theta0, theta_star, iphi)
            % INTERPDENSITYTHETASINGLE 2-point linear interpolation in theta for single target
            
            ncomp = size(dens_3d, 3);
            theta_grid = grid.theta(:);
            
            % Find two nearest theta points
            [~, sorted_idx] = mink(abs(theta_grid - theta_star), 2);
            local_it = sort(sorted_idx);
            t_pts = theta_grid(local_it);
            
            q_interp = zeros(1, ncomp);
            denom = t_pts(2) - t_pts(1);
            for c = 1:ncomp
                q_vals = squeeze(dens_3d(local_it, iphi, c));
                if abs(denom) < eps
                    q_interp(c) = q_vals(1);
                else
                    q_interp(c) = (q_vals(1) * (t_pts(2) - theta0) + ...
                                   q_vals(2) * (theta0 - t_pts(1))) / denom;
                end
            end
            
            % Handle NaN (can occur for theta0 at poles)
            q_interp(isnan(q_interp)) = 0;
        end
        
        function estTZ = computePhiEstimateDirect(geometry, grid, target, theta_star, phi0, Gphi, p, dp0fun, C_TZ, xlag, wlag, q_density)
            % COMPUTEPHIESTIMATEDIRECT Error estimate from phi direction with given density
            %
            % Like computePhiEstimate but uses the interpolated density q_density
            % instead of unit vectors. Returns a scalar.
            
            % Check both +/- imaginary parts
            phi0_candidates = [phi0; conj(phi0)];
            estTZ_both = zeros(2, 1);
            
            for ii = 1:2
                phi0_test = phi0_candidates(ii);
                
                % Kernel at phi root with actual density
                fp0t = quadest.errorest.UniformEstimateBuilder.kernelAtPhiRoot(...
                    geometry, target, theta_star, phi0_test, q_density);
                
                % Integration
                est_const = abs(2 / gamma(p) * fp0t * Gphi^(-p));
                ds = xlag / C_TZ;
                
                dtpos = phi0_test + dp0fun(ds) - dp0fun(0);
                dtneg = phi0_test + dp0fun(-ds) - dp0fun(0);
                
                % Trapezoidal error function
                int_pos = quadest.errorest.UniformEstimateBuilder.trapzErrFunc(dtpos, grid.nph, p-1);
                int_neg = quadest.errorest.UniformEstimateBuilder.trapzErrFunc(dtneg, grid.nph, p-1);
                
                estTZ_both(ii) = est_const * sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_TZ);
            end
            
            % Take minimum over +/- imaginary part
            estTZ = min(estTZ_both);
            if isnan(estTZ); estTZ = 0; end
        end
        
        function estGL = computeThetaEstimateDirect(geometry, grid, target, phi_star, t0, Gtheta, p, dt0fun, C_GL, xlag, wlag, imap, dfac, q_density)
            % COMPUTETHETAESTIMATEDIRECT Error estimate from theta direction with given density
            %
            % Like computeThetaEstimate but uses the interpolated density q_density
            % instead of unit vectors. Returns a scalar.
            
            % Check both +/- imaginary parts
            t0_candidates = [t0; conj(t0)];
            estGL_both = zeros(2, 1);
            
            for ii = 1:2
                t0_test = t0_candidates(ii);
                
                % Kernel at theta root with actual density
                ft0p = quadest.errorest.UniformEstimateBuilder.kernelAtThetaRoot(...
                    geometry, target, t0_test, phi_star, q_density, imap, dfac);
                
                % Integration
                est_const = abs(2 / gamma(p) * ft0p * Gtheta^(-p));
                ds = xlag / C_GL;
                
                dtpos = t0_test - (dt0fun(ds) - dt0fun(0));
                dtneg = t0_test - (dt0fun(-ds) - dt0fun(0));
                
                % Gauss-Legendre error function
                int_pos = quadest.errorest.UniformEstimateBuilder.glErrFunc(dtpos, grid.nth, p-1);
                int_neg = quadest.errorest.UniformEstimateBuilder.glErrFunc(dtneg, grid.nth, p-1);
                
                estGL_both(ii) = est_const * sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_GL);
            end
            
            % Take minimum over +/- imaginary part
            estGL = min(estGL_both);
        end
        
        function estTZ = computePhiEstimate(geometry, grid, target, theta_star, phi0, Gphi, p, dp0fun, C_TZ, xlag, wlag)
            % COMPUTEPHIESTIMATE Error estimate from phi direction (trapezoidal)
            
            ncomp = 3;
            
            % Check both +/- imaginary parts
            phi0_candidates = [phi0; conj(phi0)];
            estTZ_both = zeros(2, ncomp);
            
            for ii = 1:2
                phi0_test = phi0_candidates(ii);
                
                for c = 1:ncomp
                    % Unit density for precomputation
                    density = zeros(1, 3);
                    density(c) = 1.0;
                    
                    % Kernel at phi root
                    fp0t = quadest.errorest.UniformEstimateBuilder.kernelAtPhiRoot(...
                        geometry, target, theta_star, phi0_test, density);
                    
                    % Integration
                    est_const = abs(2 / gamma(p) * fp0t * Gphi^(-p));
                    ds = xlag / C_TZ;
                    
                    dtpos = phi0_test + dp0fun(ds) - dp0fun(0);
                    dtneg = phi0_test + dp0fun(-ds) - dp0fun(0);
                    
                    % Trapezoidal error function
                    int_pos = quadest.errorest.UniformEstimateBuilder.trapzErrFunc(dtpos, grid.nph, p-1);
                    int_neg = quadest.errorest.UniformEstimateBuilder.trapzErrFunc(dtneg, grid.nph, p-1);
                    
                    estTZ_both(ii, c) = est_const * sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_TZ);
                end
            end
            
            % Take minimum over +/- imaginary part
            estTZ = min(estTZ_both, [], 1);
            estTZ(isnan(estTZ)) = 0;  % Handle axis case
        end
        
        function estGL = computeThetaEstimate(geometry, grid, target, phi_star, t0, Gtheta, p, dt0fun, C_GL, xlag, wlag, imap, dfac)
            % COMPUTETHETAESTIMATE Error estimate from theta direction (Gauss-Legendre)
            
            ncomp = 3;
            
            % Check both +/- imaginary parts
            t0_candidates = [t0; conj(t0)];
            estGL_both = zeros(2, ncomp);
            
            for ii = 1:2
                t0_test = t0_candidates(ii);
                
                for c = 1:ncomp
                    % Unit density for precomputation
                    density = zeros(1, 3);
                    density(c) = 1.0;
                    
                    % Kernel at theta root
                    ft0p = quadest.errorest.UniformEstimateBuilder.kernelAtThetaRoot(...
                        geometry, target, t0_test, phi_star, density, imap, dfac);
                    
                    % Integration
                    est_const = abs(2 / gamma(p) * ft0p * Gtheta^(-p));
                    ds = xlag / C_GL;
                    
                    dtpos = t0_test - (dt0fun(ds) - dt0fun(0));
                    dtneg = t0_test - (dt0fun(-ds) - dt0fun(0));
                    
                    % Gauss-Legendre error function
                    int_pos = quadest.errorest.UniformEstimateBuilder.glErrFunc(dtpos, grid.nth, p-1);
                    int_neg = quadest.errorest.UniformEstimateBuilder.glErrFunc(dtneg, grid.nth, p-1);
                    
                    estGL_both(ii, c) = est_const * sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_GL);
                end
            end
            
            % Take minimum over +/- imaginary part
            estGL = min(estGL_both, [], 1);
        end
        
        function fp0t = kernelAtPhiRoot(geometry, target, theta_star, phi0, density)
            % KERNELATPHIROOT Evaluate stresslet kernel at phi root
            %
            % Uses the unnormalized normal n = cross(dr/dtheta, dr/dphi), matching
            % the legacy kernel_phi0 which captures the 2nd output of element_normal
            % (the raw cross product, NOT the unit normal).
            
            srfvec = geometry.evaluate(theta_star, phi0);
            
            % Unnormalized normal: n = cross(drdtheta, drdphi)
            % For axisymmetric body: n = [-dcdt*at*cos(phi), -dcdt*at*sin(phi), dadt*at]
            at_val = geometry.at(theta_star);
            dadt_val = geometry.dadt(theta_star);
            dcdt_val = geometry.dcdt(theta_star);
            
            nvec = [-dcdt_val .* at_val .* cos(phi0), ...
                    -dcdt_val .* at_val .* sin(phi0), ...
                    dadt_val .* at_val];
            
            r = srfvec - target;
            rq = sum(r .* density);
            rn = sum(r .* nvec);
            f = r * (rq * rn);
            
            fp0t = -6 * max(abs(f));
        end
        
        function ft0p = kernelAtThetaRoot(geometry, target, t0, phi_star, density, imap, dfac)
            % KERNELATTHETAROOT Evaluate stresslet kernel at theta root
            %
            % Uses the unnormalized normal n = cross(drdtheta_mapped, drdphi) which
            % includes the dfac factor from the t->theta mapping chain rule, matching
            % the legacy kernel_t0 which uses element_normal(t0, phi, a, da, dc, imap, dfac).
            %
            % legacy drdtheta_mapped = dfac * [dadt(theta0)*cos(phi), dadt(theta0)*sin(phi), dcdt(theta0)]
            % so n = dfac * cross([dadt*cos,dadt*sin,dcdt], [-at*sin,at*cos,0])
            %      = dfac * [-dcdt*at*cos(phi), -dcdt*at*sin(phi), dadt*at]
            
            theta0 = imap(t0);
            srfvec = geometry.evaluate(theta0, phi_star);
            
            at_val = geometry.at(theta0);
            dadt_val = geometry.dadt(theta0);
            dcdt_val = geometry.dcdt(theta0);
            
            % Unnormalized normal including dfac (from t-to-theta mapping)
            nvec = dfac * [-dcdt_val .* at_val .* cos(phi_star), ...
                           -dcdt_val .* at_val .* sin(phi_star), ...
                           dadt_val .* at_val];
            
            r = srfvec - target;
            rq = sum(r .* density);
            rn = sum(r .* nvec);
            f = r * (rq * rn);
            
            ft0p = -6 * max(abs(f));
        end
        
        function knq = trapzErrFunc(z, n, q)
            % TRAPZERRFUNC q-th derivative of trapezoidal error function
            b = abs(imag(z));
            knq = 2 * pi * n^q * exp(-n * b);
        end
        
        function knq = glErrFunc(z, n, q)
            % GLERRFUNC q-th derivative of Gauss-Legendre error function
            
            % Transform to first quadrant
            z = abs(real(z)) + 1i * imag(z);
            
            knq = 2 * pi ./ (z + sqrt(z.^2 - 1)).^(2*n + 1);
            if q ~= 0
                knq = knq .* (-(2*n + 1) ./ sqrt(z.^2 - 1)).^q;
            end
            
            knq = abs(knq);
        end
    end
end
