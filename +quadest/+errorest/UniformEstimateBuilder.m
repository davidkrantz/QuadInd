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
            
            % Find closest grid points using parametric approach
            [itheta_all, iphi_all] = quadest.errorest.UniformEstimateBuilder.findNearestNodes(grid, targets);
            
            % Get Gauss-Laguerre quadrature
            [xlag, wlag] = quadest.util.GaussLaguerre8();
            
            % Loop through target points
            for ii = 1:M
                target = targets(ii, :);
                
                [est, theta0] = quadest.errorest.UniformEstimateBuilder.computeSingleEstimate(...
                    geometry, grid, target, itheta_all(ii), iphi_all(ii), p, xlag, wlag);
                
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
            
            % Find closest grid points using parametric approach
            [itheta_all, iphi_all] = quadest.errorest.UniformEstimateBuilder.findNearestNodes(grid, targets);
            
            % Get Gauss-Laguerre quadrature
            [xlag, wlag] = quadest.util.GaussLaguerre8();
            
            % Loop through target points
            for ii = 1:M
                target = targets(ii, :);
                
                estimates(ii) = quadest.errorest.UniformEstimateBuilder.computeSingleDirectEstimate(...
                    geometry, grid, target, itheta_all(ii), iphi_all(ii), p, xlag, wlag, dens_3d);
            end
        end

        function [itheta, iphi] = findNearestNodes(grid, targets)
            % FINDNEARESTNODES Find closest grid node for each target
            %
            % Uses the same parametric approach as DensityModifier:
            %   - phi: nearest index via rounding (exact for uniform grid)
            %   - theta: minimize actual 3D distance within nearest phi column
            %
            % This replaces knnsearch(grid.x, targets) and avoids the
            % Statistics Toolbox dependency and O(N) KD-tree overhead.
            
            % Phi (uniform grid): nearest index via rounding
            phi_target = mod(atan2(targets(:,2), targets(:,1)), 2*pi);
            dphi = 2*pi / grid.nph;
            iphi = mod(round(phi_target / dphi), grid.nph) + 1;
            
            % Theta (GL nodes): minimize 3D distance within nearest phi column.
            % d²(θ_i) = at(θ_i)² + rxy² − 2·at(θ_i)·rxy·cos(φ_col − φ_t) + (ct(θ_i) − z)²
            at_vals = grid.geometry.at(grid.theta(:));    % [nth × 1]
            ct_vals = grid.geometry.ct(grid.theta(:));    % [nth × 1]
            rxy = sqrt(targets(:,1).^2 + targets(:,2).^2);
            cos_dphi = cos(grid.phi(iphi(:)).' - phi_target);
            d2 = at_vals.'.^2 + rxy.^2 ...
                - 2 * (at_vals.' .* rxy) .* cos_dphi ...
                + (ct_vals.' - targets(:,3)).^2;          % [M × nth]
            [~, itheta] = min(d2, [], 2);
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
        
        function [est, theta0] = computeSingleEstimate(geometry, grid, target, itheta, iphi, p, xlag, wlag)
            % COMPUTESINGLEESTIMATE Compute error estimate for single target
            %
            % Optimized: closures eliminated, geometry evaluations minimized,
            % kernel and integration folded inline, component loop vectorized.
            
            % --- Grid values and mapping constants ---
            theta_star = grid.theta(itheta);
            phi_star = grid.phi(iphi);
            dfac = pi/2;
            
            % --- Phi root (analytic) ---
            [phi0, Gphi] = quadest.errorest.RootFinder.computePhiRoot(geometry, target, theta_star);
            
            % --- Theta root ---
            if target(3) ~= 0 && isa(geometry, 'quadest.geometry.Spheroid')
                theta0 = geometry.findThetaRoot(target, phi_star);
            else
                t_star = (2/pi) * theta_star - 1;
                [t0, ~] = quadest.errorest.RootFinder.newtonSolve(...
                    geometry, target, phi_star, t_star, dfac);
                theta0 = dfac * (t0 + 1);
            end
            t0 = (2/pi) * theta0 - 1;
            
            % --- Geometric factor Gtheta (one evaluate+drdtheta call) ---
            gamma_final = geometry.evaluate(theta0, phi_star);
            dgamma_final = dfac * geometry.drdtheta(theta0, phi_star);
            r_theta = gamma_final - target;
            Gtheta = 2 * dot(r_theta, dgamma_final);
            
            % --- Geometry at closest point (shared by kernel + semiRoot) ---
            at_star = geometry.at(theta_star);
            ct_star = geometry.ct(theta_star);
            dadt_star = geometry.dadt(theta_star);
            dcdt_star = geometry.dcdt(theta_star);
            
            rvec = [at_star * cos(phi_star), at_star * sin(phi_star), ct_star] - target;
            drdphi_vec = [-at_star * sin(phi_star), at_star * cos(phi_star), 0];
            drdtheta_vec = [dadt_star * cos(phi_star), dadt_star * sin(phi_star), dcdt_star];
            
            % --- SemiRoot coefficients (avoids closures) ---
            % Phi direction: drds=drdphi, drdt=drdtheta
            dp_r2       = sum(rvec.^2);
            dp_r_drdt   = rvec * drdtheta_vec.';
            dp_drdt2    = sum(drdtheta_vec.^2);
            dp_r_drds   = rvec * drdphi_vec.';
            dp_dt_ds    = drdtheta_vec * drdphi_vec.';
            dp_drds2    = sum(drdphi_vec.^2);
            
            % Theta direction: drds=drdtheta, drdt=drdphi (swap roles)
            dt_r2       = dp_r2;
            dt_r_drdt   = dp_r_drds;
            dt_drdt2    = dp_drds2;
            dt_r_drds   = dp_r_drdt;
            dt_dt_ds    = dp_dt_ds;
            dt_drds2    = dp_drdt2;
            
            % --- Scaling constants ---
            norm_drdtheta = sqrt(dp_drdt2);
            norm_drdphi = sqrt(dp_drds2);
            C_TZ = grid.nph * norm_drdtheta / norm_drdphi;
            C_GL = max(grid.nth * norm_drdphi / norm_drdtheta, grid.nth);
            
            pm1 = p - 1;
            
            % ===============================================
            % Phi direction estimate (trapezoidal rule)
            % ===============================================
            ds_phi = xlag / C_TZ;
            dp0_ref = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dp_r2, dp_r_drdt, dp_drdt2, dp_r_drds, dp_dt_ds, dp_drds2, 0);
            dp0_pos = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dp_r2, dp_r_drdt, dp_drdt2, dp_r_drds, dp_dt_ds, dp_drds2, ds_phi);
            dp0_neg = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dp_r2, dp_r_drdt, dp_drdt2, dp_r_drds, dp_dt_ds, dp_drds2, -ds_phi);
            
            phi0_candidates = [phi0; conj(phi0)];
            estTZ_both = zeros(2, 3);
            
            for ii = 1:2
                phi0_test = phi0_candidates(ii);
                
                % Inline kernel at phi root — precomputed geometry at theta_star
                srfvec_phi = [at_star * cos(phi0_test), at_star * sin(phi0_test), ct_star];
                nvec_phi = [-dcdt_star * at_star * cos(phi0_test), ...
                            -dcdt_star * at_star * sin(phi0_test), ...
                             dadt_star * at_star];
                r_phi = srfvec_phi - target;
                rn = sum(r_phi .* nvec_phi);
                
                % With unit density e_c: rq=r(c), fp0t = -6*max(|r|)*|r(c)|*|rn|
                % So est_const(c) = base_const * |r(c)|, vectorized over all 3 components
                base_const = 6 * abs(2 / gamma(p)) * max(abs(r_phi)) * abs(rn) * abs(Gphi^(-p));
                
                % Integration (identical for all components)
                dtpos = phi0_test + dp0_pos - dp0_ref;
                dtneg = phi0_test + dp0_neg - dp0_ref;
                
                int_pos = quadest.errorest.UniformEstimateBuilder.trapzErrFunc(dtpos, grid.nph, pm1);
                int_neg = quadest.errorest.UniformEstimateBuilder.trapzErrFunc(dtneg, grid.nph, pm1);
                
                int_sum = sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_TZ);
                
                % All 3 components at once
                estTZ_both(ii, :) = base_const * int_sum * abs(r_phi);
            end
            
            estTZ = min(estTZ_both, [], 1);
            estTZ(isnan(estTZ)) = 0;
            
            % ===============================================
            % Theta direction estimate (Gauss-Legendre)
            % ===============================================
            ds_theta = xlag / C_GL;
            dt0_ref = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dt_r2, dt_r_drdt, dt_drdt2, dt_r_drds, dt_dt_ds, dt_drds2, 0);
            dt0_pos = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dt_r2, dt_r_drdt, dt_drdt2, dt_r_drds, dt_dt_ds, dt_drds2, ds_theta);
            dt0_neg = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dt_r2, dt_r_drdt, dt_drdt2, dt_r_drds, dt_dt_ds, dt_drds2, -ds_theta);
            
            t0_candidates = [t0; conj(t0)];
            estGL_both = zeros(2, 3);
            
            for ii = 1:2
                t0_test = t0_candidates(ii);
                
                % Inline kernel at theta root
                theta0_test = dfac * (t0_test + 1);
                at_t0 = geometry.at(theta0_test);
                ct_t0 = geometry.ct(theta0_test);
                dadt_t0 = geometry.dadt(theta0_test);
                dcdt_t0 = geometry.dcdt(theta0_test);
                
                srfvec_th = [at_t0 * cos(phi_star), at_t0 * sin(phi_star), ct_t0];
                nvec_th = dfac * [-dcdt_t0 * at_t0 * cos(phi_star), ...
                                  -dcdt_t0 * at_t0 * sin(phi_star), ...
                                   dadt_t0 * at_t0];
                r_th = srfvec_th - target;
                rn_th = sum(r_th .* nvec_th);
                
                % Vectorized kernel: base_const * |r_th| gives [1x3]
                base_const = 6 * abs(2 / gamma(p)) * max(abs(r_th)) * abs(rn_th) * abs(Gtheta^(-p));
                
                % Integration (identical for all components)
                dtpos = t0_test - (dt0_pos - dt0_ref);
                dtneg = t0_test - (dt0_neg - dt0_ref);
                
                int_pos = quadest.errorest.UniformEstimateBuilder.glErrFunc(dtpos, grid.nth, pm1);
                int_neg = quadest.errorest.UniformEstimateBuilder.glErrFunc(dtneg, grid.nth, pm1);
                
                int_sum = sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_GL);
                
                % All 3 components at once
                estGL_both(ii, :) = base_const * int_sum * abs(r_th);
            end
            
            estGL = min(estGL_both, [], 1);
            
            % Combined estimate (sum of both directions)
            est = estTZ + estGL;
        end
        
        function est = computeSingleDirectEstimate(geometry, grid, target, itheta, iphi, p, xlag, wlag, dens_3d)
            % COMPUTESINGLEDIRECTESTIMATE Compute error estimate for single target with density
            %
            % Optimized: all closures eliminated, geometry evaluations minimized,
            % kernel evaluation and integration folded inline.
            
            % --- Grid values and mapping constants ---
            theta_star = grid.theta(itheta);
            phi_star = grid.phi(iphi);
            dfac = pi/2;
            
            % --- Phi root (analytic) ---
            [phi0, Gphi] = quadest.errorest.RootFinder.computePhiRoot(geometry, target, theta_star);
            
            % --- Theta root ---
            if target(3) ~= 0 && isa(geometry, 'quadest.geometry.Spheroid')
                % Analytic formula for spheroid
                theta0 = geometry.findThetaRoot(target, phi_star);
                t0 = (2/pi) * theta0 - 1;
            else
                t_star = (2/pi) * theta_star - 1;
                [t0, ~] = quadest.errorest.RootFinder.newtonSolve(...
                    geometry, target, phi_star, t_star, dfac);
                theta0 = dfac * (t0 + 1);
            end
            
            % --- Geometric factor Gtheta (one evaluate+drdtheta call) ---
            gamma_final = geometry.evaluate(theta0, phi_star);
            dgamma_final = dfac * geometry.drdtheta(theta0, phi_star);
            r_theta = gamma_final - target;
            Gtheta = 2 * dot(r_theta, dgamma_final);
            
            % --- Geometry at closest point (shared by kernel + semiRoot) ---
            at_star = geometry.at(theta_star);
            ct_star = geometry.ct(theta_star);
            dadt_star = geometry.dadt(theta_star);
            dcdt_star = geometry.dcdt(theta_star);
            
            rvec = [at_star * cos(phi_star), at_star * sin(phi_star), ct_star] - target;
            drdphi_vec = [-at_star * sin(phi_star), at_star * cos(phi_star), 0];
            drdtheta_vec = [dadt_star * cos(phi_star), dadt_star * sin(phi_star), dcdt_star];
            
            % --- SemiRoot coefficients (avoids closures) ---
            % Phi direction: drds=drdphi, drdt=drdtheta
            dp_r2       = sum(rvec.^2);
            dp_r_drdt   = rvec * drdtheta_vec.';
            dp_drdt2    = sum(drdtheta_vec.^2);
            dp_r_drds   = rvec * drdphi_vec.';
            dp_dt_ds    = drdtheta_vec * drdphi_vec.';
            dp_drds2    = sum(drdphi_vec.^2);
            
            % Theta direction: drds=drdtheta, drdt=drdphi (swap roles)
            dt_r2       = dp_r2;
            dt_r_drdt   = dp_r_drds;
            dt_drdt2    = dp_drds2;
            dt_r_drds   = dp_r_drdt;
            dt_dt_ds    = dp_dt_ds;
            dt_drds2    = dp_drdt2;
            
            % --- Scaling constants ---
            norm_drdtheta = sqrt(dp_drdt2);
            norm_drdphi = sqrt(dp_drds2);
            C_TZ = grid.nph * norm_drdtheta / norm_drdphi;
            C_GL = max(grid.nth * norm_drdphi / norm_drdtheta, grid.nth);
            
            % --- Density interpolation ---
            q_phi = quadest.errorest.DensityModifier.interpDensityPhi(...
                grid, dens_3d, phi0, itheta);
            q_theta = quadest.errorest.DensityModifier.interpDensityTheta(...
                grid, dens_3d, theta0, iphi);
            
            % ===============================================
            % Phi direction estimate (trapezoidal rule)
            % ===============================================
            ds_phi = xlag / C_TZ;
            dp0_ref = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dp_r2, dp_r_drdt, dp_drdt2, dp_r_drds, dp_dt_ds, dp_drds2, 0);
            dp0_pos = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dp_r2, dp_r_drdt, dp_drdt2, dp_r_drds, dp_dt_ds, dp_drds2, ds_phi);
            dp0_neg = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dp_r2, dp_r_drdt, dp_drdt2, dp_r_drds, dp_dt_ds, dp_drds2, -ds_phi);
            
            phi0_candidates = [phi0; conj(phi0)];
            estTZ_both = zeros(2, 1);
            pm1 = p - 1;
            
            for ii = 1:2
                phi0_test = phi0_candidates(ii);
                
                % Inline kernel at phi root
                srfvec_phi = [at_star * cos(phi0_test), at_star * sin(phi0_test), ct_star];
                nvec_phi = [-dcdt_star * at_star * cos(phi0_test), ...
                            -dcdt_star * at_star * sin(phi0_test), ...
                             dadt_star * at_star];
                r_phi = srfvec_phi - target;
                rq = sum(r_phi .* q_phi);
                rn = sum(r_phi .* nvec_phi);
                fp0t = -6 * max(abs(r_phi * (rq * rn)));
                
                % Integration
                est_const = abs(2 / gamma(p) * fp0t * Gphi^(-p));
                
                dtpos = phi0_test + dp0_pos - dp0_ref;
                dtneg = phi0_test + dp0_neg - dp0_ref;
                
                int_pos = quadest.errorest.UniformEstimateBuilder.trapzErrFunc(dtpos, grid.nph, pm1);
                int_neg = quadest.errorest.UniformEstimateBuilder.trapzErrFunc(dtneg, grid.nph, pm1);
                
                estTZ_both(ii) = est_const * sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_TZ);
            end
            
            estTZ = min(estTZ_both);
            if isnan(estTZ); estTZ = 0; end
            
            % ===============================================
            % Theta direction estimate (Gauss-Legendre)
            % ===============================================
            ds_theta = xlag / C_GL;
            dt0_ref = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dt_r2, dt_r_drdt, dt_drdt2, dt_r_drds, dt_dt_ds, dt_drds2, 0);
            dt0_pos = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dt_r2, dt_r_drdt, dt_drdt2, dt_r_drds, dt_dt_ds, dt_drds2, ds_theta);
            dt0_neg = quadest.errorest.UniformEstimateBuilder.semiRootEval(...
                dt_r2, dt_r_drdt, dt_drdt2, dt_r_drds, dt_dt_ds, dt_drds2, -ds_theta);
            
            t0_candidates = [t0; conj(t0)];
            estGL_both = zeros(2, 1);
            
            for ii = 1:2
                t0_test = t0_candidates(ii);
                
                % Inline kernel at theta root
                theta0_test = dfac * (t0_test + 1);
                at_t0 = geometry.at(theta0_test);
                ct_t0 = geometry.ct(theta0_test);
                dadt_t0 = geometry.dadt(theta0_test);
                dcdt_t0 = geometry.dcdt(theta0_test);
                
                srfvec_th = [at_t0 * cos(phi_star), at_t0 * sin(phi_star), ct_t0];
                nvec_th = dfac * [-dcdt_t0 * at_t0 * cos(phi_star), ...
                                  -dcdt_t0 * at_t0 * sin(phi_star), ...
                                   dadt_t0 * at_t0];
                r_th = srfvec_th - target;
                rq_th = sum(r_th .* q_theta);
                rn_th = sum(r_th .* nvec_th);
                ft0p = -6 * max(abs(r_th * (rq_th * rn_th)));
                
                % Integration
                est_const = abs(2 / gamma(p) * ft0p * Gtheta^(-p));
                
                dtpos = t0_test - (dt0_pos - dt0_ref);
                dtneg = t0_test - (dt0_neg - dt0_ref);
                
                int_pos = quadest.errorest.UniformEstimateBuilder.glErrFunc(dtpos, grid.nth, pm1);
                int_neg = quadest.errorest.UniformEstimateBuilder.glErrFunc(dtneg, grid.nth, pm1);
                
                estGL_both(ii) = est_const * sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_GL);
            end
            
            estGL = min(estGL_both);
            
            % Combined estimate
            est = estTZ + estGL;
        end
        
        function val = semiRootEval(r2, r_drdt, drdt2, r_drds, drdt_drds, drds2, dt)
            % SEMIROOTEVAL Evaluate semi-analytical root from precomputed coefficients
            %
            % Inline replacement for RootFinder.semiRoot closures. Takes the 6
            % scalar dot-product coefficients and perturbation dt (scalar or vector).
            aa = r2 + 2*r_drdt*dt + drdt2*dt.^2;
            bb = 2*r_drds + 2*drdt_drds*dt;
            val = -bb./(2*drds2) + 1i*sqrt(aa./drds2 - (bb./(2*drds2)).^2);
        end
        
        function knq = trapzErrFunc(z, n, q)
            % TRAPZERRFUNC q-th derivative of trapezoidal error function
            b = abs(imag(z));
            knq = 2 * pi * n^q * exp(-n * b);
        end
        
        function knq = glErrFunc(z, n, q)
            % GLERRFUNC q-th derivative of Gauss-Legendre error function
            
            % Transform to first quadrant (real part)
            z = abs(real(z)) + 1i * imag(z);
            
            % Joukowski map: rho = z + sqrt(z^2 - 1)
            % For z on/near the imaginary axis with imag(z) < 0, the
            % principal branch can give |rho| < 1, leading to a spurious
            % blow-up.  The two branches rho and 1/rho are reciprocals;
            % always choose the one with |rho| >= 1.
            sq = sqrt(z.^2 - 1);
            rho = z + sq;
            flip = abs(rho) < 1;
            rho(flip) = z(flip) - sq(flip);
            
            knq = 2 * pi ./ rho.^(2*n + 1);
            if q ~= 0
                knq = knq .* (-(2*n + 1) ./ sq).^q;
            end
            
            knq = abs(knq);
        end
    end
end
