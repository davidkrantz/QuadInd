classdef UniformIndicatorBuilder
    % UNIFORMINDICATORBUILDER Precompute uniform indicators on a tabulation grid
    %
    % Builds gridded interpolants for error indicators assuming unit density.
    % The actual error indicator is obtained by multiplying by the density
    % magnitude at complex roots (see DensityModifier).
    %
    % This is the most computationally expensive part of the precomputation
    % and should only be done once per geometry/discretization.
    %
    % See also: quadind.IndicatorEvaluator, quadind.indicator.DensityModifier
    
    properties (Constant, Access = private)
        % Singularity order for stresslet kernel
        STRESSLET_P = 5/2
    end
    
    methods (Static)
        function [interpolants, rootInterp, tabGrid] = build(geometry, grid, kernel, upsampFactors, interpolateRoots, config)
            % BUILD Create gridded interpolants for error indicators
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
            %                  for log10(error indicator)
            %   rootInterp   - Struct with theta root interpolants (if interpolateRoots)
            %   tabGrid      - Struct with tabulation grid info (for debugging)
            
            arguments
                geometry (1,1) quadind.geometry.AxsymGeometry
                grid (1,1) quadind.grid.AxsymGrid
                kernel (1,1) quadind.kernel.Kernel
                upsampFactors (1,:) {mustBePositive, mustBeInteger}
                interpolateRoots (1,1) logical = true
                config (1,1) quadind.util.Config = quadind.util.Config()
            end
            
            quadind.util.Diagnostics.info('Building uniform indicator interpolants...');
            
            % Ensure upsampling factors include 1 (direct quadrature)
            upsampFactors = unique([1, upsampFactors]);
            nfac = length(upsampFactors);
            ncomp = kernel.numComponents();
            p = kernel.singularityOrder();

            if ~isa(kernel, 'quadind.kernel.StokesStresslet') || ncomp ~= 3
                error('quadind:UniformIndicatorBuilder:unsupportedKernel', ...
                    'UniformIndicatorBuilder currently supports only StokesStresslet.');
            end
            
            % Build tabulation grid
            [tabGrid, evalPoints, evalMask] = quadind.indicator.UniformIndicatorBuilder.buildTabGrid(geometry, config);
            
            % Initialize storage
            interpolants = cell(max(upsampFactors), ncomp);
            if interpolateRoots
                rootInterp.Re = cell(max(upsampFactors), 1);
                rootInterp.Im = cell(max(upsampFactors), 1);
            else
                rootInterp = [];
            end
            
            % Find axis points that need special handling
            axisPoints = quadind.indicator.UniformIndicatorBuilder.findAxisPoints(evalPoints, evalMask);
            
            % Loop through upsampling factors
            for j = 1:nfac
                upfac = upsampFactors(j);
                quadind.util.Diagnostics.info('  Upsampling factor %d/%d (kappa=%d)', j, nfac, upfac);
                
                % Create upsampled grid
                if upfac == 1
                    refinedGrid = grid;
                else
                    refinedGrid = grid.upsample(upfac);
                end
                
                % Compute error indicators at exterior tabulation points
                [indicators, thetaRoots] = quadind.indicator.UniformIndicatorBuilder.computeIndicators(...
                    geometry, refinedGrid, kernel, evalPoints(evalMask, :), p);
                
                % Build full indicator array (interior points get large values)
                fullIndicators = ones(size(evalPoints, 1), ncomp) * config.interiorIndicator;
                fullIndicators(evalMask, :) = indicators;
                
                % Handle NaN values
                fullIndicators(isnan(fullIndicators)) = config.interiorIndicator;
                
                % Safety: axis point closest to surface
                if ~isempty(axisPoints)
                    fullIndicators(axisPoints, :) = config.interiorIndicator;
                end
                
                % Add lower bound
                fullIndicators = fullIndicators + config.indicatorLowerBound;
                
                % Build interpolants for each component
                for c = 1:ncomp
                    estMat = reshape(fullIndicators(:, c), config.ntab, config.nztab);
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
            
            quadind.util.Diagnostics.info('Done building interpolants.');
        end
        
        function [indicators, thetaRoots] = computeIndicators(geometry, grid, kernel, targets, p)
            % COMPUTEINDICATORS Compute uniform error indicators at target points
            %
            %   [indicators, thetaRoots] = computeIndicators(geometry, grid, kernel, targets, p)
            %
            % Computes per-component error indicators assuming unit density e_c
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
            %   indicators  - [M×ncomp] per-component error indicators
            %   thetaRoots - [M×1] complex theta roots
            
            M = size(targets, 1);
            ncomp = kernel.numComponents();
            
            indicators = zeros(M, ncomp);
            thetaRoots = zeros(M, 1);
            
            % Find closest grid points using parametric approach
            [itheta_all, iphi_all] = quadind.indicator.UniformIndicatorBuilder.findNearestNodes(grid, targets);
            
            % Get Gauss-Laguerre quadrature
            [xlag, wlag] = quadind.util.GaussLaguerre8();
            
            % Loop through target points
            for ii = 1:M
                target = targets(ii, :);
                
                [est, theta0] = quadind.indicator.UniformIndicatorBuilder.computeSingleIndicator(...
                    geometry, grid, target, itheta_all(ii), iphi_all(ii), p, xlag, wlag);

                indicators(ii, :) = est;
                thetaRoots(ii) = theta0;
            end
        end
        
        function indicators = computeDirectIndicators(geometry, grid, kernel, targets, density)
            % COMPUTEDIRECTINDICATORS Compute error indicators with given density
            %
            %   indicators = computeDirectIndicators(geometry, grid, kernel, targets, density)
            %
            % Computes the same per-component uniform indicators used to build
            % the tabulated interpolants, but directly at each target point.
            % The supplied density is then applied through the same root-density
            % modifier used by IndicatorEvaluator.evaluate().
            %
            % Inputs:
            %   geometry - AxsymGeometry object
            %   grid     - AxsymGrid object
            %   kernel   - Kernel object
            %   targets  - [M×3] target points (must be exterior)
            %   density  - [N×ncomp] density on grid (N = nth*nph)
            %
            % Outputs:
            %   indicators - [M×1] componentwise infinity-norm error indicators

            ncomp = kernel.numComponents();
            p = kernel.singularityOrder();

            if size(density, 1) ~= grid.numPoints()
                error('quadind:UniformIndicatorBuilder:invalidDensity', ...
                    'Density must have %d rows (got %d)', grid.numPoints(), size(density, 1));
            end
            if size(density, 2) ~= ncomp
                error('quadind:UniformIndicatorBuilder:invalidDensity', ...
                    'Density must have %d columns (got %d)', ncomp, size(density, 2));
            end

            [uniformIndicators, thetaRoots] = ...
                quadind.indicator.UniformIndicatorBuilder.computeIndicators(...
                    geometry, grid, kernel, targets, p);

            [qphi, qtheta] = quadind.indicator.DensityModifier.computeDensityAtRoots(...
                grid, targets, density, thetaRoots);
            modifier = quadind.indicator.DensityModifier.computeModifier(qphi, qtheta, targets);

            indicators = sum(uniformIndicators .* modifier, 2);
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
        
        function axisPoints = findAxisPoints(evalPoints, evalMask)
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
        
        function [est, theta0] = computeSingleIndicator(geometry, grid, target, itheta, iphi, p, xlag, wlag)
            % COMPUTESINGLEINDICATOR Compute error indicator for single target
            %
            % Optimized: closures eliminated, geometry evaluations minimized,
            % kernel and integration folded inline, component loop vectorized.
            
            % --- Grid values and mapping constants ---
            theta_star = grid.theta(itheta);
            phi_star = grid.phi(iphi);
            dfac = pi/2;
            
            % --- Phi root (analytic) ---
            [phi0, Gphi] = quadind.indicator.RootFinder.computePhiRoot(geometry, target, theta_star);
            
            % --- Theta root ---
            if isa(geometry, 'quadind.geometry.Spheroid')
                theta0 = geometry.findThetaRoot(target, phi_star);
            else
                t_star = (2/pi) * theta_star - 1;
                [t0, converged] = quadind.indicator.RootFinder.newtonSolve(...
                    geometry, target, phi_star, t_star, dfac);
                if ~converged
                    error('quadind:UniformIndicatorBuilder:thetaRootFailure', ...
                        'Theta-root Newton iteration failed for target [%g %g %g].', target);
                end
                theta0 = dfac * (t0 + 1);
            end
            t0 = (2/pi) * theta0 - 1;
            
            % --- Geometric factor Gtheta (one evaluate+drdtheta call) ---
            gamma_final = geometry.evaluate(theta0, phi_star);
            dgamma_final = dfac * geometry.drdtheta(theta0, phi_star);
            r_theta = gamma_final - target;
            Gtheta = 2 * sum(r_theta .* dgamma_final);
            
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
            C_GL = 2 * grid.nth * norm_drdphi / norm_drdtheta;
            
            pm1 = p - 1;
            
            % ===============================================
            % Phi direction indicator (trapezoidal rule)
            % ===============================================
            onAxis = hypot(target(1), target(2)) == 0;
            ds_phi = xlag / C_TZ;
            dp0_ref = quadind.indicator.UniformIndicatorBuilder.semiRootEval(...
                dp_r2, dp_r_drdt, dp_drdt2, dp_r_drds, dp_dt_ds, dp_drds2, 0);
            dp0_pos = quadind.indicator.UniformIndicatorBuilder.semiRootEval(...
                dp_r2, dp_r_drdt, dp_drdt2, dp_r_drds, dp_dt_ds, dp_drds2, ds_phi);
            dp0_neg = quadind.indicator.UniformIndicatorBuilder.semiRootEval(...
                dp_r2, dp_r_drdt, dp_drdt2, dp_r_drds, dp_dt_ds, dp_drds2, -ds_phi);
            
            phi0_candidates = [phi0; conj(phi0)];
            estTZ_both = zeros(2, 3);
            
            for ii = 1:2
                if onAxis
                    continue;
                end
                phi0_test = phi0_candidates(ii);
                
                % Inline kernel at phi root — precomputed geometry at theta_star
                srfvec_phi = [at_star * cos(phi0_test), at_star * sin(phi0_test), ct_star];
                nvec_phi = [-dcdt_star * at_star * cos(phi0_test), ...
                            -dcdt_star * at_star * sin(phi0_test), ...
                             dadt_star * at_star];
                
                % Integration (identical for all components)
                dtpos = phi0_test - dp0_pos + dp0_ref;
                dtneg = phi0_test - dp0_neg + dp0_ref;
                
                int_pos = quadind.indicator.UniformIndicatorBuilder.trapzErrFunc(dtpos, grid.nph, pm1);
                int_neg = quadind.indicator.UniformIndicatorBuilder.trapzErrFunc(dtneg, grid.nph, pm1);
                
                int_sum = sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_TZ);
                
                for c = 1:3
                    density = zeros(1, 3);
                    density(c) = 1.0;
                    fp0t = quadind.indicator.UniformIndicatorBuilder.stressletRootValue( ...
                        srfvec_phi, target, density, nvec_phi);
                    est_const = abs(2 / gamma(p) * fp0t * Gphi^(-p));
                    estTZ_both(ii, c) = est_const * int_sum;
                end
            end
            
            estTZ = min(estTZ_both, [], 1);
            estTZ(isnan(estTZ)) = 0;
            
            % ===============================================
            % Theta direction indicator (Gauss-Legendre)
            % ===============================================
            ds_theta = xlag / C_GL;
            dt0_ref = quadind.indicator.UniformIndicatorBuilder.semiRootEval(...
                dt_r2, dt_r_drdt, dt_drdt2, dt_r_drds, dt_dt_ds, dt_drds2, 0);
            dt0_pos = quadind.indicator.UniformIndicatorBuilder.semiRootEval(...
                dt_r2, dt_r_drdt, dt_drdt2, dt_r_drds, dt_dt_ds, dt_drds2, ds_theta);
            dt0_neg = quadind.indicator.UniformIndicatorBuilder.semiRootEval(...
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
                
                dtpos = t0_test - dt0_pos + dt0_ref;
                dtneg = t0_test - dt0_neg + dt0_ref;
                
                int_pos = quadind.indicator.UniformIndicatorBuilder.glErrFunc(dtpos, grid.nth, pm1);
                int_neg = quadind.indicator.UniformIndicatorBuilder.glErrFunc(dtneg, grid.nth, pm1);
                
                int_sum = sum((int_pos + int_neg) .* exp(xlag) .* wlag / C_GL);
                
                for c = 1:3
                    density = zeros(1, 3);
                    density(c) = 1.0;
                    ft0p = quadind.indicator.UniformIndicatorBuilder.stressletRootValue( ...
                        srfvec_th, target, density, nvec_th);
                    est_const = abs(2 / gamma(p) * ft0p * Gtheta^(-p));
                    estGL_both(ii, c) = est_const * int_sum;
                end
            end
            
            estGL = min(estGL_both, [], 1);
            
            % Combined indicator (sum of both directions)
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

        function val = stressletRootValue(source, target, density, normal)
            % STRESSLETROOTVALUE Stresslet numerator at a complex root.
            %
            % Use the analytic bilinear continuation
            % r_i (r_j q_j) (r_k n_k). Do not conjugate r in the normal
            % contraction; conjugation would break analyticity at the root.

            r = zeros(size(target, 1), size(source, 1), 3);
            for d = 1:3
                xd = source(:, d).';
                yd = target(:, d);
                r(:, :, d) = bsxfun(@minus, xd, yd);
            end

            q = reshape(density, [size(density, 1), 1, 3]);
            n = reshape(normal, [size(normal, 1), 1, 3]);
            rq = sum(bsxfun(@times, r, q), 3);
            rn = sum(bsxfun(@times, r, n), 3);
            f = bsxfun(@times, r, rq .* rn);
            val = -6 * max(abs(f));
        end
        
        function knq = trapzErrFunc(z, n, q)
            % TRAPZERRFUNC q-th derivative of trapezoidal error function
            b = abs(imag(z));
            knq = 2 * pi * n^q * exp(-n * b);
        end
        
        function knq = glErrFunc(z, n, q)
            % GLERRFUNC q-th derivative of Gauss-Legendre error function
            
            % Branch cut handling according to paper
            sq = sqrt(z+1).*sqrt(z-1);
            rho = z + sq;
            knq = 2 * pi ./ rho.^(2*n + 1);

            if q ~= 0
                knq = knq .* (-(2*n + 1) ./ sq).^q;
            end
            
            knq = abs(knq);
        end
    end
end
