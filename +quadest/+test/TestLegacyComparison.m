classdef TestLegacyComparison < quadest.test.TestBase
%TESTLEGACYCOMPARISON Validation tests comparing new vs legacy implementation
%
%   These tests verify numerical equivalence between the new modular
%   quadest package and the original research code in legacy/.
%
%   IMPORTANT: Run from project root with legacy/ on path (init.m does this)
%
%   Example:
%       init;  % adds legacy/ to path
%       results = quadest.test.TestLegacyComparison().run();

    properties (Access = private)
        legacyAvailable = false
    end

    methods
        function setUp(obj)
            %SETUP Check if legacy code is available
            obj.legacyAvailable = exist('init_axsymbody', 'file') == 2;
            if ~obj.legacyAvailable
                warning('Legacy code not on path. Run: addpath(genpath(''legacy''))');
            end
        end

        %% ---- Gauss-Legendre quadrature ----

        function test_gauss_legendre_matches_legacy(obj)
            %TEST_GAUSS_LEGENDRE_MATCHES_LEGACY GL nodes/weights match legacy
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            n = 40;
            [x_new, w_new] = quadest.util.GaussLegendre(n, 0, pi);

            % lgwt is a local subfunction inside init_axsymbody and cannot
            % be called directly. Extract GL nodes from a legacy body.
            geometry_leg = struct('shape', 'spheroid', 'a', 0.05, 'c', 0.1);
            body_leg = init_axsymbody(geometry_leg, 2, n);
            x_leg = body_leg.theta;

            obj.assertAlmostEqual(x_new, x_leg, 1e-12, ...
                'GL nodes should match legacy');

            % Additionally verify via known integrals (40-pt GL is exact
            % for polynomials up to degree 79, so these should be exact)
            err_sin = abs(sum(w_new .* sin(x_new)) - 2);
            err_lin = abs(sum(w_new .* x_new) - pi^2/2);
            obj.assertTrue(max(err_sin, err_lin) < 1e-12, ...
                'GL weights fail known-integral check');
        end

        %% ---- Geometry parameterization ----

        function test_spheroid_geometry_matches(obj)
            %TEST_SPHEROID_GEOMETRY_MATCHES Spheroid parameterization
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            geom_new = quadest.geometry.Spheroid('a', a, 'c', c);

            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body_leg = init_axsymbody(geometry_leg, 60, 40);

            theta_test = [0.1, pi/4, pi/2, 3*pi/4, pi-0.1]';

            for i = 1:length(theta_test)
                th = theta_test(i);
                obj.assertAlmostEqual(geom_new.at(th), body_leg.at(th), 1e-14, ...
                    sprintf('a(theta) mismatch at theta=%.4f', th));
                obj.assertAlmostEqual(geom_new.ct(th), body_leg.ct(th), 1e-14, ...
                    sprintf('c(theta) mismatch at theta=%.4f', th));
                obj.assertAlmostEqual(geom_new.dadt(th), body_leg.dadt(th), 1e-14, ...
                    sprintf('da/dtheta mismatch at theta=%.4f', th));
                obj.assertAlmostEqual(geom_new.dcdt(th), body_leg.dcdt(th), 1e-14, ...
                    sprintf('dc/dtheta mismatch at theta=%.4f', th));
            end
        end

        %% ---- Grid discretization ----

        function test_grid_points_match(obj)
            %TEST_GRID_POINTS_MATCH Grid x/n/w at production resolution
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            nth = 40; nph = 60;

            geom = quadest.geometry.Spheroid('a', a, 'c', c);
            grid_new = quadest.grid.AxsymGrid(geom, 'nth', nth, 'nph', nph);

            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body_leg = init_axsymbody(geometry_leg, nph, nth);

            obj.assertAlmostEqual(grid_new.x, body_leg.x, 1e-12, ...
                'Grid points should match');
            obj.assertAlmostEqual(grid_new.n, body_leg.n, 1e-12, ...
                'Grid normals should match');
            obj.assertAlmostEqual(grid_new.w, body_leg.w, 1e-12, ...
                'Grid weights should match');
        end

        %% ---- Stresslet kernel ----

        function test_stresslet_kernel_matches(obj)
            %TEST_STRESSLET_KERNEL_MATCHES Direct stresslet at 3 targets
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            nth = 40; nph = 60;

            geom = quadest.geometry.Spheroid('a', a, 'c', c);
            grid_new = quadest.grid.AxsymGrid(geom, 'nth', nth, 'nph', nph);

            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body_leg = init_axsymbody(geometry_leg, nph, nth);

            q = randn(nth*nph, 3);
            kernel = quadest.kernel.StokesStresslet();
            targets = [0.2, 0.1, 0.15; 0.3, 0.2, 0.1; 0.1, 0.1, 0.3];

            n_weighted = bsxfun(@times, body_leg.n, body_leg.w);
            for t = 1:size(targets, 1)
                tgt = targets(t, :);
                u_new = kernel.evaluate(tgt, grid_new.x, grid_new.n, q, grid_new.w);
                u_leg = stresslet(tgt, body_leg.x, n_weighted, q);
                obj.assertAlmostEqual(u_new, u_leg, 1e-12, ...
                    sprintf('Stresslet mismatch at target %d', t));
            end
        end

        %% ---- evaluateOnGrid convenience method ----

        function test_evaluate_on_grid(obj)
            %TEST_EVALUATE_ON_GRID Kernel.evaluateOnGrid vs legacy stresslet
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            nth = 10; nph = 16;

            geom = quadest.geometry.Spheroid('a', a, 'c', c);
            grid_new = quadest.grid.AxsymGrid(geom, 'nth', nth, 'nph', nph);

            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body_leg = init_axsymbody(geometry_leg, nph, nth);

            q = randn(nth*nph, 3);
            kernel = quadest.kernel.StokesStresslet();
            target = [0.2, 0.1, 0.15];

            % evaluateOnGrid internally combines n.*w
            u_new = kernel.evaluateOnGrid(target, grid_new, q);

            % Legacy: stresslet expects pre-weighted normals
            n_weighted = bsxfun(@times, body_leg.n, body_leg.w);
            u_leg = stresslet(target, body_leg.x, n_weighted, q);

            obj.assertAlmostEqual(u_new, u_leg, 1e-12, ...
                'evaluateOnGrid should match legacy stresslet');
        end

        %% ---- Upsampled grid geometry ----

        function test_upsampling_matches(obj)
            %TEST_UPSAMPLING_MATCHES Upsampled x/n/w match legacy
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            nth = 10; nph = 16;
            upfac = 2;

            geom = quadest.geometry.Spheroid('a', a, 'c', c);
            grid_new = quadest.grid.AxsymGrid(geom, 'nth', nth, 'nph', nph);

            % New: upsample via AxsymGrid.upsample → creates a new grid
            fineGrid = grid_new.upsample(upfac);

            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body_leg = init_axsymbody(geometry_leg, nph, nth);
            body_up_leg = upsample(body_leg, upfac);

            obj.assertAlmostEqual(fineGrid.x, body_up_leg.x, 1e-11, ...
                'Upsampled points should match');
            obj.assertAlmostEqual(fineGrid.n, body_up_leg.n, 1e-11, ...
                'Upsampled normals should match');
            obj.assertAlmostEqual(fineGrid.w, body_up_leg.w, 1e-11, ...
                'Upsampled weights should match');
        end

        %% ---- Density interpolation ----

        function test_density_upsampling_matches(obj)
            %TEST_DENSITY_UPSAMPLING_MATCHES Upsampler.apply vs legacy interpolate
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            nth = 10; nph = 16;
            upfac = 2;

            geom = quadest.geometry.Spheroid('a', a, 'c', c);
            grid_new = quadest.grid.AxsymGrid(geom, 'nth', nth, 'nph', nph);
            fineGrid = grid_new.upsample(upfac);
            upsampler = quadest.grid.Upsampler(grid_new, fineGrid);

            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body_leg = init_axsymbody(geometry_leg, nph, nth);

            q = grid_new.x;  % Smooth density (position)

            q_up_new = upsampler.apply(q);

            body_up_leg = upsample(body_leg, upfac);
            q_up_leg = body_up_leg.interpolate(q);

            obj.assertAlmostEqual(q_up_new, q_up_leg, 1e-10, ...
                'Upsampled density should match');
        end

        %% ---- Stresslet with upsampling ----

        function test_stresslet_with_upsampling(obj)
            %TEST_STRESSLET_WITH_UPSAMPLING Manual upsample+evaluate vs legacy
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            nth = 40; nph = 60;
            upfac = 2;

            geom = quadest.geometry.Spheroid('a', a, 'c', c);
            grid_new = quadest.grid.AxsymGrid(geom, 'nth', nth, 'nph', nph);
            fineGrid = grid_new.upsample(upfac);
            upsampler = quadest.grid.Upsampler(grid_new, fineGrid);

            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body_leg = init_axsymbody(geometry_leg, nph, nth);

            kernel = quadest.kernel.StokesStresslet();
            q = ones(nth*nph, 3);
            tgt = [0.06, 0.01, 0.02];

            % New: manual upsample + evaluate
            q_up = upsampler.apply(q);
            u_new = kernel.evaluate(tgt, fineGrid.x, fineGrid.n, q_up, fineGrid.w);

            % Legacy: stresslet_with_upsamp
            u_leg = stresslet_with_upsamp(tgt, q, body_leg, upfac);

            obj.assertAlmostEqual(u_new, u_leg, 1e-10, ...
                'Stresslet with upsampling should match legacy');
        end

        %% ---- evaluateUpsampled convenience method ----

        function test_evaluate_upsampled(obj)
            %TEST_EVALUATE_UPSAMPLED Kernel.evaluateUpsampled vs legacy
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            nth = 10; nph = 16;
            upfac = 2;

            geom = quadest.geometry.Spheroid('a', a, 'c', c);
            grid_new = quadest.grid.AxsymGrid(geom, 'nth', nth, 'nph', nph);

            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body_leg = init_axsymbody(geometry_leg, nph, nth);

            kernel = quadest.kernel.StokesStresslet();
            q = ones(nth*nph, 3);
            tgt = [0.06, 0.01, 0.02];

            % New: one-call upsampled evaluation
            u_new = kernel.evaluateUpsampled(tgt, grid_new, q, upfac);

            % Legacy
            u_leg = stresslet_with_upsamp(tgt, q, body_leg, upfac);

            obj.assertAlmostEqual(u_new, u_leg, 1e-10, ...
                'evaluateUpsampled should match legacy stresslet_with_upsamp');
        end

        %% ---- Theta root finding ----

        function test_theta_root_spheroid_matches(obj)
            %TEST_THETA_ROOT_SPHEROID_MATCHES Quartic root vs legacy
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            geom = quadest.geometry.Spheroid('a', a, 'c', c);

            % Test several exterior target points (columns: rxy, z_abs)
            test_cases = [
                0.06, 0.08;
                0.04, 0.05;
                0.055, 0.01;
                0.03, 0.07;
            ];

            for i = 1:size(test_cases, 1)
                rxy   = test_cases(i, 1);
                z_abs = test_cases(i, 2);

                % Both APIs: 3D target point + phi value
                target_3d = [rxy, 0, z_abs];
                phi_val   = 0;

                theta_new = geom.findThetaRoot(target_3d, phi_val);
                theta_leg = find_theta_root_spheroid(target_3d, a, c, phi_val);

                if ~isnan(theta_new) && ~isnan(theta_leg)
                    obj.assertAlmostEqual(theta_new, theta_leg, 1e-10, ...
                        sprintf('Theta root mismatch at case %d', i));
                end
            end
        end

        %% ---- Exterior point classification ----

        function test_is_exterior_matches(obj)
            %TEST_IS_EXTERIOR_MATCHES isExterior vs find_exterior_points
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            geom = quadest.geometry.Spheroid('a', a, 'c', c);

            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body_leg = init_axsymbody(geometry_leg, 60, 40);

            % Mix of interior, exterior, and near-surface points
            points = [
                0,    0,    0;         % interior (origin)
                0.03, 0,    0;         % interior
                0.2,  0.1,  0.05;      % exterior (far)
                0.06, 0.01, 0.02;      % exterior (near equator)
                0.01, 0,    0.12;      % exterior (near pole)
                0,    0,    0.05;      % interior (on axis)
                0,    0,    0.15;      % exterior (on axis)
                0.049, 0,   0;         % interior (just inside equator)
                0.051, 0,   0;         % exterior (just outside equator)
            ];

            mask_new = geom.isExterior(points);
            [~, mask_leg] = find_exterior_points(body_leg, points);

            obj.assertEqual(mask_new, mask_leg, ...
                'isExterior should match legacy find_exterior_points');
        end

        function test_uniform_single_estimates_match_legacy_reference(obj)
            %TEST_UNIFORM_SINGLE_ESTIMATES_MATCH_LEGACY_REFERENCE Direct point estimates
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            nth = 20; nph = 30;

            geom = quadest.geometry.Spheroid('a', a, 'c', c);
            grid = quadest.grid.AxsymGrid(geom, 'nth', nth, 'nph', nph);
            kernel = quadest.kernel.StokesStresslet();

            targets = [
                0.0183673469388, 0, 0.1469387755100;
                0.0704081632653, 0, 0.1224489795920;
                0.0979591836735, 0, 0.0183673469388;
            ];
            expected = [
                1.30870479819e-06, 3.38102441411e-39, 1.30870479819e-06;
                1.49072312607e-05, 2.40675263166e-15, 1.49072312617e-05;
                4.84002446179e-03, 2.13371363950e-06, 4.83829323757e-03;
            ];

            estimates = quadest.errorest.UniformEstimateBuilder.computeEstimates( ...
                geom, grid, kernel, targets, kernel.singularityOrder());

            absTol = 1e-13;
            relTol = 1e-10;
            err = abs(estimates - expected);
            allowed = absTol + relTol * abs(expected);
            obj.assertTrue(all(err(:) <= allowed(:)), ...
                sprintf('Uniform single-point estimates differ from legacy reference (max scaled err=%.3g)', ...
                    max(err(:) ./ allowed(:))));
        end

        %% ---- Full error estimation workflow ----

        function test_full_error_estimate_workflow(obj)
            %TEST_FULL_ERROR_ESTIMATE_WORKFLOW Uniform estimates + end-to-end
            if ~obj.legacyAvailable; obj.skip('Legacy code not available'); end

            a = 0.05; c = 0.1;
            nth = 20; nph = 30;
            upsamp_fac = 1:2;

            % New implementation - precompute
            geom = quadest.geometry.Spheroid('a', a, 'c', c);
            kernel = quadest.kernel.StokesStresslet();
            estimator = quadest.errorest.ErrorEstimator(geom, kernel, ...
                'nth', nth, 'nph', nph, 'upsampFactors', upsamp_fac);

            % Legacy implementation - precompute
            geometry_leg = struct('shape', 'spheroid', 'a', a, 'c', c);
            body = init_axsymbody(geometry_leg, nph, nth);
            errestinterp = uniform_error_estimates(body, upsamp_fac, true);

            % --- Part A: Compare uniform estimate interpolants ---
            % Query both interpolant sets at a grid of sample points
            rxy_test = linspace(0.02, 0.2, 8).';
            z_test   = linspace(0.02, 0.2, 8).';
            [Rxy, Zabs] = ndgrid(rxy_test, z_test);
            rxy_flat = Rxy(:);
            z_flat   = Zabs(:);

            % Upsampling factor indices differ between new (indexed by
            % actual factor value) and legacy (indexed sequentially).
            % Legacy upsampfacs = unique([1, upsamp_fac]) → [1 2 3 4].
            upfacvals = unique([1, upsamp_fac]);
            for k = 1:length(upfacvals)
                upfac = upfacvals(k);
                for comp = 1:3
                    % Legacy: errest_INT{k, comp} is log10(estimate)
                    val_leg = zeros(size(rxy_flat));
                    for ii = 1:length(rxy_flat)
                        val_leg(ii) = errestinterp.errest_INT{k, comp}(rxy_flat(ii), z_flat(ii));
                    end

                    % New: interpolants{upfac, comp} is log10(estimate)
                    val_new = zeros(size(rxy_flat));
                    for ii = 1:length(rxy_flat)
                        val_new(ii) = estimator.interpolants{upfac, comp}(rxy_flat(ii), z_flat(ii));
                    end

                    rel_diff = max(abs(10.^val_new-10.^val_leg)./abs(10.^val_leg));
                    obj.assertTrue(rel_diff < 1e-10, ...
                        sprintf('Interpolant mismatch for upfac=%d, comp=%d (rel_diff=%.2g)', ...
                        upfac, comp, rel_diff));
                end
            end

            % --- Part B: End-to-end evaluate with stresslet identity ---
            targets = [
                0.2,  0.1,  0.05;   % Far field
                0.06, 0.01, 0.02;   % Near equator
                0.02, 0.01, 0.12;   % Near pole
            ];
            q = ones(nth*nph, 3);

            for t = 1:size(targets, 1)
                tgt = targets(t, :);
                est_new = estimator.evaluate(tgt, q);

                % Replicate legacy workflow for this target
                xycoord = sqrt(tgt(1)^2 + tgt(2)^2);
                zcoord  = abs(tgt(3));
                esttmp  = cellfun(@(c) 10.^(c(xycoord, zcoord)), ...
                    {errestinterp.errest_INT{1, :}}, 'UniformOutput', false);
                estval_leg = cell2mat(esttmp);

                % Density modifier via legacy Qintp
                tmpR = knnsearch(body.x, tgt);
                dens_3d = reshape(q, body.nth, body.nph, 3);
                if errestinterp.flag
                    th_intp = errestinterp.thRe_INT{1}(xycoord, zcoord) ...
                              + 1i * errestinterp.thIm_INT{1}(xycoord, zcoord);
                    [qp, qt] = Qintp(body, tgt, body.a, body.c, dens_3d, tmpR, th_intp);
                else
                    [qp, qt] = Qintp(body, tgt, body.a, body.c, dens_3d, tmpR);
                end
                Mqpqt = max(qp, qt);
                est_leg = max(estval_leg .* Mqpqt);

                % Estimates should agree within reasonable tolerance
                % (interpolation grids and tabulation may differ slightly)
                if est_leg > 0 && est_new > 0
                    rel_err = abs(est_new - est_leg)/abs(est_leg);
                    obj.assertTrue(rel_err < 1e-10, ...
                        sprintf('Estimate mismatch at target %d: new=%.4g, leg=%.4g (rel diff=%.2g)', ...
                            t, est_new, est_leg, rel_err));
                else
                    obj.assertTrue(est_new >= 0, ...
                        sprintf('Estimate should be non-negative at target %d', t));
                end
            end
        end
    end
end
