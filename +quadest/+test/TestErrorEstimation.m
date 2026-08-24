classdef TestErrorEstimation < quadest.test.TestBase
%TESTERRORESTIMATION Unit tests for +quadest/+errorest module
%
%   Tests for ErrorEstimator, RootFinder, DensityModifier, UniformEstimateBuilder

    properties (Access = private)
        geom
        kernel
    end
    
    methods
        function setUp(obj)
            %SETUP Create test fixtures
            obj.geom = quadest.geometry.Spheroid('a', 0.05, 'c', 0.1);
            obj.kernel = quadest.kernel.StokesStresslet();
        end
        
        function test_densitymodifier_interpolate(obj)
            %TEST_DENSITYMODIFIER_INTERPOLATE Verify density at roots computation
            % Create grid
            grid = quadest.grid.AxsymGrid(obj.geom, 'nth', 10, 'nph', 16);
            
            % Create smooth density
            q = ones(grid.numPoints(), 3);
            
            % Target point
            targets = [0.06, 0, 0.05];
            
            % Compute density at roots
            [qphi, qtheta, ~, ~] = quadest.errorest.DensityModifier.computeDensityAtRoots(...
                grid, targets, q);
            
            obj.assertSize(qphi, [1, 3]);
            obj.assertSize(qtheta, [1, 3]);
        end

        function test_densitymodifier_matches_local_linear_continuation(obj)
            grid = quadest.grid.AxsymGrid(obj.geom, 'nth', 18, 'nph', 32);
            [theta, phi] = grid.meshgrid();
            density = [phi(:), theta(:), 2*phi(:) - 3*theta(:)];
            alpha = 0.7;
            target = [0.09*cos(alpha), 0.09*sin(alpha), 0.04];

            [qphi, qtheta, phi0] = ...
                quadest.errorest.DensityModifier.computeDensityAtRoots( ...
                    grid, target, density);
            [itheta, iphi] = ...
                quadest.errorest.UniformEstimateBuilder.findNearestNodes(grid, target);
            thetaStar = grid.theta(itheta);
            phiStar = grid.phi(iphi);
            theta0 = obj.geom.findThetaRoot(target, phiStar);

            expectedPhi = [phi0, thetaStar, 2*phi0 - 3*thetaStar];
            expectedTheta = [phiStar, theta0, 2*phiStar - 3*theta0];
            obj.assertAlmostEqual(qphi, expectedPhi, 1e-11);
            obj.assertAlmostEqual(qtheta, expectedTheta, 1e-10);
        end

        function test_densitymodifier_smooth_continuation_converges(obj)
            target = [0.065*cos(2*pi-0.03), 0.065*sin(2*pi-0.03), 0.035];
            errors = zeros(2,2);
            resolutions = [16, 24; 32, 48];
            for k = 1:2
                grid = quadest.grid.AxsymGrid(obj.geom, ...
                    'nth', resolutions(k,1), 'nph', resolutions(k,2));
                [theta, phi] = grid.meshgrid();
                density = [cos(2*phi(:)), cos(2*theta(:)), ...
                    cos(2*phi(:)).*cos(2*theta(:))];
                [qphi, qtheta, phi0] = ...
                    quadest.errorest.DensityModifier.computeDensityAtRoots( ...
                        grid, target, density);
                [~, iphi] = ...
                    quadest.errorest.UniformEstimateBuilder.findNearestNodes(grid, target);
                theta0 = obj.geom.findThetaRoot(target, grid.phi(iphi));
                errors(k,1) = abs(qphi(1) - cos(2*phi0));
                errors(k,2) = abs(qtheta(2) - cos(2*theta0));
            end
            obj.assertTrue(all(errors(2,:) < errors(1,:)), ...
                sprintf('Complex density interpolation did not converge: %s', mat2str(errors)));
        end

        function test_complex_geometry_factors_are_bilinear(obj)
            target = [0.08, 0.03, 0.04];
            thetaStar = 1.1;
            [phi0, Gphi] = quadest.errorest.RootFinder.computePhiRoot( ...
                obj.geom, target, thetaStar);
            rphi = obj.geom.evaluate(thetaStar, phi0) - target;
            expectedPhi = 2 * sum(rphi .* obj.geom.drdphi(thetaStar, phi0), 2);
            obj.assertAlmostEqual(Gphi, expectedPhi, 1e-12);

            theta0 = obj.geom.findThetaRoot(target, 0);
            [~, Gtheta] = quadest.errorest.RootFinder.computeThetaRoot( ...
                obj.geom, target, 0, 1); %#ok<ASGLU>
            rtheta = obj.geom.evaluate(theta0, 0) - target;
            expectedTheta = pi * sum(rtheta .* obj.geom.drdtheta(theta0, 0), 2);
            obj.assertAlmostEqual(Gtheta, expectedTheta, 1e-12);
        end

        function test_spheroid_roots_satisfy_analytic_equation(obj)
            geometries = {obj.geom, ...
                quadest.geometry.Spheroid('a', 0.1, 'c', 0.04), ...
                quadest.geometry.Spheroid('a', 0.08, 'c', 0.08)};
            targets = [0.06, 0, 0.04; 0.03, 0.02, 0.09; 0, 0, 0.12];
            phi = mod(atan2(targets(:,2), targets(:,1)), 2*pi);
            for k = 1:numel(geometries)
                roots = geometries{k}.findThetaRoot(targets, phi);
                displacement = geometries{k}.evaluate(roots, phi) - targets;
                residual = abs(sum(displacement.^2, 2)) ./ ...
                    max(sum(abs(displacement).^2, 2), realmin);
                obj.assertTrue(all(isfinite(roots)));
                obj.assertTrue(max(residual) < 1e-9, ...
                    sprintf('Maximum normalized root residual was %.3e.', max(residual)));
            end
        end

        function test_densitymodifier_axis_is_explicit_and_finite(obj)
            grid = quadest.grid.AxsymGrid(obj.geom, 'nth', 12, 'nph', 16);
            density = ones(grid.numPoints(), 3);
            [qphi, qtheta, phi0] = ...
                quadest.errorest.DensityModifier.computeDensityAtRoots( ...
                    grid, [0, 0, 0.2], density);
            obj.assertAlmostEqual(qphi, zeros(1,3), 0);
            obj.assertAlmostEqual(qtheta, ones(1,3), 1e-12);
            obj.assertTrue(isinf(imag(phi0)));

            estimate = quadest.errorest.UniformEstimateBuilder.computeDirectEstimates( ...
                obj.geom, grid, obj.kernel, [0, 0, 0.2], density);
            obj.assertTrue(isfinite(estimate) && estimate >= 0);
        end

        function test_estimator_rejects_wrong_density_component_count(obj)
            cfg = quadest.util.Config('ntab', 12, 'nztab', 12);
            estimator = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1, 'config', cfg);
            density = ones(10*16, 2);
            obj.assertError(@() estimator.evaluate([0.2, 0, 0], density), ...
                'quadest:ErrorEstimator:invalidDensity');
        end

        function test_tabulation_safety_factor_scales_indicator(obj)
            cfg1 = quadest.util.Config('ntab', 20, 'nztab', 20, ...
                'tabulationSafetyFactor', 1);
            cfg2 = quadest.util.Config('ntab', 20, 'nztab', 20, ...
                'tabulationSafetyFactor', 2);
            estimator1 = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1, 'config', cfg1);
            estimator2 = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1, 'config', cfg2);
            target = [0.12, 0.03, 0.04];
            density = ones(10*16, 3);
            estimate1 = estimator1.evaluate(target, density);
            estimate2 = estimator2.evaluate(target, density);
            obj.assertAlmostEqual(estimate2, 2*estimate1, ...
                1e-12*max(1, estimate2));
        end

        function test_densitymodifier_rotates_to_target_frame(obj)
            % Root densities are global Cartesian vectors; modifiers use
            % radial, azimuthal, and axial components at each target.
            qphi = [1+2i, 3-1i, 4; 2, -1i, 5];
            qtheta = [2-1i, -1+3i, 1; -3i, 4, -2];
            targets = [0, 2, 1; -1, 0, -1];

            modifier = quadest.errorest.DensityModifier.computeModifier( ...
                qphi, qtheta, targets);

            alpha = atan2(targets(:,2), targets(:,1));
            c = cos(alpha);
            s = sin(alpha);
            sz = [1; -1];
            qphiExpected = [c.*qphi(:,1) + s.*qphi(:,2), ...
                           -s.*qphi(:,1) + c.*qphi(:,2), sz.*qphi(:,3)];
            qthetaExpected = [c.*qtheta(:,1) + s.*qtheta(:,2), ...
                             -s.*qtheta(:,1) + c.*qtheta(:,2), sz.*qtheta(:,3)];
            expected = max(abs(qphiExpected), abs(qthetaExpected));

            obj.assertAlmostEqual(modifier, expected, 1e-14);
        end

        function test_densitymodifier_rotation_covariance(obj)
            qphi = [1+2i, 3-1i, 4];
            qtheta = [2-1i, -1+3i, 1];
            target = [0.2, -0.1, 0.3];
            beta = pi/2;
            R = [cos(beta), -sin(beta), 0; ...
                 sin(beta),  cos(beta), 0; 0, 0, 1];

            base = quadest.errorest.DensityModifier.computeModifier( ...
                qphi, qtheta, target);
            rotated = quadest.errorest.DensityModifier.computeModifier( ...
                qphi * R.', qtheta * R.', target * R.');

            obj.assertAlmostEqual(rotated, base, 1e-14);
        end

        function test_root_interpolation_rotation_covariance(obj)
            % A grid-aligned active rotation of the target and vector field
            % must commute with both complex-root interpolations.
            grid = quadest.grid.AxsymGrid(obj.geom, 'nth', 10, 'nph', 16);
            beta = pi/2;
            R = [cos(beta), -sin(beta), 0; ...
                 sin(beta),  cos(beta), 0; 0, 0, 1];
            target = [0.08, 0, 0.04];

            x = grid.x;
            density = [x(:,1) + 2*x(:,3), ...
                       x(:,2) - x(:,3), x(:,1) - x(:,2)];

            % At fixed grid locations x', the actively rotated field is
            % q'(x') = R*q(R^T*x'). Rows therefore transform with R^T.
            xOriginal = x * R;
            densityOriginal = [xOriginal(:,1) + 2*xOriginal(:,3), ...
                               xOriginal(:,2) - xOriginal(:,3), ...
                               xOriginal(:,1) - xOriginal(:,2)];
            densityRotated = densityOriginal * R.';

            [qphi, qtheta] = quadest.errorest.DensityModifier.computeDensityAtRoots( ...
                grid, target, density);
            [qphiRot, qthetaRot] = quadest.errorest.DensityModifier.computeDensityAtRoots( ...
                grid, target * R.', densityRotated);

            obj.assertAlmostEqual(qphiRot, qphi * R.', 1e-11);
            obj.assertAlmostEqual(qthetaRot, qtheta * R.', 1e-11);
        end

        function test_interpolated_root_residual_diagnostic(obj)
            target = [0.06, 0, 0.04];
            phi = 0;
            theta0 = obj.geom.findThetaRoot(target, phi);

            [valid, residual] = quadest.util.Diagnostics.checkInterpolatedThetaRoots( ...
                obj.geom, theta0, phi, target, 1e-6);
            obj.assertTrue(valid);
            obj.assertTrue(residual <= 1e-6);

            mirroredTarget = target;
            mirroredTarget(3) = -mirroredTarget(3);
            [validMirrored, residualMirrored] = ...
                quadest.util.Diagnostics.checkInterpolatedThetaRoots( ...
                    obj.geom, pi - theta0, phi, mirroredTarget, 1e-6);
            obj.assertTrue(validMirrored);
            obj.assertTrue(residualMirrored <= 1e-6);

            warningId = 'quadest:Diagnostics:interpolatedRootResidual';
            oldWarning = warning('query', warningId);
            cleanup = onCleanup(@() warning(oldWarning.state, warningId)); %#ok<NASGU>
            warning('off', warningId);
            [validBad, residualBad] = quadest.util.Diagnostics.checkInterpolatedThetaRoots( ...
                obj.geom, theta0 + 0.05, phi, target, 1e-6);
            obj.assertFalse(validBad);
            obj.assertTrue(residualBad > 1e-6);
        end
        
        function test_uniformbuilder_create(obj)
            %TEST_UNIFORMBUILDER_CREATE Verify builder builds interpolants
            grid = quadest.grid.AxsymGrid(obj.geom, 'nth', 10, 'nph', 16);
            
            upsampFactors = [1, 2];
            config = quadest.util.Config();
            [interpolants, ~, ~] = quadest.errorest.UniformEstimateBuilder.build(...
                obj.geom, grid, obj.kernel, upsampFactors, false, config);
            
            obj.assertTrue(iscell(interpolants), 'Should return cell array');
        end
        
        function test_uniformbuilder_interpolants(obj)
            %TEST_UNIFORMBUILDER_INTERPOLANTS Verify interpolant properties
            grid = quadest.grid.AxsymGrid(obj.geom, 'nth', 10, 'nph', 16);
            
            upsampFactors = [1, 2];
            config = quadest.util.Config();
            [interpolants, ~, ~] = quadest.errorest.UniformEstimateBuilder.build(...
                obj.geom, grid, obj.kernel, upsampFactors, false, config);
            
            % Should have interpolants for each factor and component
            % The cell array is sized [max(upsampFac) x ncomp]
            obj.assertTrue(~isempty(interpolants{1,1}), 'Should have interpolant for factor 1');
            obj.assertTrue(~isempty(interpolants{2,1}), 'Should have interpolant for factor 2');
            
            % Each interpolant should be a griddedInterpolant
            obj.assertTrue(isa(interpolants{1,1}, 'griddedInterpolant'));
        end
        
        function test_estimator_creation(obj)
            %TEST_ESTIMATOR_CREATION Verify estimator construction
            estimator = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', [1, 2]);
            
            obj.assertTrue(isa(estimator, 'quadest.errorest.ErrorEstimator'));
            obj.assertEqual(estimator.grid.nth, 10);
            obj.assertEqual(estimator.grid.nph, 16);
        end
        
        function test_estimator_evaluate_farfield(obj)
            %TEST_ESTIMATOR_EVALUATE_FARFIELD Estimate for far target
            estimator = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:3);
            
            % Target far from surface
            target = [0.5, 0, 0];
            density = ones(10*16, 3);
            
            est = estimator.evaluate(target, density);
            
            obj.assertSize(est, [1, 1]);
            obj.assertTrue(est > 0, 'Estimate should be positive');
            obj.assertTrue(est < 1e-3, 'Far-field estimate should be small');
        end
        
        function test_estimator_evaluate_nearfield(obj)
            %TEST_ESTIMATOR_EVALUATE_NEARFIELD Estimate for near target
            estimator = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:3);
            
            % Target close to surface
            target = [0.051, 0, 0];  % Just outside equator (a=0.05)
            density = ones(10*16, 3);
            
            est = estimator.evaluate(target, density);
            
            obj.assertSize(est, [1, 1]);
            obj.assertTrue(est > 0, 'Estimate should be positive');
        end
        
        function test_estimator_evaluate_multiple_targets(obj)
            %TEST_ESTIMATOR_EVALUATE_MULTIPLE_TARGETS Multiple target evaluation
            estimator = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:2);
            
            targets = [0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5; 0.2, 0.2, 0.2];
            density = ones(10*16, 3);
            
            est = estimator.evaluate(targets, density);
            
            obj.assertSize(est, [4, 1]);
            obj.assertTrue(all(est > 0), 'All estimates should be positive');
        end
        
        function test_estimator_classification(obj)
            %TEST_ESTIMATOR_CLASSIFICATION Verify classification with tolerance
            estimator = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:2);
            
            % Mix of near and far targets
            targets = [0.051, 0, 0; 0.5, 0, 0];
            density = ones(10*16, 3);
            
            [est, classification] = estimator.evaluate(targets, density, 'tol', 1e-6);
            
            obj.assertTrue(isstruct(classification), 'Should return classification struct');
            obj.assertTrue(isfield(classification, 'upsampfac'), 'Should have upsampfac field');
            obj.assertTrue(isfield(classification, 'isSQ'), 'Should have isSQ field');
            
            obj.assertSize(classification.upsampfac, [2, 1]);
            obj.assertSize(classification.isSQ, [2, 1]);
            
            % Far target should need less upsampling
            obj.assertTrue(classification.upsampfac(2) <= classification.upsampfac(1), ...
                'Far target should need less upsampling');
        end
        
        function test_estimator_handles_axis(obj)
            %TEST_ESTIMATOR_HANDLES_AXIS Target on z-axis
            estimator = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:2);
            
            % Target on z-axis (rxy = 0)
            target = [0, 0, 0.2];
            density = ones(10*16, 3);
            
            % Should not error
            obj.assertNoError(@() estimator.evaluate(target, density));
        end
        
        function test_evaluateDirect_vs_evaluate(obj)
            %TEST_EVALUATEDIRECT_VS_EVALUATE Direct and tabulated methods agree for unit density
            estimator = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:2);
            
            % Exterior targets whose direct indicators remain above the
            % floating-point lower bound used by the tabulation.
            targets = [0.15, 0, 0; 0, 0.15, 0; 0, 0, 0.25; 0.12, 0.08, 0.15];
            density = ones(10*16, 3);
            
            est_tab = estimator.evaluate(targets, density);
            est_dir = estimator.evaluateDirect(targets, density);
            
            obj.assertSize(est_dir, [4, 1]);
            obj.assertTrue(all(est_dir > 0), 'Direct estimates should be positive');
            obj.assertTrue(all(isfinite(est_dir)), 'Direct estimates should be finite');
            
            % Both methods should agree within one order of magnitude for unit density
            ratio = log10(est_tab) - log10(est_dir);
            meaningful = est_dir > 1e-10;
            obj.assertTrue(any(meaningful), 'Expected meaningful direct indicators.');
            obj.assertTrue(all(abs(ratio(meaningful)) < 1), ...
                sprintf(['Tabulated and direct should be within 1 order of magnitude ' ...
                'above the numerical floor (max ratio: %.1f)'], ...
                max(abs(ratio(meaningful)))));
        end
        
        function test_evaluateDirect_smooth_density(obj)
            %TEST_EVALUATEDIRECT_SMOOTH_DENSITY Non-trivial density produces valid output
            estimator = quadest.errorest.ErrorEstimator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:2);
            
            grid = estimator.getGrid();
            [theta_mat, ~] = grid.meshgrid();
            
            % Smooth non-uniform density: cos(theta) pattern
            density_pattern = cos(theta_mat(:));
            density = [density_pattern, 0.5*density_pattern, 2*density_pattern];
            
            targets = [0.5, 0, 0; 0.2, 0.2, 0.2];
            
            est = estimator.evaluateDirect(targets, density);
            
            obj.assertSize(est, [2, 1]);
            obj.assertTrue(all(est >= 0), 'Estimates should be non-negative');
            obj.assertTrue(all(isfinite(est)), 'Estimates should be finite');
            
            % Should differ from unit density case
            est_unit = estimator.evaluateDirect(targets, ones(10*16, 3));
            obj.assertTrue(any(abs(est - est_unit) > 0), ...
                'Smooth density should differ from unit density');
        end

        function test_computeEstimates_public(obj)
            %TEST_COMPUTEESTIMATES_PUBLIC computeEstimates is accessible as public static
            grid = quadest.grid.AxsymGrid(obj.geom, 'nth', 10, 'nph', 16);
            p = obj.kernel.singularityOrder();
            
            targets = [0.5, 0, 0];
            [estimates, thetaRoots] = quadest.errorest.UniformEstimateBuilder.computeEstimates(...
                obj.geom, grid, obj.kernel, targets, p);
            
            obj.assertSize(estimates, [1, 3]);
            obj.assertSize(thetaRoots, [1, 1]);
            obj.assertTrue(all(estimates > 0), 'Estimates should be positive');
        end
    end
end
