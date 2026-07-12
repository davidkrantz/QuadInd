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
            
            % Far-field targets where estimates are well-behaved
            targets = [0.5, 0, 0; 0, 0.5, 0; 0, 0, 0.5; 0.2, 0.2, 0.2];
            density = ones(10*16, 3);
            
            est_tab = estimator.evaluate(targets, density);
            est_dir = estimator.evaluateDirect(targets, density);
            
            obj.assertSize(est_dir, [4, 1]);
            obj.assertTrue(all(est_dir > 0), 'Direct estimates should be positive');
            obj.assertTrue(all(isfinite(est_dir)), 'Direct estimates should be finite');
            
            % Both methods should agree within ~1 order of magnitude for unit density
            ratio = log10(est_tab) - log10(est_dir);
            obj.assertTrue(all(abs(ratio) < 2), ...
                sprintf('Tabulated and direct should be within 2 orders of magnitude (max ratio: %.1f)', max(abs(ratio))));
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
