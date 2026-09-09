classdef TestIndicatorEvaluation < QuadIndTestCase
%TESTINDICATOREVALUATION Tests for indicator construction and evaluation
%
%   Tests for IndicatorEvaluator, RootFinder, DensityModifier, UniformIndicatorBuilder

    properties (Access = private)
        geom
        kernel
    end
    
    methods (TestMethodSetup)
        function setUp(obj)
            %SETUP Create test fixtures
            obj.geom = quadind.geometry.Spheroid('a', 0.05, 'c', 0.1);
            obj.kernel = quadind.kernel.StokesStresslet();
        end
    end

    methods (Test)
        function test_densitymodifier_matches_local_linear_continuation(obj)
            grid = quadind.grid.AxsymGrid(obj.geom, 'nth', 18, 'nph', 32);
            [theta, phi] = grid.meshgrid();
            density = [phi(:), theta(:), 2*phi(:) - 3*theta(:)];
            alpha = 0.7;
            target = [0.09*cos(alpha), 0.09*sin(alpha), 0.04];

            [qphi, qtheta, phi0] = ...
                quadind.indicator.DensityModifier.computeDensityAtRoots( ...
                    grid, target, density);
            [itheta, iphi] = ...
                quadind.indicator.UniformIndicatorBuilder.findNearestNodes(grid, target);
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
                grid = quadind.grid.AxsymGrid(obj.geom, ...
                    'nth', resolutions(k,1), 'nph', resolutions(k,2));
                [theta, phi] = grid.meshgrid();
                density = [cos(2*phi(:)), cos(2*theta(:)), ...
                    cos(2*phi(:)).*cos(2*theta(:))];
                [qphi, qtheta, phi0] = ...
                    quadind.indicator.DensityModifier.computeDensityAtRoots( ...
                        grid, target, density);
                [~, iphi] = ...
                    quadind.indicator.UniformIndicatorBuilder.findNearestNodes(grid, target);
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
            [phi0, Gphi] = quadind.indicator.RootFinder.computePhiRoot( ...
                obj.geom, target, thetaStar);
            rphi = obj.geom.evaluate(thetaStar, phi0) - target;
            expectedPhi = 2 * sum(rphi .* obj.geom.drdphi(thetaStar, phi0), 2);
            obj.assertAlmostEqual(Gphi, expectedPhi, 1e-12);

            theta0 = obj.geom.findThetaRoot(target, 0);
            [~, Gtheta] = quadind.indicator.RootFinder.computeThetaRoot( ...
                obj.geom, target, 0);
            rtheta = obj.geom.evaluate(theta0, 0) - target;
            expectedTheta = pi * sum(rtheta .* obj.geom.drdtheta(theta0, 0), 2);
            obj.assertAlmostEqual(Gtheta, expectedTheta, 1e-12);
        end

        function test_spheroid_roots_satisfy_analytic_equation(obj)
            geometries = {obj.geom, ...
                quadind.geometry.Spheroid('a', 0.1, 'c', 0.04), ...
                quadind.geometry.Spheroid('a', 0.08, 'c', 0.08)};
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
            grid = quadind.grid.AxsymGrid(obj.geom, 'nth', 12, 'nph', 16);
            density = ones(grid.numPoints(), 3);
            [qphi, qtheta, phi0] = ...
                quadind.indicator.DensityModifier.computeDensityAtRoots( ...
                    grid, [0, 0, 0.2], density);
            obj.assertAlmostEqual(qphi, zeros(1,3), 0);
            obj.assertAlmostEqual(qtheta, ones(1,3), 1e-12);
            obj.assertTrue(isinf(imag(phi0)));

            indicator = quadind.indicator.UniformIndicatorBuilder.computeDirectIndicators( ...
                obj.geom, grid, obj.kernel, [0, 0, 0.2], density);
            obj.assertTrue(isfinite(indicator) && indicator >= 0);
        end

        function test_evaluator_rejects_wrong_density_component_count(obj)
            cfg = quadind.util.Config('ntab', 12, 'nztab', 12);
            evaluator = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1, 'config', cfg);
            density = ones(10*16, 2);
            obj.assertError(@() evaluator.evaluate([0.2, 0, 0], density), ...
                'quadind:IndicatorEvaluator:invalidDensity');
        end

        function test_tabulation_safety_factor_scales_indicator(obj)
            cfg1 = quadind.util.Config('ntab', 20, 'nztab', 20, ...
                'tabulationSafetyFactor', 1);
            cfg2 = quadind.util.Config('ntab', 20, 'nztab', 20, ...
                'tabulationSafetyFactor', 2);
            evaluator1 = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1, 'config', cfg1);
            evaluator2 = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1, 'config', cfg2);
            target = [0.12, 0.03, 0.04];
            density = ones(10*16, 3);
            indicator1 = evaluator1.evaluate(target, density);
            indicator2 = evaluator2.evaluate(target, density);
            obj.assertAlmostEqual(indicator2, 2*indicator1, ...
                1e-12*max(1, indicator2));
        end

        function test_interpolants_use_linear_extrapolation(obj)
            evaluator = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 8, 'nph', 12, 'upsampFactors', 1:2, ...
                'config', quadind.util.Config('ntab', 8, 'nztab', 8));

            for factor = evaluator.upsampFactors
                for component = 1:obj.kernel.numComponents()
                    obj.assertEqual(evaluator.interpolants{factor,component}.Method, ...
                        'linear');
                    obj.assertEqual(evaluator.interpolants{factor,component}.ExtrapolationMethod, ...
                        'linear');
                end
                obj.assertEqual(evaluator.rootInterpolants.Re{factor}.ExtrapolationMethod, ...
                    'linear');
                obj.assertEqual(evaluator.rootInterpolants.Im{factor}.ExtrapolationMethod, ...
                    'linear');
            end
        end

        function test_log_linear_extrapolation_across_all_boundaries(obj)
            cfg = quadind.util.Config('ntab', 8, 'nztab', 8);
            evaluator = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 8, 'nph', 12, 'upsampFactors', 1:2, ...
                'interpolateRoots', false, 'config', cfg);
            grid = evaluator.getGrid();
            extent = 2 * (obj.geom.maxRadius() + obj.geom.maxHeight());
            h = extent / (cfg.ntab - 2);
            targets = [ ...
                extent-h, 0, extent-h; ... % corner: inner-inner
                extent,   0, extent-h; ... % corner: edge-inner
                extent-h, 0, extent;   ... % corner: inner-edge
                extent,   0, extent;   ... % corner: edge-edge
                extent+h, 0, extent+h; ... % corner: extrapolated
                extent-h, 0, 0.1; ...      % radial triplet
                extent,   0, 0.1; ...
                extent+h, 0, 0.1; ...
                0.1, 0, -(extent-h); ...   % negative-z triplet
                0.1, 0, -extent; ...
                0.1, 0, -(extent+h)];

            for component = 1:obj.kernel.numComponents()
                density = zeros(grid.numPoints(), obj.kernel.numComponents());
                density(:,component) = 1;
                [~, classification] = evaluator.evaluate( ...
                    targets, density, 'tol', realmin);

                for factorIndex = 1:numel(classification.indicators)
                    values = classification.indicators{factorIndex};
                    logs = log10(values);
                    obj.verifyEqual(logs(8), 2*logs(7) - logs(6), ...
                        'AbsTol', 5e-12);
                    obj.verifyEqual(logs(11), 2*logs(10) - logs(9), ...
                        'AbsTol', 5e-12);
                    expectedCorner = logs(1) - 2*logs(2) - 2*logs(3) + 4*logs(4);
                    obj.verifyEqual(logs(5), expectedCorner, 'AbsTol', 5e-12);
                end
            end
        end

        function test_near_boundary_extrapolation_matches_direct(obj)
            cfg = quadind.util.Config('ntab', 16, 'nztab', 16);
            evaluator = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1, ...
                'interpolateRoots', false, 'config', cfg);
            grid = evaluator.getGrid();
            density = ones(grid.numPoints(), obj.kernel.numComponents());
            extent = 2 * (obj.geom.maxRadius() + obj.geom.maxHeight());
            halfCell = 0.5 * extent / (cfg.ntab - 2);
            targets = [extent+halfCell, 0, 0.1; ...
                0.1, 0, extent+halfCell; ...
                extent+halfCell, 0, extent+halfCell];

            extrapolated = evaluator.evaluate(targets, density);
            direct = evaluator.evaluateDirect(targets, density);
            obj.verifyTrue(all(abs(log10(extrapolated ./ direct)) < 0.25));
        end

        function test_extrapolated_root_residual_fallback(obj)
            base = {'ntab', 8, 'nztab', 8};
            guarded = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 8, 'nph', 12, 'upsampFactors', 1, ...
                'config', quadind.util.Config(base{:}, ...
                    'interpolatedRootResidualTolerance', 1e-6));
            directRoots = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 8, 'nph', 12, 'upsampFactors', 1, ...
                'interpolateRoots', false, 'config', quadind.util.Config(base{:}));
            looseGuard = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 8, 'nph', 12, 'upsampFactors', 1, ...
                'config', quadind.util.Config(base{:}, ...
                    'interpolatedRootResidualTolerance', 1));
            rawRoots = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 8, 'nph', 12, 'upsampFactors', 1, ...
                'config', quadind.util.Config(base{:}, ...
                    'interpolatedRootResidualTolerance', NaN));

            grid = guarded.getGrid();
            [theta, phi] = grid.meshgrid();
            density = [theta(:), 2*theta(:)-phi(:), cos(theta(:))+sin(phi(:))];
            target = [0.5, 0, 0.1];

            guardedValue = guarded.evaluate(target, density);
            directRootValue = directRoots.evaluate(target, density);
            looseValue = looseGuard.evaluate(target, density);
            rawValue = rawRoots.evaluate(target, density);

            obj.verifyEqual(guardedValue, directRootValue, 'RelTol', 5e-13);
            obj.verifyEqual(looseValue, rawValue, 'RelTol', 5e-13);
            obj.verifyGreaterThan(abs(log10(guardedValue / rawValue)), 1e-4);
        end

        function test_densitymodifier_rotates_to_target_frame(obj)
            % Root densities are global Cartesian vectors; modifiers use
            % radial, azimuthal, and axial components at each target.
            qphi = [1+2i, 3-1i, 4; 2, -1i, 5];
            qtheta = [2-1i, -1+3i, 1; -3i, 4, -2];
            targets = [0, 2, 1; -1, 0, -1];

            modifier = quadind.indicator.DensityModifier.computeModifier( ...
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

        function test_root_interpolation_rotation_covariance(obj)
            % A grid-aligned active rotation of the target and vector field
            % must commute with both complex-root interpolations.
            grid = quadind.grid.AxsymGrid(obj.geom, 'nth', 10, 'nph', 16);
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

            [qphi, qtheta] = quadind.indicator.DensityModifier.computeDensityAtRoots( ...
                grid, target, density);
            [qphiRot, qthetaRot] = quadind.indicator.DensityModifier.computeDensityAtRoots( ...
                grid, target * R.', densityRotated);

            obj.assertAlmostEqual(qphiRot, qphi * R.', 1e-11);
            obj.assertAlmostEqual(qthetaRot, qtheta * R.', 1e-11);
        end

        function test_evaluator_classification(obj)
            %TEST_EVALUATOR_CLASSIFICATION Verify classification with tolerance
            evaluator = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:2, ...
                'config', quadind.util.Config('ntab', 12, 'nztab', 12));
            
            % Mix of near and far targets
            targets = [0.051, 0, 0; 0.5, 0, 0];
            density = ones(10*16, 3);
            
            [~, classification] = evaluator.evaluate(targets, density, 'tol', 1e-6);
            
            obj.assertTrue(isstruct(classification), 'Should return classification struct');
            obj.assertTrue(isfield(classification, 'upsamplingFactor'));
            obj.assertTrue(isfield(classification, 'requiresSpecialQuadrature'));
            
            obj.assertSize(classification.upsamplingFactor, [2, 1]);
            obj.assertSize(classification.requiresSpecialQuadrature, [2, 1]);
            
            % The far target should be accepted without special quadrature;
            % the near target may require a higher factor or special quadrature.
            obj.assertFalse(classification.requiresSpecialQuadrature(2));
            obj.assertTrue(classification.requiresSpecialQuadrature(1) || ...
                classification.upsamplingFactor(2) <= classification.upsamplingFactor(1));
        end

        function test_evaluator_classification_stops_after_acceptance(obj)
            evaluator = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 8, 'nph', 12, 'upsampFactors', 1:3, ...
                'config', quadind.util.Config('ntab', 8, 'nztab', 8));
            density = ones(8*12, 3);

            [indicators, classification] = evaluator.evaluate( ...
                [0.5, 0, 0; 0, 0, 0.5], density, 'tol', realmax);

            obj.assertEqual(classification.indicators{1}, indicators);
            obj.assertTrue(all(classification.masks{1}));
            obj.assertTrue(all(isnan(classification.indicators{2})));
            obj.assertTrue(all(isnan(classification.indicators{3})));
            obj.assertFalse(any(classification.masks{2}));
            obj.assertFalse(any(classification.masks{3}));
            obj.assertFalse(any(classification.requiresSpecialQuadrature));
        end
        
        function test_evaluator_handles_axis(obj)
            %TEST_EVALUATOR_HANDLES_AXIS Target on z-axis
            evaluator = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:2, ...
                'config', quadind.util.Config('ntab', 12, 'nztab', 12));
            
            % Target on z-axis (rxy = 0)
            target = [0, 0, 0.2];
            density = ones(10*16, 3);
            
            % Should not error
            obj.assertNoError(@() evaluator.evaluate(target, density));
        end
        
        function test_evaluateDirect_vs_evaluate(obj)
            %TEST_EVALUATEDIRECT_VS_EVALUATE Direct and tabulated methods agree for unit density
            evaluator = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 10, 'nph', 16, 'upsampFactors', 1:2, ...
                'config', quadind.util.Config('ntab', 16, 'nztab', 16));
            
            % Exterior targets whose direct indicators remain above the
            % floating-point lower bound used by the tabulation.
            targets = [0.15, 0, 0; 0, 0.15, 0; 0, 0, 0.25; 0.12, 0.08, 0.15];
            density = ones(10*16, 3);
            
            est_tab = evaluator.evaluate(targets, density);
            est_dir = evaluator.evaluateDirect(targets, density);
            
            obj.assertSize(est_dir, [4, 1]);
            obj.assertTrue(all(est_dir > 0), 'Direct indicators should be positive');
            obj.assertTrue(all(isfinite(est_dir)), 'Direct indicators should be finite');
            
            % Both methods should agree within one order of magnitude for unit density
            ratio = log10(est_tab) - log10(est_dir);
            meaningful = est_dir > 1e-10;
            obj.assertTrue(any(meaningful), 'Expected meaningful direct indicators.');
            obj.assertTrue(all(abs(ratio(meaningful)) < 1), ...
                sprintf(['Tabulated and direct should be within 1 order of magnitude ' ...
                'above the numerical floor (max ratio: %.1f)'], ...
                max(abs(ratio(meaningful)))));
        end
        
        function test_public_api_quick_start(obj)
            config = quadind.util.Config('ntab', 10, 'nztab', 10);
            evaluator = quadind.IndicatorEvaluator(obj.geom, obj.kernel, ...
                'nth', 8, 'nph', 12, 'upsampFactors', 1:2, 'config', config);
            grid = evaluator.getGrid();
            density = ones(grid.numPoints(), 3);
            [indicator, classification] = evaluator.evaluate(...
                [0.12, 0, 0.04], density, 'tol', 1e-6);
            obj.assertTrue(isfinite(indicator) && indicator >= 0);
            obj.assertTrue(isfield(classification, 'indicators'));
            obj.assertTrue(isfield(classification, 'upsamplingFactor'));
            obj.assertTrue(isfield(classification, 'requiresSpecialQuadrature'));
        end

        function test_asymmetric_geometry_is_rejected(obj)
            geometry = quadind.geometry.CustomAxsym(...
                @(theta) sin(theta), @(theta) cos(theta) + 0.05*sin(theta), ...
                @(theta) cos(theta), @(theta) -sin(theta) + 0.05*cos(theta));
            obj.assertError(@() quadind.IndicatorEvaluator(geometry, obj.kernel, ...
                'config', quadind.util.Config('ntab', 8, 'nztab', 8)), ...
                'quadind:IndicatorEvaluator:asymmetricGeometry');
        end
    end
end
