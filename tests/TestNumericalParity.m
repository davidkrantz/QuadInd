classdef TestNumericalParity < matlab.unittest.TestCase
    % Regression values generated with previous code version.
    methods (Test)
        function test_reference_results(testCase)
            geometry = quadind.geometry.Spheroid('a', 0.05, 'c', 0.1);
            kernel = quadind.kernel.StokesStresslet();
            config = quadind.util.Config('ntab', 8, 'nztab', 8);
            evaluator = quadind.IndicatorEvaluator(geometry, kernel, ...
                'nth', 8, 'nph', 12, 'upsampFactors', 1:2, 'config', config);
            grid = evaluator.getGrid();
            density = [1 + grid.x(:,1) + 0.5*grid.x(:,3), ...
                2 + grid.x(:,2) - 0.25*grid.x(:,3), ...
                0.5 + grid.x(:,1).*grid.x(:,2) + grid.x(:,3).^2];
            targets = [0.12,0,0.04; 0.08*cos(0.7),0.08*sin(0.7),0.04; 1,0.1,0.2];

            [indicator, classification] = evaluator.evaluate(...
                targets, density, 'tol', 1e-6);
            directIndicator = evaluator.evaluateDirect(targets, density);
            potential = kernel.evaluateOnGrid(targets, grid, density);

            expectedIndicator = [0.22094725865433557; ...
                3.9351058459220134; 2.11669788592363e-09];
            expectedDirectIndicator = [0.21633251960760247; ...
                6.2964097728816197; 1.5143326808568707e-07];
            expectedPotential = [ ...
                -0.46446902534572942, -0.0091871402226986823, -0.039654841421988357; ...
                -3.0995856244027093, -3.6081356007284864, -2.6749809911835767; ...
                -0.0061152009426978481, -0.00061080886978322722, -0.0012174938013666909];
            expectedFactorIndicators = {expectedIndicator, [ ...
                0.00093511041461729065; 0.25021206824970149; NaN]};

            testCase.verifyEqual(indicator, expectedIndicator, 'RelTol', 5e-14);
            testCase.verifyEqual(directIndicator, expectedDirectIndicator, 'RelTol', 5e-14);
            testCase.verifyEqual(potential, expectedPotential, 'RelTol', 5e-14);
            for factorIndex = 1:numel(expectedFactorIndicators)
                testCase.verifyEqual(classification.indicators{factorIndex}, ...
                    expectedFactorIndicators{factorIndex}, 'RelTol', 5e-14);
            end
            testCase.verifyEqual(classification.upsamplingFactor, [3;3;1]);
            testCase.verifyEqual(classification.requiresSpecialQuadrature, ...
                [true;true;false]);
        end
    end
end
