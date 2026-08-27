classdef QuadIndTestCase < matlab.unittest.TestCase
    methods
        function assertAlmostEqual(testCase, actual, expected, tolerance, diagnostic)
            if nargin < 4 || isempty(tolerance)
                tolerance = 1e-12;
            end
            if nargin < 5
                diagnostic = '';
            end
            testCase.assertEqual(isnan(actual), isnan(expected), ...
                'NaN locations differ between actual and expected values.');
            testCase.assertEqual(actual, expected, 'AbsTol', tolerance, diagnostic);
        end

        function assertNoError(~, functionHandle, ~)
            functionHandle();
        end
    end
end
