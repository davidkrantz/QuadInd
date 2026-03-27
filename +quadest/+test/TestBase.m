classdef TestBase < handle
%TESTBASE Base class for unit tests in the quadest package
%
%   Provides common testing infrastructure including assertions,
%   test discovery, and result aggregation.
%
%   Subclasses should define methods starting with 'test_' which
%   will be automatically discovered and run.
%
%   Example:
%       classdef TestMyModule < quadest.test.TestBase
%           methods
%               function test_something(obj)
%                   obj.assertEqual(1+1, 2);
%               end
%           end
%       end

    properties (Access = protected)
        tolerance = 1e-12;  % Default tolerance for floating point comparisons
    end
    
    methods
        function results = run(obj, varargin)
            %RUN Execute all test methods in this class
            
            p = inputParser;
            addParameter(p, 'verbose', true, @islogical);
            parse(p, varargin{:});
            opts = p.Results;
            
            % Initialize results
            results = struct();
            results.tests = {};
            results.passed = 0;
            results.failed = 0;
            results.errors = 0;
            results.skipped = 0;
            
            % Get test methods
            testMethods = obj.getTestMethods();
            
            % Run each test
            for i = 1:length(testMethods)
                methodName = testMethods{i};
                testResult = struct();
                testResult.name = sprintf('%s.%s', class(obj), methodName);
                
                try
                    % Setup
                    if ismethod(obj, 'setUp')
                        obj.setUp();
                    end
                    
                    % Run test
                    tic;
                    obj.(methodName)();
                    testResult.time = toc;
                    testResult.status = 'PASSED';
                    testResult.message = '';
                    results.passed = results.passed + 1;
                    
                    if opts.verbose
                        fprintf('  ✓ %s (%.3fs)\n', methodName, testResult.time);
                    end
                    
                catch ME
                    testResult.time = toc;
                    
                    if contains(ME.identifier, 'TestBase:AssertionFailed')
                        testResult.status = 'FAILED';
                        results.failed = results.failed + 1;
                    elseif contains(ME.identifier, 'TestBase:Skipped')
                        testResult.status = 'SKIPPED';
                        results.skipped = results.skipped + 1;
                    else
                        testResult.status = 'ERROR';
                        results.errors = results.errors + 1;
                    end
                    
                    testResult.message = ME.message;
                    
                    if opts.verbose
                        if strcmp(testResult.status, 'SKIPPED')
                            fprintf('  ○ %s (skipped: %s)\n', methodName, ME.message);
                        else
                            fprintf('  ✗ %s [%s]: %s\n', methodName, testResult.status, ME.message);
                        end
                    end
                end
                
                % Teardown
                try
                    if ismethod(obj, 'tearDown')
                        obj.tearDown();
                    end
                catch
                    % Ignore teardown errors
                end
                
                results.tests{end+1} = testResult;
            end
        end
    end
    
    methods (Access = protected)
        function assertEqual(obj, actual, expected, msg)
            %ASSERTEQUAL Assert that actual equals expected
            if nargin < 4
                msg = sprintf('Expected %s but got %s', mat2str(expected), mat2str(actual));
            end
            
            if isnumeric(actual) && isnumeric(expected)
                if ~isequal(size(actual), size(expected))
                    error('TestBase:AssertionFailed', ...
                        'Size mismatch: expected %s, got %s', ...
                        mat2str(size(expected)), mat2str(size(actual)));
                end
                if any(abs(actual(:) - expected(:)) > obj.tolerance)
                    error('TestBase:AssertionFailed', msg);
                end
            elseif ~isequal(actual, expected)
                error('TestBase:AssertionFailed', msg);
            end
        end
        
        function assertAlmostEqual(obj, actual, expected, tol, msg)
            %ASSERTALMOSTEQUAL Assert values are equal within tolerance
            if nargin < 4 || isempty(tol)
                tol = obj.tolerance;
            end
            if nargin < 5
                msg = sprintf('Values differ by more than %g', tol);
            end
            
            diff = abs(actual - expected);
            if any(diff(:) > tol)
                maxDiff = max(diff(:));
                error('TestBase:AssertionFailed', '%s (max diff: %g)', msg, maxDiff);
            end
        end
        
        function assertTrue(obj, condition, msg)
            %ASSERTTRUE Assert that condition is true
            if nargin < 3
                msg = 'Condition is not true';
            end
            if ~condition
                error('TestBase:AssertionFailed', msg);
            end
        end
        
        function assertFalse(obj, condition, msg)
            %ASSERTFALSE Assert that condition is false
            if nargin < 3
                msg = 'Condition is not false';
            end
            if condition
                error('TestBase:AssertionFailed', msg);
            end
        end
        
        function assertSize(obj, actual, expectedSize, msg)
            %ASSERTSIZE Assert array has expected size
            if nargin < 4
                msg = sprintf('Expected size %s but got %s', ...
                    mat2str(expectedSize), mat2str(size(actual)));
            end
            if ~isequal(size(actual), expectedSize)
                error('TestBase:AssertionFailed', msg);
            end
        end
        
        function assertError(obj, func, expectedId)
            %ASSERTERROR Assert that function throws an error
            try
                func();
                error('TestBase:AssertionFailed', 'Expected error was not thrown');
            catch ME
                if nargin >= 3 && ~isempty(expectedId)
                    if ~contains(ME.identifier, expectedId)
                        error('TestBase:AssertionFailed', ...
                            'Expected error ID containing "%s" but got "%s"', ...
                            expectedId, ME.identifier);
                    end
                end
                % Error was thrown as expected
            end
        end
        
        function assertNoError(obj, func, msg)
            %ASSERTNOERROR Assert that function does not throw
            if nargin < 3
                msg = 'Unexpected error';
            end
            try
                func();
            catch ME
                error('TestBase:AssertionFailed', '%s: %s', msg, ME.message);
            end
        end
        
        function skip(obj, reason)
            %SKIP Skip this test with given reason
            error('TestBase:Skipped', reason);
        end
    end
    
    methods (Access = private)
        function methods = getTestMethods(obj)
            %GETTESTMETHODS Get all methods starting with 'test_'
            meta = metaclass(obj);
            allMethods = {meta.MethodList.Name};
            methods = allMethods(startsWith(allMethods, 'test_'));
        end
    end
end
