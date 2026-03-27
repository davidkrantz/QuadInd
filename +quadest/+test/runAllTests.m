function results = runAllTests(varargin)
%RUNALLTESTS Run all unit tests for the quadest package
%
%   results = runAllTests() runs all tests and returns results struct
%   results = runAllTests('verbose', true) enables verbose output
%   results = runAllTests('stopOnFailure', true) stops on first failure
%
%   Example:
%       results = quadest.test.runAllTests();
%       disp(results.summary);
%
%   See also: quadest.test.TestGeometry, quadest.test.TestGrid

    % Parse inputs
    p = inputParser;
    addParameter(p, 'verbose', true, @islogical);
    addParameter(p, 'stopOnFailure', false, @islogical);
    parse(p, varargin{:});
    opts = p.Results;
    
    % Initialize results
    results = struct();
    results.tests = {};
    results.passed = 0;
    results.failed = 0;
    results.errors = 0;
    results.skipped = 0;
    results.startTime = datetime('now');
    
    if opts.verbose
        fprintf('\n========================================\n');
        fprintf('  QuadEst Test Suite\n');
        fprintf('  %s\n', char(results.startTime));
        fprintf('========================================\n\n');
    end
    
    % List of test classes to run
    testClasses = {
        'quadest.test.TestUtil'
        'quadest.test.TestGeometry'
        'quadest.test.TestGrid'
        'quadest.test.TestKernel'
        'quadest.test.TestErrorEstimation'
    };
    
    % Run each test class
    for i = 1:length(testClasses)
        className = testClasses{i};
        
        if opts.verbose
            fprintf('Running %s...\n', className);
        end
        
        try
            % Create test instance and run
            testObj = feval(className);
            classResults = testObj.run('verbose', opts.verbose);
            
            % Aggregate results
            results.tests = [results.tests, classResults.tests];  % horizontal concatenation
            results.passed = results.passed + classResults.passed;
            results.failed = results.failed + classResults.failed;
            results.errors = results.errors + classResults.errors;
            results.skipped = results.skipped + classResults.skipped;
            
            if opts.stopOnFailure && (classResults.failed > 0 || classResults.errors > 0)
                if opts.verbose
                    fprintf('\nStopping on failure.\n');
                end
                break;
            end
            
        catch ME
            results.errors = results.errors + 1;
            results.tests = [results.tests, {struct('name', className, 'status', 'ERROR', ...
                'message', ME.message, 'time', 0)}];  % horizontal concatenation
            if opts.verbose
                fprintf('  ERROR loading test class: %s\n', ME.message);
            end
        end
        
        if opts.verbose
            fprintf('\n');
        end
    end
    
    % Compute summary
    results.endTime = datetime('now');
    results.totalTime = seconds(results.endTime - results.startTime);
    results.total = results.passed + results.failed + results.errors + results.skipped;
    
    % Generate summary string
    if results.failed == 0 && results.errors == 0
        statusStr = 'PASSED';
    else
        statusStr = 'FAILED';
    end
    
    results.summary = sprintf('%s: %d/%d tests passed (%.2f s)', ...
        statusStr, results.passed, results.total, results.totalTime);
    
    if opts.verbose
        fprintf('========================================\n');
        fprintf('  %s\n', results.summary);
        if results.failed > 0
            fprintf('  Failed: %d\n', results.failed);
        end
        if results.errors > 0
            fprintf('  Errors: %d\n', results.errors);
        end
        if results.skipped > 0
            fprintf('  Skipped: %d\n', results.skipped);
        end
        fprintf('========================================\n\n');
    end
end
