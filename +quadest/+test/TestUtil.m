classdef TestUtil < quadest.test.TestBase
%TESTUTIL Unit tests for +quadest/+util module
%
%   Tests for Config, GaussLegendre, GaussLaguerre8, and Diagnostics

    methods
        function test_config_defaults(obj)
            %TEST_CONFIG_DEFAULTS Verify default configuration values
            cfg = quadest.util.Config();  % Create Config object
            
            obj.assertTrue(isa(cfg, 'quadest.util.Config'), 'Config should return Config object');
            obj.assertEqual(cfg.nth, 40);
            obj.assertEqual(cfg.nph, 60);
            obj.assertAlmostEqual(cfg.tol, 1e-6);
        end
        
        function test_config_get(obj)
            %TEST_CONFIG_GET Test Config with overrides
            cfg = quadest.util.Config('nth', 100);
            obj.assertEqual(cfg.nth, 100);
            
            % Default when not specified
            cfg2 = quadest.util.Config();
            obj.assertEqual(cfg2.nth, 40);
        end
        
        function test_gauss_legendre_nodes(obj)
            %TEST_GAUSS_LEGENDRE_NODES Verify GL nodes are in [-1, 1]
            [x, w] = quadest.util.GaussLegendre(10, -1, 1);
            
            obj.assertSize(x, [10, 1]);
            obj.assertSize(w, [10, 1]);
            obj.assertTrue(all(x >= -1 & x <= 1), 'Nodes should be in [-1, 1]');
            obj.assertTrue(all(w > 0), 'Weights should be positive');
        end
        
        function test_gauss_legendre_integration(obj)
            %TEST_GAUSS_LEGENDRE_INTEGRATION Verify integration accuracy
            % Integrate x^2 from -1 to 1 (exact = 2/3)
            [x, w] = quadest.util.GaussLegendre(5, -1, 1);
            integral = sum(w .* x.^2);
            obj.assertAlmostEqual(integral, 2/3, 1e-14);
            
            % Integrate cos(x) from 0 to pi (exact = 0)
            [x, w] = quadest.util.GaussLegendre(20, 0, pi);
            integral = sum(w .* cos(x));
            obj.assertAlmostEqual(integral, 0, 1e-10);
        end
        
        function test_gauss_legendre_symmetry(obj)
            %TEST_GAUSS_LEGENDRE_SYMMETRY Nodes and weights should be symmetric
            [x, w] = quadest.util.GaussLegendre(8, -1, 1);
            
            obj.assertAlmostEqual(x, -flipud(x), 1e-14, 'Nodes should be symmetric');
            obj.assertAlmostEqual(w, flipud(w), 1e-14, 'Weights should be symmetric');
        end
        
        function test_gauss_laguerre_hardcoded(obj)
            %TEST_GAUSS_LAGUERRE_HARDCODED Verify 8-point Gauss-Laguerre
            [x, w] = quadest.util.GaussLaguerre8();
            
            obj.assertSize(x, [1, 8]);  % Returns row vectors
            obj.assertSize(w, [1, 8]);
            obj.assertTrue(all(x > 0), 'Laguerre nodes should be positive');
            
            % Integrate exp(-x) from 0 to inf (exact = 1)
            % GL8 uses weights that include exp(x) factor
            % So we compute sum(w) which should integrate 1*exp(-x)
            integral = sum(w);
            obj.assertAlmostEqual(integral, 1, 1e-10);
        end
        
        function test_diagnostics_checkDensity(obj)
            %TEST_DIAGNOSTICS_CHECKDENSITY Verify density resolution check
            % Well-resolved density (constant)
            q = ones(40*60, 3);
            obj.assertNoError(@() quadest.util.Diagnostics.checkDensityResolution(q, 60));
            
            % Note: Under-resolved would issue warning, not error
        end
    end
end
