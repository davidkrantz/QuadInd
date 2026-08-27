classdef TestKernel < QuadIndTestCase
%TESTKERNEL Unit tests for +quadind/+kernel module
%
%   Tests for Kernel and StokesStresslet

    properties (Access = private)
        kernel
    end
    
    methods (TestMethodSetup)
        function setUp(obj)
            %SETUP Create test fixtures
            obj.kernel = quadind.kernel.StokesStresslet();
        end
    end

    methods (Test)
        function test_stresslet_creation(obj)
            %TEST_STRESSLET_CREATION Verify stresslet construction
            k = quadind.kernel.StokesStresslet();
            obj.assertTrue(isa(k, 'quadind.kernel.StokesStresslet'));
            obj.assertTrue(isa(k, 'quadind.kernel.Kernel'));
        end
        
        function test_stresslet_singularity_order(obj)
            %TEST_STRESSLET_SINGULARITY_ORDER Verify p = 5/2
            obj.assertEqual(obj.kernel.singularityOrder(), 2.5);
        end
        
        function test_stresslet_components(obj)
            %TEST_STRESSLET_COMPONENTS Verify 3 velocity components
            obj.assertEqual(obj.kernel.numComponents(), 3);
        end
        
        function test_stresslet_evaluate_single(obj)
            %TEST_STRESSLET_EVALUATE_SINGLE Evaluate at single target
            % Source at origin, target at (1, 0, 0)
            x_src = [0, 0, 0];
            n_src = [1, 0, 0];  % Normal pointing in x
            x_tgt = [1, 0, 0];
            q = [1, 0, 0];  % Density in x direction
            w = 1;  % Unit weight
            
            u = obj.kernel.evaluate(x_tgt, x_src, n_src, q, w);
            
            obj.assertSize(u, [1, 3]);
            
            % Manual calculation:
            % r = x_tgt - x_src = [1, 0, 0]
            % |r| = 1
            % r_dot_n = 1
            % r_dot_q = 1
            % u_j = -6 * q_i * r_i * r_j * r_k * n_k / |r|^5
            %     = -6 * 1 * 1 * r_j * 1 / 1
            %     = -6 * [1, 0, 0]
            expected = -6 * [1, 0, 0];
            obj.assertAlmostEqual(u, expected, 1e-10);
        end

        function test_stresslet_multiple_sources(obj)
            %TEST_STRESSLET_MULTIPLE_SOURCES Evaluate with multiple sources
            x_src = [0, 0, 0; 1, 0, 0; 0, 1, 0];  % 3 sources
            n_src = [0, 0, 1; 0, 0, 1; 0, 0, 1];
            q = [1, 0, 0; 0, 1, 0; 0, 0, 1];      % 3 densities
            w = ones(3, 1);  % Unit weights
            
            x_tgt = [2, 2, 2];  % Single target
            
            u = obj.kernel.evaluate(x_tgt, x_src, n_src, q, w);
            
            obj.assertSize(u, [1, 3]);
        end
        
        function test_stresslet_multiple_targets(obj)
            %TEST_STRESSLET_MULTIPLE_TARGETS Evaluate at multiple targets
            x_src = [0, 0, 0];
            n_src = [0, 0, 1];
            q = [1, 0, 0];
            w = 1;
            
            x_tgt = [1, 0, 0; 0, 1, 0; 0, 0, 1; 2, 2, 2];  % 4 targets
            
            u = obj.kernel.evaluate(x_tgt, x_src, n_src, q, w);
            
            obj.assertSize(u, [4, 3]);
        end
        
        function test_stresslet_symmetry(obj)
            %TEST_STRESSLET_SYMMETRY Verify symmetry properties
            x_src = [0, 0, 0];
            n_src = [0, 0, 1];
            q = [1, 0, 0];
            w = 1;
            
            % Two targets symmetric about z-axis
            u1 = obj.kernel.evaluate([1, 0, 0], x_src, n_src, q, w);
            u2 = obj.kernel.evaluate([-1, 0, 0], x_src, n_src, q, w);
            
            % u_x should have opposite signs, u_y same, u_z same
            obj.assertAlmostEqual(u1(1), -u2(1), 1e-10);
        end
        
        function test_stresslet_with_weights(obj)
            %TEST_STRESSLET_WITH_WEIGHTS Verify weighted evaluation
            x_src = [0, 0, 0; 1, 0, 0];
            n_src = [0, 0, 1; 0, 0, 1];
            q = [1, 0, 0; 1, 0, 0];
            w = [0.5; 0.5];
            
            x_tgt = [2, 2, 2];
            
            u = obj.kernel.evaluate(x_tgt, x_src, n_src, q, w);
            
            obj.assertSize(u, [1, 3]);
        end
        
        function test_stresslet_zero_density(obj)
            %TEST_STRESSLET_ZERO_DENSITY Zero density gives zero velocity
            x_src = [0, 0, 0];
            n_src = [0, 0, 1];
            q = [0, 0, 0];
            w = 1;
            
            u = obj.kernel.evaluate([1, 0, 0], x_src, n_src, q, w);
            
            obj.assertAlmostEqual(u, [0, 0, 0], 1e-14);
        end
    end
end
