classdef TestGrid < QuadIndTestCase
%TESTGRID Unit tests for +quadind/+grid module
%
%   Tests for AxsymGrid and Upsampler

    properties (Access = private)
        geom
        grid
    end
    
    methods (TestMethodSetup)
        function setUp(obj)
            %SETUP Create test fixtures
            obj.geom = quadind.geometry.Spheroid('a', 1, 'c', 2);
            obj.grid = quadind.grid.AxsymGrid(obj.geom, 'nth', 10, 'nph', 16);
        end
    end

    methods (Test)
        function test_grid_creation(obj)
            %TEST_GRID_CREATION Verify grid construction
            obj.assertEqual(obj.grid.nth, 10);
            obj.assertEqual(obj.grid.nph, 16);
            obj.assertEqual(obj.grid.numPoints(), 160);
        end
        
        function test_grid_points(obj)
            %TEST_GRID_POINTS Verify grid point coordinates
            x = obj.grid.x;
            
            obj.assertSize(x, [160, 3]);
            
            % All points should be on the spheroid surface
            % (x/a)^2 + (y/a)^2 + (z/c)^2 = 1
            a = obj.geom.a;
            c = obj.geom.c;
            surface_eq = (x(:,1).^2 + x(:,2).^2) / a^2 + x(:,3).^2 / c^2;
            obj.assertAlmostEqual(surface_eq, ones(160, 1), 1e-10);
        end
        
        function test_grid_normals(obj)
            %TEST_GRID_NORMALS Verify normals are unit vectors pointing outward
            n = obj.grid.n;
            
            obj.assertSize(n, [160, 3]);
            
            % All normals should be unit vectors
            n_mag = sqrt(sum(n.^2, 2));
            obj.assertAlmostEqual(n_mag, ones(160, 1), 1e-10);
            
            % For spheroid, normals should point outward (positive dot with position)
            x = obj.grid.x;
            dot_prod = sum(x .* n, 2);
            obj.assertTrue(all(dot_prod > 0), 'Normals should point outward');
        end
        
        function test_grid_weights(obj)
            %TEST_GRID_WEIGHTS Verify quadrature weights are positive
            w = obj.grid.w;
            
            obj.assertSize(w, [160, 1]);
            obj.assertTrue(all(w > 0), 'Weights should be positive');
            
            % Sum of weights should approximate surface area
            % Surface area of prolate spheroid with a=1, c=2
            e = sqrt(1 - (1/2)^2);  % eccentricity
            expected_area = 2*pi*1^2 * (1 + (2/1)/(e) * asin(e));
            computed_area = sum(w);
            
            % Should be within 1% for 10x16 grid
            obj.assertAlmostEqual(computed_area, expected_area, 0.01 * expected_area);
        end
        
        function test_grid_theta_phi(obj)
            %TEST_GRID_THETA_PHI Verify theta/phi arrays
            theta = obj.grid.theta;
            phi = obj.grid.phi;
            
            obj.assertSize(theta, [10, 1]);
            obj.assertSize(phi, [1, 16]);  % phi is row vector
            
            % Theta should be in (0, pi) - GL nodes
            obj.assertTrue(all(theta > 0 & theta < pi));
            
            % Phi should be in [0, 2*pi)
            obj.assertTrue(all(phi >= 0 & phi < 2*pi));
            
            % Phi should be uniformly spaced
            dphi = diff(phi);
            obj.assertAlmostEqual(dphi, dphi(1) * ones(1, 15), 1e-14);
        end
        
        function test_upsampler_factor2(obj)
            %TEST_UPSAMPLER_FACTOR2 Verify 2x upsampling
            targetGrid = obj.grid.upsample(2);  % Create 2x upsampled grid
            upsampler = quadind.grid.Upsampler(obj.grid, targetGrid);
            
            % Upsample a density and check target grid size
            q = ones(160, 3);  % Constant density
            q_up = upsampler.apply(q);
            
            % Should have 4x as many points
            obj.assertSize(q_up, [640, 3]);
            
            % Get target grid points
            x_up = targetGrid.x;
            
            % Points should still be on surface
            a = obj.geom.a;
            c = obj.geom.c;
            surface_eq = (x_up(:,1).^2 + x_up(:,2).^2) / a^2 + x_up(:,3).^2 / c^2;
            obj.assertAlmostEqual(surface_eq, ones(640, 1), 1e-10);
        end
        
        function test_upsampler_density(obj)
            %TEST_UPSAMPLER_DENSITY Verify density upsampling (interpolation)
            targetGrid = obj.grid.upsample(2);
            upsampler = quadind.grid.Upsampler(obj.grid, targetGrid);
            
            % Create smooth density: q = position (for testing)
            q = obj.grid.x;  % [160 x 3]
            
            % Upsample density
            q_up = upsampler.apply(q);
            obj.assertSize(q_up, [640, 3]);
        end
        
        function test_upsampler_preserves_integral(obj)
            %TEST_UPSAMPLER_PRESERVES_INTEGRAL Upsampling should preserve integral
            targetGrid = obj.grid.upsample(2);
            upsampler = quadind.grid.Upsampler(obj.grid, targetGrid);
            
            % Constant density
            q = ones(160, 1);
            
            % Original integral
            I_orig = sum(obj.grid.w .* q);
            
            % Upsampled integral
            q_up = upsampler.apply(q);
            I_up = sum(targetGrid.w .* q_up);
            
            % Should be approximately equal (relax tolerance for spectral accuracy)
            obj.assertAlmostEqual(I_orig, I_up, 1e-5);
        end
        
        function test_upsampler_high_factor(obj)
            %TEST_UPSAMPLER_HIGH_FACTOR Test higher upsampling factor
            targetGrid = obj.grid.upsample(4);
            upsampler = quadind.grid.Upsampler(obj.grid, targetGrid);
            
            q = ones(160, 1);
            q_up = upsampler.apply(q);
            
            % Should have 16x as many points
            obj.assertSize(q_up, [2560, 1]);
        end
    end
end
