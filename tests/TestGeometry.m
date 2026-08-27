classdef TestGeometry < QuadIndTestCase
%TESTGEOMETRY Unit tests for +quadind/+geometry module
%
%   Tests for AxsymGeometry, Spheroid, Peanut, and CustomAxsym

    properties (Access = private)
        spheroid
        peanut
        capsule
    end
    
    methods (TestMethodSetup)
        function setUp(obj)
            %SETUP Create test fixtures
            obj.spheroid = quadind.geometry.Spheroid('a', 0.05, 'c', 0.1);
            obj.peanut = quadind.geometry.Peanut('amplitude', 2.25);
            obj.capsule = quadind.geometry.Capsule('R', 1, 'L', 6, 'kappa', 4);
        end
    end

    methods (Test)
        function test_spheroid_creation(obj)
            %TEST_SPHEROID_CREATION Verify spheroid construction
            geom = quadind.geometry.Spheroid('a', 1, 'c', 2);
            obj.assertEqual(geom.a, 1);
            obj.assertEqual(geom.c, 2);
        end
        
        function test_spheroid_parameterization(obj)
            %TEST_SPHEROID_PARAMETERIZATION Verify a(theta), c(theta) formulas
            theta = [0; pi/4; pi/2; pi];
            
            % At theta=0: point is at north pole (0, 0, c)
            % At theta=pi/2: point is at equator (a, 0, 0)
            % At theta=pi: point at south pole (0, 0, -c)
            
            a_vals = obj.spheroid.at(theta);
            c_vals = obj.spheroid.ct(theta);
            
            % At poles, a should be 0
            obj.assertAlmostEqual(a_vals(1), 0, 1e-10);
            obj.assertAlmostEqual(a_vals(4), 0, 1e-10);
            
            % At equator, a should be 0.05, c should be 0
            obj.assertAlmostEqual(a_vals(3), 0.05, 1e-10);
            obj.assertAlmostEqual(c_vals(3), 0, 1e-10);
            
            % c at poles
            obj.assertAlmostEqual(c_vals(1), 0.1, 1e-10);
            obj.assertAlmostEqual(c_vals(4), -0.1, 1e-10);
        end
        
        function test_spheroid_evaluate(obj)
            %TEST_SPHEROID_EVALUATE Verify surface point generation
            theta = [pi/2; pi/2; pi/2];  % Same theta for all
            phi = [0; pi/2; pi];
            
            % evaluate only returns points [N x 3]
            x = obj.spheroid.evaluate(theta, phi);
            
            obj.assertSize(x, [3, 3]);  % 3 points in 3D
            
            % At equator, z should be 0
            obj.assertAlmostEqual(x(:, 3), [0; 0; 0], 1e-10);
            
            % Points should be at distance a from z-axis
            r_xy = sqrt(x(:,1).^2 + x(:,2).^2);
            obj.assertAlmostEqual(r_xy, [0.05; 0.05; 0.05], 1e-10);
            
            % Get normals via jacobian
            [~, n] = obj.spheroid.jacobian(theta, phi);
            
            % Normals should be unit vectors
            n_mag = sqrt(sum(n.^2, 2));
            obj.assertAlmostEqual(n_mag, [1; 1; 1], 1e-10);
        end
        
        function test_spheroid_exterior_test(obj)
            %TEST_SPHEROID_EXTERIOR_TEST Verify inside/outside classification
            % isExterior expects [N x 3] format with points as rows
            % Point clearly outside
            obj.assertTrue(obj.spheroid.isExterior([0.2, 0, 0]));
            obj.assertTrue(obj.spheroid.isExterior([0, 0, 0.2]));
            
            % Point clearly inside
            obj.assertFalse(obj.spheroid.isExterior([0.01, 0, 0]));
            obj.assertFalse(obj.spheroid.isExterior([0, 0, 0.01]));
            
            % Multiple points at once
            pts = [0.2, 0, 0; 0.01, 0, 0];  % First outside, second inside
            mask = obj.spheroid.isExterior(pts);
            obj.assertTrue(mask(1));
            obj.assertFalse(mask(2));
        end
        
        function test_peanut_creation(obj)
            %TEST_PEANUT_CREATION Verify peanut geometry construction
            geom = quadind.geometry.Peanut('amplitude', 2.3);
            obj.assertEqual(geom.amplitude, 2.3);
        end
        
        function test_peanut_parameterization(obj)
            %TEST_PEANUT_PARAMETERIZATION Verify peanut shape values
            theta = [0; pi/2; pi];
            
            a_vals = obj.peanut.at(theta);
            c_vals = obj.peanut.ct(theta);
            
            % At poles, a should be 0
            obj.assertAlmostEqual(a_vals(1), 0, 1e-10);
            obj.assertAlmostEqual(a_vals(3), 0, 1e-10);
            
            % Peanut has c(0) = (amplitude+1)*cos(0) = amplitude+1
            % c(pi) = (amplitude+1)*cos(pi) = -(amplitude+1)
            amp = obj.peanut.amplitude;
            obj.assertAlmostEqual(c_vals(1), amp + 1, 1e-10);
            obj.assertAlmostEqual(c_vals(3), -(amp + 1), 1e-10);
        end
        
        function test_peanut_derivatives(obj)
            %TEST_PEANUT_DERIVATIVES Verify derivatives via finite difference
            theta = pi/3;
            h = 1e-7;
            
            % Numerical derivative of a(theta)
            da_num = (obj.peanut.at(theta + h) - obj.peanut.at(theta - h)) / (2*h);
            da_exact = obj.peanut.dadt(theta);
            obj.assertAlmostEqual(da_num, da_exact, 1e-5);
            
            % Numerical derivative of c(theta)
            dc_num = (obj.peanut.ct(theta + h) - obj.peanut.ct(theta - h)) / (2*h);
            dc_exact = obj.peanut.dcdt(theta);
            obj.assertAlmostEqual(dc_num, dc_exact, 1e-5);
        end
        
        function test_custom_geometry(obj)
            %TEST_CUSTOM_GEOMETRY Verify custom geometry with function handles
            % Define a sphere of radius 1
            at = @(th) sin(th);
            ct = @(th) cos(th);
            dadt = @(th) cos(th);
            dcdt = @(th) -sin(th);
            
            geom = quadind.geometry.CustomAxsym(at, ct, dadt, dcdt);
            
            theta = pi/4;
            obj.assertAlmostEqual(geom.at(theta), sin(pi/4), 1e-10);
            obj.assertAlmostEqual(geom.ct(theta), cos(pi/4), 1e-10);
        end
        
        function test_capsule_creation(obj)
            %TEST_CAPSULE_CREATION Verify capsule geometry construction
            geom = quadind.geometry.Capsule('R', 0.5, 'L', 4, 'kappa', 6);
            obj.assertEqual(geom.R, 0.5);
            obj.assertEqual(geom.L, 4);
            obj.assertEqual(geom.kappa, 6);
        end
        
        function test_capsule_parameterization(obj)
            %TEST_CAPSULE_PARAMETERIZATION Verify capsule shape values
            theta = [0; pi/2; pi];
            
            a_vals = obj.capsule.at(theta);
            c_vals = obj.capsule.ct(theta);
            
            % At poles, a should be 0
            obj.assertAlmostEqual(a_vals(1), 0, 1e-10);
            obj.assertAlmostEqual(a_vals(3), 0, 1e-10);
            
            % At equator, a should be R=1, c should be 0
            obj.assertAlmostEqual(a_vals(2), 1, 1e-10);
            obj.assertAlmostEqual(c_vals(2), 0, 1e-10);
            
            % c at poles: +/- L/2 = +/- 3
            obj.assertAlmostEqual(c_vals(1), 3, 1e-10);
            obj.assertAlmostEqual(c_vals(3), -3, 1e-10);
        end
        
        function test_capsule_derivatives(obj)
            %TEST_CAPSULE_DERIVATIVES Verify derivatives via finite difference
            theta = pi/3;
            h = 1e-7;
            
            % Numerical derivative of at(theta)
            da_num = (obj.capsule.at(theta + h) - obj.capsule.at(theta - h)) / (2*h);
            da_exact = obj.capsule.dadt(theta);
            obj.assertAlmostEqual(da_num, da_exact, 1e-5);
            
            % Numerical derivative of ct(theta)
            dc_num = (obj.capsule.ct(theta + h) - obj.capsule.ct(theta - h)) / (2*h);
            dc_exact = obj.capsule.dcdt(theta);
            obj.assertAlmostEqual(dc_num, dc_exact, 1e-5);
        end
        
        function test_capsule_derivatives_at_poles(obj)
            %TEST_CAPSULE_DERIVATIVES_AT_POLES Verify derivatives at endpoints
            % dadt(0) = R*kappa/tanh(kappa), dcdt(0) = 0
            R = obj.capsule.R;
            kap = obj.capsule.kappa;
            
            da0_expected = R * kap / tanh(kap);
            obj.assertAlmostEqual(obj.capsule.dadt(0), da0_expected, 1e-10);
            obj.assertAlmostEqual(obj.capsule.dcdt(0), 0, 1e-10);
            
            % dadt(pi) = -R*kappa/tanh(kappa), dcdt(pi) = 0
            obj.assertAlmostEqual(obj.capsule.dadt(pi), -da0_expected, 1e-10);
            obj.assertAlmostEqual(obj.capsule.dcdt(pi), 0, 1e-10);
        end
        
        function test_capsule_extents(obj)
            %TEST_CAPSULE_EXTENTS Verify maxRadius and maxHeight
            obj.assertEqual(obj.capsule.maxRadius(), 1);
            obj.assertEqual(obj.capsule.maxHeight(), 3);
        end
        
        function test_capsule_exterior_test(obj)
            %TEST_CAPSULE_EXTERIOR_TEST Verify inside/outside classification
            % Point clearly outside (beyond radius)
            obj.assertTrue(obj.capsule.isExterior([2, 0, 0]));
            % Point clearly outside (beyond length)
            obj.assertTrue(obj.capsule.isExterior([0, 0, 4]));
            
            % Point clearly inside
            obj.assertFalse(obj.capsule.isExterior([0.1, 0, 0]));
            obj.assertFalse(obj.capsule.isExterior([0, 0, 0.1]));
        end
        
        function test_capsule_evaluate(obj)
            %TEST_CAPSULE_EVALUATE Verify surface points
            theta = [0; pi/2; pi];
            phi = [0; 0; 0];
            
            pts = obj.capsule.evaluate(theta, phi);
            
            % North pole: (0, 0, L/2)
            obj.assertAlmostEqual(pts(1,:), [0, 0, 3], 1e-10);
            % Equator: (R, 0, 0)
            obj.assertAlmostEqual(pts(2,:), [1, 0, 0], 1e-10);
            % South pole: (0, 0, -L/2)
            obj.assertAlmostEqual(pts(3,:), [0, 0, -3], 1e-10);
        end
        
        function test_geometry_jacobian(obj)
            %TEST_GEOMETRY_JACOBIAN Verify Jacobian calculation
            theta = pi/3;
            phi = 0;  % jacobian requires both theta and phi
            [J, normals] = obj.spheroid.jacobian(theta, phi);
            
            obj.assertTrue(J > 0, 'Jacobian should be positive');
            
            % For spheroid, can verify analytically
            a = obj.spheroid.at(theta);
            dadt = obj.spheroid.dadt(theta);
            dcdt = obj.spheroid.dcdt(theta);
            J_expected = a * sqrt(dadt^2 + dcdt^2);
            
            obj.assertAlmostEqual(J, J_expected, 1e-10);
            
            % Normal should be unit vector
            n_mag = sqrt(sum(normals.^2, 2));
            obj.assertAlmostEqual(n_mag, 1, 1e-10);
        end
    end
end
