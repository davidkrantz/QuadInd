classdef TestGeometry < quadest.test.TestBase
%TESTGEOMETRY Unit tests for +quadest/+geometry module
%
%   Tests for AxsymGeometry, Spheroid, Peanut, and CustomAxsym

    properties (Access = private)
        spheroid
        peanut
    end
    
    methods
        function setUp(obj)
            %SETUP Create test fixtures
            obj.spheroid = quadest.geometry.Spheroid('a', 0.05, 'c', 0.1);
            obj.peanut = quadest.geometry.Peanut('amplitude', 2.25);
        end
        
        function test_spheroid_creation(obj)
            %TEST_SPHEROID_CREATION Verify spheroid construction
            geom = quadest.geometry.Spheroid('a', 1, 'c', 2);
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
            geom = quadest.geometry.Peanut('amplitude', 2.3);
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
            
            geom = quadest.geometry.CustomAxsym(at, ct, dadt, dcdt);
            
            theta = pi/4;
            obj.assertAlmostEqual(geom.at(theta), sin(pi/4), 1e-10);
            obj.assertAlmostEqual(geom.ct(theta), cos(pi/4), 1e-10);
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
