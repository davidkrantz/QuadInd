# QuadEst: Quadrature Error Estimation for Axisymmetric Geometries

A MATLAB package for computing error estimates for layer potentials on axisymmetric 3D geometries.

**Authors**: Pritpal 'Pip' Matharu, David Krantz (MPI MiS / KTH)

## Installation

```matlab
% Option A: From the QuadEst directory in MATLAB
addpath('.')          % or: addpath(pwd)
init                  % also adds legacy/ to path

% Option B: From any location
addpath('/path/to/QuadEst');
```

## Quick Start

```matlab
% 1. Define geometry
geom = quadest.geometry.Spheroid('a', 0.05, 'c', 0.1);

% 2. Select kernel
kernel = quadest.kernel.StokesStresslet();

% 3. Precompute error estimates
estimator = quadest.errorest.ErrorEstimator(geom, kernel, ...
    'nth', 40, 'nph', 60, 'upsampFactors', 1:6);

% 4. Define density and targets
density = ones(40*60, 3);  % [N×3] on grid
targets = [0.1, 0, 0];     % [M×3] evaluation points

% 5. Evaluate
estimates = estimator.evaluate(targets, density);

% With classification for tolerance-based method selection:
[estimates, classification] = estimator.evaluate(targets, density, 'tol', 1e-6);
```

## Package Structure

```
+quadest/
├── +geometry/          # Axisymmetric geometry definitions
│   ├── AxsymGeometry.m # Abstract base class
│   ├── Spheroid.m      # Prolate/oblate spheroid
│   ├── Capsule.m       # Capsule shape
│   ├── Peanut.m        # Peanut shape
│   └── CustomAxsym.m   # User-defined parameterization (function handles)
├── +kernel/            # Layer potential kernels
│   ├── Kernel.m        # Abstract kernel interface
│   └── StokesStresslet.m  # Stokes stresslet (p = 5/2)
├── +grid/              # Quadrature discretization
│   ├── AxsymGrid.m     # GL×Trapezoidal grid
│   └── Upsampler.m     # Trig + barycentric Lagrange interpolation
├── +errorest/          # Error estimation (core module)
│   ├── ErrorEstimator.m      # Main API
│   ├── UniformEstimateBuilder.m
│   ├── RootFinder.m
│   └── DensityModifier.m
├── +util/              # Utilities
│   ├── Config.m        # Default parameters
│   ├── Diagnostics.m   # Warnings, logging
│   ├── GaussLegendre.m # GL nodes/weights
│   ├── GaussLaguerre8.m
│   └── Plotting.m      # Static visualization utilities
└── +test/              # Unit tests
    ├── runAllTests.m   # Test runner
    └── Test*.m         # Test classes (Util, Geometry, Grid, Kernel,
                        #   ErrorEstimation, LegacyComparison)

examples/               # Demo scripts
legacy/                 # Original research code (for regression testing)
init.m                  # Path initialization
```

## Examples

Run the example scripts in `examples/`:

```matlab
% Minimal example
run('examples/demo_minimal.m');

% Near-field example
run('examples/demo_nearfield.m');

% Numerical validation against measured quadrature errors
run('examples/demo_validation.m');
```

## Running Tests

```matlab
% Run all unit tests
results = quadest.test.runAllTests();
disp(results.summary);

% With verbose output or stop on first failure
results = quadest.test.runAllTests('verbose', true, 'stopOnFailure', true);

% Run a single test class
results = quadest.test.TestGeometry().run();

% Include legacy regression tests (requires legacy/ on path)
addpath(genpath('legacy'));
results = quadest.test.TestLegacyComparison().run();
```

## Requirements

- Statistics and Machine Learning Toolbox (`knnsearch`) — used only in `TestLegacyComparison.m`; the core error estimation and density interpolation paths use a parametric nearest-node approach and do **not** require this toolbox

## Adding New Geometries

Subclass `quadest.geometry.AxsymGeometry` and implement:

```matlab
classdef MyGeometry < quadest.geometry.AxsymGeometry
    methods
        function val = at(obj, theta)
            % Radial distance in xy-plane
        end
        function val = ct(obj, theta)
            % Z-coordinate
        end
        function val = dadt(obj, theta)
            % d(at)/d(theta)
        end
        function val = dcdt(obj, theta)
            % d(ct)/d(theta)
        end
        function val = maxRadius(obj), end
        function val = maxHeight(obj), end
        function mask = isExterior(obj, points), end
        function theta = findThetaRoot(obj, targets, phi), end
    end
end
```

## Adding New Kernels

Subclass `quadest.kernel.Kernel` and implement:

```matlab
classdef MyKernel < quadest.kernel.Kernel
    methods
        function u = evaluate(obj, targets, sources, normals, density, weights)
            % Compute kernel
        end
        function p = singularityOrder(obj)
            % Return singularity order (e.g., 5/2 for stresslet)
        end
        function n = numComponents(obj)
            % Return number of output components
        end
    end
end
```

## References
If you find **QuadEst** useful in your work, we ask that you please cite the following works:
* This software. TODO.
* Our paper. TODO.
* L. af Klinteberg, C. Sorgentone, and A.-K. Tornberg, *Quadrature error estimates for layer potentials evaluated near curved surfaces in three dimensions*, Computers & Mathematics with Applications, 111 (2022), pp. 1–19, https://doi.org/10.1016/j.camwa.2022.02.001.
* Sorgentone and A.-K. Tornberg, *Estimation of quadrature errors for layer potentials evaluated near surfaces with spherical topology*, Advances in Computational Mathematics, 49 (2023), p. 87, https://doi.org/10.1007/s10444-023-10083-7.
