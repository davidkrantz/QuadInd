classdef (Abstract) Kernel
    % KERNEL Abstract base class for integral equation kernels
    %
    % Defines the interface that all kernels must implement for use with
    % the error-indicator evaluation framework.
    %
    % To add a new kernel:
    %   1. Create a subclass in +quadind/+kernel/
    %   2. Implement evaluate(), singularityOrder(), and numComponents()
    %
    % See also: quadind.kernel.StokesStresslet
    
    methods (Abstract)
        % EVALUATE Compute kernel potential at target points
        %
        %   u = evaluate(obj, targets, sources, normals, density, weights)
        %
        % Inputs:
        %   targets - [M×3] target points
        %   sources - [N×3] source points (on surface)
        %   normals - [N×3] unit normal vectors at sources
        %   density - [N×ncomp] density values at sources
        %   weights - [N×1] quadrature weights
        %
        % Output:
        %   u - [M×ncomp] potential values at targets
        u = evaluate(obj, targets, sources, normals, density, weights)
        
        % SINGULARITYORDER Return the singularity order p
        %
        p = singularityOrder(obj)
        
        % NUMCOMPONENTS Number of density/output components
        %
        % For scalar kernels: ncomp = 1
        % For vector kernels (Stokes): ncomp = 3
        ncomp = numComponents(obj)
    end
    
    methods
        function u = evaluateOnGrid(obj, targets, grid, density)
            % EVALUATEONGRID Evaluate kernel using an AxsymGrid
            %
            %   u = evaluateOnGrid(obj, targets, grid, density)
            %
            % Convenience method that extracts sources, normals, weights
            % from the grid object.
            %
            % Inputs:
            %   targets - [M×3] target points
            %   grid    - AxsymGrid object
            %   density - [N×ncomp] density on grid
            %
            % Output:
            %   u - [M×ncomp] potential values
            
            arguments
                obj
                targets (:,3) {mustBeNumeric}
                grid (1,1) quadind.grid.AxsymGrid
                density (:,:) {mustBeNumeric}
            end
            
            if size(density, 1) ~= grid.numPoints() || size(density, 2) ~= obj.numComponents()
                error('quadind:Kernel:invalidDensity', ...
                    'Density must have size %d-by-%d.', grid.numPoints(), obj.numComponents());
            end
            % Combine normals and weights for efficiency
            normals_weighted = bsxfun(@times, grid.n, grid.w);
            
            u = obj.evaluate(targets, grid.x, normals_weighted, density, ones(size(grid.w)));
        end
        
        function u = evaluateUpsampled(obj, targets, grid, density, factor)
            % EVALUATEUPSAMPLED Evaluate with upsampled grid
            %
            %   u = evaluateUpsampled(obj, targets, grid, density, factor)
            %
            % Upsamples the grid and density by the given factor before
            % evaluating the kernel.
            %
            % Inputs:
            %   targets - [M×3] target points
            %   grid    - AxsymGrid object (coarse)
            %   density - [N×ncomp] density on coarse grid
            %   factor  - Upsampling factor (integer >= 1)
            %
            % Output:
            %   u - [M×ncomp] potential values
            
            arguments
                obj
                targets (:,3) {mustBeNumeric}
                grid (1,1) quadind.grid.AxsymGrid
                density (:,:) {mustBeNumeric}
                factor (1,1) {mustBePositive, mustBeInteger} = 1
            end
            
            if factor == 1
                u = obj.evaluateOnGrid(targets, grid, density);
                return;
            end
            
            % Create upsampled grid
            fineGrid = grid.upsample(factor);
            
            % Interpolate density
            interpolator = grid.getInterpolator(fineGrid);
            density_fine = interpolator.apply(density);
            
            % Evaluate on fine grid
            u = obj.evaluateOnGrid(targets, fineGrid, density_fine);
        end
        
        function name = kernelName(obj)
            % KERNELNAME Return human-readable kernel name
            %
            % Default implementation returns class name. Override for
            % custom names.
            
            name = class(obj);
            % Remove package prefix
            parts = strsplit(name, '.');
            name = parts{end};
        end
    end
end
