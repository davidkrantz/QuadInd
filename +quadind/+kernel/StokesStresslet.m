classdef StokesStresslet < quadind.kernel.Kernel
    % STOKESSTRESSLET Stokes stresslet (double-layer) kernel
    %
    % The stresslet kernel computes:
    %
    %   u_j(x) = -6 * integral[ q_i(y) * r_i * r_j * r_k * n_k(y) / |r|^5 ] dS(y)
    %
    % where r = x - y is the displacement vector.
    %
    % Properties:
    %   - Singularity order: p = 5/2
    %   - Number of components: 3 (vector density and output)
    %
    % Usage:
    %   kernel = quadind.kernel.StokesStresslet();
    %   u = kernel.evaluate(targets, sources, normals, density, weights);
    %
    % See also: quadind.kernel.Kernel
    
    methods
        function u = evaluate(obj, targets, sources, normals, density, weights)
            % EVALUATE Compute stresslet potential
            %
            %   u = evaluate(obj, targets, sources, normals, density, weights)
            %
            % Inputs:
            %   targets - [M×3] target points
            %   sources - [N×3] source points
            %   normals - [N×3] weighted normals (n * w already multiplied)
            %   density - [N×3] density values
            %   weights - [N×1] additional weights (typically ones if normals are pre-weighted)
            %
            % Output:
            %   u - [M×3] potential values
            %
            % Formula:
            %   u_j = -6 * sum_y [ q_i * r_i * r_j * (r_k * n_k) / |r|^5 ] * w
            
            arguments
                obj
                targets (:,3) {mustBeNumeric}
                sources (:,3) {mustBeNumeric}
                normals (:,3) {mustBeNumeric}
                density (:,3) {mustBeNumeric}
                weights (:,1) {mustBeNumeric}
            end
            
            N = size(sources, 1);
            if size(normals,1) ~= N || size(density,1) ~= N || numel(weights) ~= N
                error('quadind:StokesStresslet:inconsistentSourceArrays', ...
                    'Sources, normals, density, and weights must have the same row count.');
            end
            
            if N == 1
                % Vectorize over targets (single source)
                u = obj.evaluateSingleSource(targets, sources, normals, density, weights);
            else
                % Vectorize over sources (iterate over targets)
                u = obj.evaluateMultipleSources(targets, sources, normals, density, weights);
            end
        end
        
        function p = singularityOrder(~)
            % SINGULARITYORDER Return 5/2 for stresslet
            p = 5/2;
        end
        
        function ncomp = numComponents(~)
            % NUMCOMPONENTS Return 3 for vector kernel
            ncomp = 3;
        end
        
        function name = kernelName(~)
            % KERNELNAME Return descriptive name
            name = 'Stokes Stresslet';
        end
        
        function disp(obj)
            % DISP Display kernel info
            fprintf('StokesStresslet kernel (p = %.1f, ncomp = %d)\n', ...
                obj.singularityOrder(), obj.numComponents());
        end
    end
    
    methods (Access = private)
        function u = evaluateSingleSource(~, targets, source, normal, density, weight)
            % EVALUATESINGLESOURCE Vectorized over targets for single source
            
            % r = target - source
            r = bsxfun(@minus, targets, source);
            
            % |r|^{-5}
            r5inv = sum(r.^2, 2).^(-5/2);
            
            % q_i * r_i (dot product)
            qr = r * density';
            
            % r_k * n_k (dot product)
            rn = r * normal';
            
            % u_j = -6 * q_i * r_i * r_j * r_k * n_k / |r|^5 * w
            u = -6 * bsxfun(@times, r, qr .* rn .* r5inv) * weight;
        end
        
        function u = evaluateMultipleSources(~, targets, sources, normals, density, weights)
            % EVALUATEMULTIPLESOURCES Iterate over targets, vectorize over sources
            
            M = size(targets, 1);
            u = zeros(M, 3);
            
            for i = 1:M
                % r = target(i) - sources
                r = bsxfun(@minus, targets(i, :), sources);
                
                % |r|^{-5}
                r5inv = sum(r.^2, 2).^(-5/2);
                
                % q_i * r_i (row-wise dot product)
                qr = sum(density .* r, 2);
                
                % r_k * n_k (row-wise dot product)
                rn = sum(r .* normals, 2);
                
                % u_j = -6 * sum [ q_i * r_i * r_j * r_k * n_k / |r|^5 * w ]
                u(i, :) = -6 * sum(bsxfun(@times, r, qr .* rn .* r5inv .* weights), 1);
            end
        end
    end
end
