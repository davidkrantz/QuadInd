classdef DensityModifier
    % DENSITYMODIFIER Evaluate density at complex roots for error modification
    %
    % The uniform error estimate assumes unit density. This class computes
    % the modification factor based on the actual density evaluated at
    % complex singularity roots.
    %
    % The modified estimate is:
    %   E_modified = E_uniform * max(|q(theta*, phi0)|, |q(theta0, phi*)|)
    %
    % where (theta*, phi*) is the closest grid point and (theta0, phi0) are
    % the complex roots.
    %
    % See also: quadest.errorest.ErrorEstimator, quadest.errorest.RootFinder
    
    methods (Static)
        function [qphi, qtheta, phi0, intpind] = computeDensityAtRoots(grid, targets, density, theta0)
            % COMPUTEDENSITYATROOTS Evaluate density magnitude at complex roots
            %
            %   [qphi, qtheta, phi0, intpind] = computeDensityAtRoots(grid, targets, density, theta0)
            %
            % Uses 2-point linear interpolation to evaluate density at complex roots.
            %
            % Inputs:
            %   grid    - AxsymGrid object
            %   targets - [M×3] target points
            %   density - [N×3] density on grid (N = nth*nph)
            %   theta0  - [M×1] (optional) precomputed theta roots
            %
            % Outputs:
            %   qphi   - [M×3] density magnitude at phi roots
            %   qtheta - [M×3] density magnitude at theta roots
            %   phi0   - [M×1] computed phi roots
            %   intpind - [M×2] indices of interpolation neighbors (for potential reuse)
            
            arguments
                grid (1,1) quadest.grid.AxsymGrid
                targets (:,3) {mustBeNumeric}
                density (:,:) {mustBeNumeric}
                theta0 (:,1) = []
            end
            
            M = size(targets, 1);
            ncomp = size(density, 2);
            
            % Reshape density for easier indexing: [nth × nph × ncomp]
            dens_3d = reshape(density, grid.nth, grid.nph, ncomp);
            
            % Find closest quadrature node for each target
            idx_closest = knnsearch(grid.x, targets);
            [itheta_star, iphi_star] = grid.ind2sub(idx_closest);
            
            theta_star = grid.theta(itheta_star);
            phi_star = grid.phi(iphi_star);
            
            % Compute phi roots analytically
            phi0 = quadest.errorest.RootFinder.computePhiRoot(grid.geometry, targets, theta_star);
            
            % Evaluate density at phi roots using 2-point linear interpolation
            [qphi, intpind] = quadest.errorest.DensityModifier.interpDensityPhi(grid, dens_3d, phi0, itheta_star);
            
            % Compute theta roots if not provided
            if isempty(theta0)
                theta0 = grid.geometry.findThetaRoot(targets, phi_star(:));
            end
            
            % Evaluate density at theta roots
            qtheta = quadest.errorest.DensityModifier.interpDensityTheta(grid, dens_3d, theta0, theta_star, iphi_star);
        end
        
        function modifier = computeModifier(qphi, qtheta)
            % COMPUTEMODIFIER Compute the density modification factor
            %
            %   modifier = computeModifier(qphi, qtheta)
            %
            % Returns max(|qphi|, |qtheta|) per component.
            
            modifier = max(abs(qphi), abs(qtheta));
        end
    end
    
    methods (Static, Access = private)
        function [q_interp, intpind] = interpDensityPhi(grid, dens_3d, phi0, itheta)
            % INTERPDENSITYPHI 2-point linear interpolation in phi direction
            %
            % Evaluates density at complex phi values using linear interpolation
            % from the two nearest phi grid points.
            
            M = size(phi0, 1);
            ncomp = size(dens_3d, 3);
            n_intp = 2;  % Number of interpolation points
            
            q_interp = zeros(M, ncomp);
            intpind = zeros(M, n_intp);
            
            phi_grid = grid.phi(:);
            nph = grid.nph;
            
            for ii = 1:M
                phi_val = real(phi0(ii));
                theta_idx = itheta(ii);
                
                % Find two nearest phi points (accounting for periodicity)
                tmpphi = phi_grid;
                [~, sorted_idx] = sort(abs(tmpphi - phi_val));
                
                % Handle periodicity at boundaries
                if sorted_idx(1) < n_intp
                    % Near phi = 0, wrap around
                    tmpphi_wrapped = tmpphi;
                    tmpphi_wrapped(end-n_intp+1:end) = tmpphi(end-n_intp+1:end) - 2*pi;
                    [~, sorted_idx] = sort(abs(tmpphi_wrapped - phi_val));
                elseif sorted_idx(1) > nph - n_intp
                    % Near phi = 2*pi, wrap around
                    tmpphi_wrapped = tmpphi;
                    tmpphi_wrapped(1:n_intp) = tmpphi(1:n_intp) + 2*pi;
                    [~, sorted_idx] = sort(abs(tmpphi_wrapped - phi_val));
                end
                
                local_ip = sort(sorted_idx(1:n_intp));
                intpind(ii, :) = local_ip;
                
                % Get phi values for interpolation
                p_pts = phi_grid(local_ip);
                
                % Handle wrapped values for interpolation formula
                if abs(p_pts(2) - p_pts(1)) > pi
                    % Wrapped case
                    if p_pts(1) < pi
                        p_pts(2) = p_pts(2) - 2*pi;
                    else
                        p_pts(1) = p_pts(1) - 2*pi;
                    end
                end
                
                % 2-point linear interpolation for each component
                for c = 1:ncomp
                    q_vals = squeeze(dens_3d(theta_idx, local_ip, c));
                    % Linear interpolation at complex phi0
                    q_interp(ii, c) = (q_vals(1) * (p_pts(2) - phi0(ii)) + ...
                                       q_vals(2) * (phi0(ii) - p_pts(1))) / ...
                                      (p_pts(2) - p_pts(1));
                end
            end
        end
        
        function q_interp = interpDensityTheta(grid, dens_3d, theta0, theta_star, iphi)
            % INTERPDENSITYTHETA 2-point linear interpolation in theta direction
            
            M = size(theta0, 1);
            ncomp = size(dens_3d, 3);
            n_intp = 2;
            
            q_interp = zeros(M, ncomp);
            theta_grid = grid.theta(:);
            
            for ii = 1:M
                theta_val = theta_star(ii);
                phi_idx = iphi(ii);
                
                % Find two nearest theta points
                [~, sorted_idx] = mink(abs(theta_grid - theta_val), n_intp);
                local_it = sort(sorted_idx);
                
                t_pts = theta_grid(local_it);
                
                % 2-point linear interpolation for each component
                for c = 1:ncomp
                    q_vals = squeeze(dens_3d(local_it, phi_idx, c));
                    
                    % Linear interpolation at complex theta0
                    denom = t_pts(2) - t_pts(1);
                    if abs(denom) < eps
                        q_interp(ii, c) = q_vals(1);
                    else
                        q_interp(ii, c) = (q_vals(1) * (t_pts(2) - theta0(ii)) + ...
                                           q_vals(2) * (theta0(ii) - t_pts(1))) / denom;
                    end
                end
                
                % Handle NaN (can occur for theta0 at poles)
                q_interp(isnan(q_interp)) = 0;
            end
        end
    end
end
