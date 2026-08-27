classdef DensityModifier
    % DENSITYMODIFIER Evaluate density at complex roots for error modification
    %
    % The uniform error indicator assumes unit density. This class computes
    % the modification factor based on the actual density evaluated at
    % complex singularity roots.
    %
    % The modified indicator is:
    %   E_modified = E_uniform * max(|q(theta*, phi0)|, |q(theta0, phi*)|)
    %
    % where (theta*, phi*) is the closest grid point and (theta0, phi0) are
    % the complex roots.
    %
    % See also: quadind.IndicatorEvaluator, quadind.indicator.RootFinder
    
    methods (Static)
        function [qphi, qtheta, phi0, intpind] = computeDensityAtRoots(grid, targets, density, theta0, rootCheckTolerance)
            % COMPUTEDENSITYATROOTS Evaluate density at complex roots
            %
            %   [qphi, qtheta, phi0, intpind] = computeDensityAtRoots(grid, targets, density, theta0, rootCheckTolerance)
            %
            % Uses 2-point linear interpolation to evaluate density at complex roots.
            %
            % Inputs:
            %   grid    - AxsymGrid object
            %   targets - [M×3] target points
            %   density - [N×3] density on grid (N = nth*nph)
            %   theta0  - [M×1] (optional) precomputed theta roots
            %   rootCheckTolerance - optional normalized residual tolerance;
            %                        NaN disables the interpolated-root check
            %
            % Outputs:
            %   qphi   - [M×3] density at phi roots
            %   qtheta - [M×3] density at theta roots
            %   phi0   - [M×1] computed phi roots
            %   intpind - [M×2] indices of interpolation neighbors (for potential reuse)
            
            arguments
                grid (1,1) quadind.grid.AxsymGrid
                targets (:,3) {mustBeNumeric}
                density (:,:) {mustBeNumeric}
                theta0 (:,1) = []
                rootCheckTolerance (1,1) double = NaN
            end
            
            ncomp = size(density, 2);
            if size(density, 1) ~= grid.numPoints()
                error('quadind:DensityModifier:invalidDensity', ...
                    'Density must have %d rows (got %d).', ...
                    grid.numPoints(), size(density, 1));
            end

            % Reshape density for easier indexing: [nth × nph × ncomp]
            dens_3d = reshape(density, grid.nth, grid.nph, ncomp);

            % Find closest quadrature node using axisymmetric grid structure.
            [itheta_star,iphi_star] = quadind.indicator.UniformIndicatorBuilder.findNearestNodes(grid, targets);
            theta_star = grid.theta(itheta_star);
            phi_star = grid.phi(iphi_star);
            
            % Compute phi roots analytically
            phi0 = quadind.indicator.RootFinder.computePhiRoot(grid.geometry, targets, theta_star);
            
            % Evaluate density at phi roots using 2-point linear interpolation
            onAxis = hypot(targets(:,1), targets(:,2)) == 0;
            qphi = zeros(size(targets,1), ncomp);
            intpind = NaN(size(targets,1), 2);
            if any(~onAxis)
                [qphi(~onAxis,:), intpind(~onAxis,:)] = ...
                    quadind.indicator.DensityModifier.interpDensityPhi( ...
                        grid, dens_3d, phi0(~onAxis), itheta_star(~onAxis));
            end
            
            % Compute theta roots if not provided
            if isempty(theta0)
                theta0 = grid.geometry.findThetaRoot(targets, phi_star(:));
            elseif isfinite(rootCheckTolerance)
                quadind.util.Diagnostics.checkInterpolatedThetaRoots( ...
                    grid.geometry, theta0, phi_star(:), targets, rootCheckTolerance);
            end
            
            % Evaluate density at theta roots
            qtheta = quadind.indicator.DensityModifier.interpDensityTheta(grid, dens_3d, theta0, iphi_star);
        end
        
        function modifier = computeModifier(qphi, qtheta, targets)
            % COMPUTEMODIFIER Compute the density modification factor
            %
            %   modifier = computeModifier(qphi, qtheta, targets)
            %
            % Rotates global Cartesian density components into the target's
            % meridional tabulation frame, then returns the maximum root
            % magnitude per component.

            arguments
                qphi (:,3) {mustBeNumeric}
                qtheta (:,3) {mustBeNumeric}
                targets (:,3) {mustBeNumeric}
            end

            if (all(targets(:,1) == 0) && all(targets(:,2) == 0)) || all(targets(:,2) == 0)
                % Special case 1: target on z-axis, rotation is undefined.
                % Special case 2: target in xz-plane, rotation is trivial (identity).
                % In both cases, use unrotated global Cartesian components.
                modifier = max(abs(qphi), abs(qtheta));
                return;
            end

            alpha = atan2(targets(:,2), targets(:,1));
            cosalpha = cos(alpha);
            sinalpha = sin(alpha);
            zsign = ones(size(targets, 1), 1);
            zsign(targets(:,3) < 0) = -1;

            qphiRot = [ ...
                cosalpha .* qphi(:,1) + sinalpha .* qphi(:,2), ...
               -sinalpha .* qphi(:,1) + cosalpha .* qphi(:,2), ...
                zsign .* qphi(:,3)];
            qthetaRot = [ ...
                cosalpha .* qtheta(:,1) + sinalpha .* qtheta(:,2), ...
               -sinalpha .* qtheta(:,1) + cosalpha .* qtheta(:,2), ...
                zsign .* qtheta(:,3)];

            modifier = max(abs(qphiRot), abs(qthetaRot));
        end
    end
    
    methods (Static, Access = {?quadind.indicator.DensityModifier, ?quadind.indicator.UniformIndicatorBuilder})
        function [q_interp, intpind] = interpDensityPhi(grid, dens_3d, phi0, itheta)
            % INTERPDENSITYPHI Vectorized 2-point linear interpolation in phi direction
            %
            % Evaluates density at complex phi values using linear interpolation
            % from the two bracketing phi grid points. Exploits uniform phi
            % spacing for O(1) bracket finding and uses dphi as the
            % denominator to correctly handle the periodic boundary.
            
            M = size(phi0, 1);
            ncomp = size(dens_3d, 3);
            nph = grid.nph;
            dphi = 2*pi / nph;
            
            % Find lower bracket index for real(phi0) in the uniform phi grid.
            % For phi grid [0, dphi, 2*dphi, ..., (nph-1)*dphi]:
            iphi_lo = mod(floor(real(phi0) / dphi), nph) + 1;
            iphi_hi = mod(iphi_lo, nph) + 1;  % next index with periodic wrap
            intpind = [iphi_lo, iphi_hi];
            
            % Fractional position within bracket. Using dphi as denominator
            % avoids the boundary bug where raw phi(1)-phi(nph) ≈ 6.18
            % instead of dphi.
            t = (phi0 - grid.phi(iphi_lo).') ./ dphi;
            w_lo = 1 - t;
            w_hi = t;
            
            % Extract density values at bracket points for all targets and components
            q_interp = zeros(M, ncomp);
            for c = 1:ncomp
                lin_lo = sub2ind([grid.nth, nph], itheta, iphi_lo);
                lin_hi = sub2ind([grid.nth, nph], itheta, iphi_hi);
                
                dens_c = dens_3d(:,:,c);  % [nth × nph]
                q_interp(:, c) = w_lo .* dens_c(lin_lo) + w_hi .* dens_c(lin_hi);
            end
        end
        
        function q_interp = interpDensityTheta(grid, dens_3d, theta0, iphi)
            % INTERPDENSITYTHETA Vectorized 2-point linear interpolation in theta direction
            %
            % Finds the true bracketing interval for real(theta0) in the GL
            % grid and interpolates at the complex theta0 value.
            
            M = size(theta0, 1);
            ncomp = size(dens_3d, 3);
            nth = grid.nth;
            theta_grid = grid.theta(:);  % descending: theta(1) > ... > theta(nth)
            
            % Find bracket: index k such that theta(k) >= real(theta0) >= theta(k+1).
            % Since GL nodes descend, count how many exceed real(theta0).
            k = sum(theta_grid.' > real(theta0), 2);  % [M×1]
            k = max(k, 1);
            k = min(k, nth - 1);
            idx_hi = k;        % theta_grid(idx_hi) >= real(theta0)
            idx_lo = k + 1;    % theta_grid(idx_lo) <= real(theta0)
            
            % Theta values at bracket endpoints
            t_hi = theta_grid(idx_hi);  % larger theta value
            t_lo = theta_grid(idx_lo);  % smaller theta value
            
            % Interpolation weights
            denom = t_hi - t_lo;  % always positive for valid brackets
            small_denom = abs(denom) < eps;
            denom(small_denom) = 1;  % avoid division by zero
            
            w_hi = (theta0 - t_lo) ./ denom;
            w_lo = (t_hi - theta0) ./ denom;
            
            % For degenerate intervals, use the hi-index value
            w_hi(small_denom) = 1;
            w_lo(small_denom) = 0;
            
            % Extract density values at bracket points for all targets
            q_interp = zeros(M, ncomp);
            for c = 1:ncomp
                dens_c = dens_3d(:,:,c);  % [nth × nph]
                lin_hi = sub2ind([nth, grid.nph], idx_hi, iphi);
                lin_lo = sub2ind([nth, grid.nph], idx_lo, iphi);
                
                q_interp(:, c) = w_hi .* dens_c(lin_hi) + w_lo .* dens_c(lin_lo);
            end
            
            % Handle NaN (can occur for theta0 at poles)
            q_interp(isnan(q_interp)) = 0;
        end
    end
end
