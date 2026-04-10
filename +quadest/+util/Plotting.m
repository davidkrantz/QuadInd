classdef Plotting
    % PLOTTING Static utility class for visualization
    %
    % Provides plotting functions for error estimation validation:
    %   - Quadrature region maps
    %   - Error vs distance scatter plots
    %   - Error contours with estimate overlays
    %   - Diverging colormaps
    %
    % Usage:
    %   quadest.util.Plotting.plotQuadRegions(...)
    %   quadest.util.Plotting.plotErrorContour(...)
    %
    % See also: quadest.errorest.ErrorEstimator
    
    methods (Static)
        function plotQuadRegions(xv, zv, classification, mask_ext, geom, varargin)
            % PLOTQUADREGIONS Plot quadrature regions colored by upsampling factor
            %
            %   plotQuadRegions(xv, zv, classification, mask_ext, geom)
            %   plotQuadRegions(..., 'FontSize', 16)
            %
            % Inputs:
            %   xv, zv       - Coordinate vectors for meshgrid
            %   classification - Struct from ErrorEstimator.evaluate with .masks, .isSQ
            %   mask_ext     - Logical mask for exterior points (size [numel(xv)*numel(zv), 1])
            %   geom         - AxsymGeometry object for surface outline
            %
            % Options:
            %   'FontSize'   - Font size for labels (default: 16)
            
            p = inputParser;
            addParameter(p, 'FontSize', 16);
            parse(p, varargin{:});
            FS = p.Results.FontSize;
            
            nfac = length(classification.masks);
            M = numel(xv);
            N = numel(zv);
            
            % Build region matrix
            rquad = NaN(M, N);
            
            % Direct quadrature (factor = 1) gets highest value
            directMask = false(numel(mask_ext), 1);
            directMask(mask_ext) = classification.masks{1};
            rquad(reshape(directMask, M, N)) = nfac;
            
            % Upsampling regions get decreasing values
            for k = 2:nfac
                upMask = false(numel(mask_ext), 1);
                upMask(mask_ext) = classification.masks{k};
                rquad(reshape(upMask, M, N)) = nfac - k + 1;
            end
            
            % Special quadrature gets 0
            sqMask = false(numel(mask_ext), 1);
            sqMask(mask_ext) = classification.isSQ;
            rquad(reshape(sqMask, M, N)) = 0;
            
            % Negate for colormap
            rquad = -abs(rquad) - 0.5;
            
            % Create contour levels
            levels = -0.5:-1:(-nfac - 0.5);
            
            % Plot
            [X, Z] = meshgrid(xv, zv);
            contourf(X, Z, rquad, levels, 'LineColor', 'none');
            hold on;
            
            % Draw geometry outline
            quadest.util.Plotting.drawGeometryOutline(geom);
            
            view(0, 90);
            axis equal;
            clim([-length(levels), 0]);
            xlim([min(xv), max(xv)]);
            ylim([min(zv), max(zv)]);
            
            % Colorbar with legend
            cbar = colorbar('Ticks', fliplr(levels), 'TickLabelInterpreter', 'latex', 'FontSize', FS);
            
            % Build legend entries
            legents = cell(1, nfac + 1);
            legents{1} = 'Direct Quadrature';
            for k = 2:nfac
                legents{k} = sprintf('Upsampling = %d', k);
            end
            legents{nfac + 1} = 'Special Quadrature';
            cbar.TickLabels = legents;
            
            % Apply colormap
            colormap(gca, quadest.util.Plotting.divergingColormap(nfac + 1));
            
            xlabel('$x$', 'Interpreter', 'latex', 'FontSize', FS);
            ylabel('$z$', 'Interpreter', 'latex', 'FontSize', FS);
            set(gca, 'FontSize', FS);
            hold off;
        end
        
        function plotErrorVsDistance(u, uref, targets, gridSurf, classification, tol, varargin)
            % PLOTERRORVSDISTANCE Scatter plot of error vs distance from surface
            %
            %   plotErrorVsDistance(u, uref, targets, gridSurf, classification, tol)
            %
            % Inputs:
            %   u, uref       - Computed and reference solutions [M×3]
            %   targets       - Target points [M×3]
            %   gridSurf      - Grid object containing surface points
            %   classification - Struct with .masks, .isSQ
            %   tol           - Tolerance for horizontal line
            %
            % Options:
            %   'FontSize'    - Font size for labels (default: 16)
            
            p = inputParser;
            addParameter(p, 'FontSize', 16);
            parse(p, varargin{:});
            FS = p.Results.FontSize;
            
            % Compute error at each target
            err = sqrt(sum((u - uref).^2, 2));
            
            % Compute distance to nearest surface point
            [itheta,iphi] = quadest.errorest.UniformEstimateBuilder.findNearestNodes(gridSurf, targets);
            % Linear index into grid.x:
            idx = (iphi - 1) * gridSurf.nth + itheta;
            dist = sqrt(sum((targets - gridSurf.x(idx, :)).^2, 2));
            
            nfac = length(classification.masks);
            
            % Build legend
            legents = cell(1, nfac + 1);
            legents{1} = 'Direct Quadrature';
            for k = 2:nfac
                legents{k} = sprintf('Upsampling = %d', k);
            end
            legents{nfac + 1} = 'Special Quadrature';
            
            % Plot each region
            ax = gca;
            hold on;
            
            % Direct quadrature (black)
            if any(classification.masks{1})
                plot(dist(classification.masks{1}), err(classification.masks{1}), 'k.', 'MarkerSize', 8);
            end
            
            % Upsampling regions (colors from ColorOrder)
            ax.ColorOrderIndex = 2;
            for k = 2:nfac
                if any(classification.masks{k})
                    plot(dist(classification.masks{k}), err(classification.masks{k}), '.', 'MarkerSize', 8);
                end
            end
            
            % Special quadrature (first color)
            ax.ColorOrderIndex = 1;
            if any(classification.isSQ)
                plot(dist(classification.isSQ), err(classification.isSQ), '.', 'MarkerSize', 8);
            end
            
            % Tolerance line
            yline(tol, 'k--', 'Tolerance', 'LineWidth', 1.5, 'FontSize', FS);
            
            set(gca, 'XScale', 'log', 'YScale', 'log');
            xlabel('Distance from surface', 'FontSize', FS);
            ylabel('Error', 'FontSize', FS);
            title('Solution Error', 'FontSize', FS);
            grid on;
            legend(legents, 'Interpreter', 'latex', 'Location', 'best');
            set(gca, 'FontSize', FS);
            hold off;
        end
        
        function plotErrorContour(xv, zv, error, estimates, mask_ext, geom, varargin)
            % PLOTERRORCONTOUR Filled contour of error with estimate overlay
            %
            %   Single upsampling factor:
            %   plotErrorContour(xv, zv, error, estimates, mask_ext, geom)
            %   plotErrorContour(..., 'uniformEstimates', estuni)
            %   plotErrorContour(..., 'levels', [0:-2:-10])
            %
            %   Multiple upsampling factors (one figure per factor):
            %   plotErrorContour(xv, zv, errors_cell, estimates_cell, mask_ext, geom, ...
            %       'upsampFactors', [1, 2, 3, 4])
            %
            % Inputs:
            %   xv, zv     - Coordinate vectors
            %   error      - [M×1] error values at exterior points, OR
            %                cell array {k} of [M×1] vectors (one per upsampling factor)
            %   estimates  - [M×1] error estimates at exterior points, OR
            %                cell array {k} of [M×1] vectors (one per upsampling factor)
            %   mask_ext   - Logical mask for exterior points
            %   geom       - AxsymGeometry object
            %
            % Options:
            %   'upsampFactors'    - Upsampling factor values for subplot titles;
            %                        required when error/estimates are cell arrays
            %   'uniformEstimates' - [M×1] uniform (unmodified) estimates for overlay
            %                        (single factor only)
            %   'levels'           - Contour levels (default: [0:-2:-10])
            %   'FontSize'         - Font size (default: 16)
            
            p = inputParser;
            addParameter(p, 'uniformEstimates', []);
            addParameter(p, 'levels', 0:-2:-10);
            addParameter(p, 'FontSize', 16);
            addParameter(p, 'upsampFactors', []);
            parse(p, varargin{:});
            
            levels        = p.Results.levels;
            FS            = p.Results.FontSize;
            upsampFactors = p.Results.upsampFactors;
            
            if iscell(error)
                % Multi-factor mode: one subplot per upsampling factor
                nfac  = length(error);
                
                for k = 1:nfac
                    figure('Name', sprintf('Error Contour - Factor %d', k));
                    
                    % Build subplot title from factor value
                    if ~isempty(upsampFactors) && k <= length(upsampFactors)
                        fac = upsampFactors(k);
                        if fac == 1
                            titleStr = 'Direct (factor 1)';
                        else
                            titleStr = sprintf('Upsampling \\times%d', fac);
                        end
                    else
                        titleStr = sprintf('Factor %d', k);
                    end
                    
                    % Per-factor estimates (cell or shared single array)
                    if iscell(estimates)
                        est_k = estimates{k};
                    else
                        est_k = estimates;
                    end
                    
                    quadest.util.Plotting.plotSingleFactor_(xv, zv, error{k}, est_k, ...
                        mask_ext, geom, levels, FS);
                    title(titleStr, 'FontSize', FS);
                end
                return;
            end
            
            % Single-factor mode (original behaviour)
            quadest.util.Plotting.plotSingleFactor_(xv, zv, error, estimates, ...
                mask_ext, geom, levels, FS, p.Results.uniformEstimates);
        end
        
        function plotSurfaceMesh(gridOrGeom, varargin)
            % PLOTSURFACEMESH 3-D surface mesh of the particle surface
            %
            %   plotSurfaceMesh(grid)
            %   plotSurfaceMesh(geom)
            %   plotSurfaceMesh(geom, 'nth', 40, 'nph', 60)
            %   plotSurfaceMesh(..., 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k')
            %
            % Inputs:
            %   gridOrGeom - AxsymGrid object, OR any AxsymGeometry object.
            %                When a geometry is supplied a default GL×trap grid
            %                with nth=40, nph=60 (or values passed via 'nth'/'nph')
            %                is created internally.
            %
            % Options:
            %   'nth'       - Theta points when gridOrGeom is a geometry (default: 40)
            %   'nph'       - Phi points when gridOrGeom is a geometry (default: 60)
            %   'FaceColor' - surf FaceColor (default: 0.8*[1 1 1])
            %   'EdgeColor' - surf EdgeColor (default: [0.3 0.3 0.3])
            %   'FaceAlpha' - surf FaceAlpha (default: 1)

            p = inputParser;
            addParameter(p, 'nth',       40);
            addParameter(p, 'nph',       60);
            addParameter(p, 'FaceColor', 0.8*[1 1 1]);
            addParameter(p, 'EdgeColor', [0.3 0.3 0.3]);
            addParameter(p, 'FaceAlpha', 1);
            parse(p, varargin{:});

            % Resolve grid
            if isa(gridOrGeom, 'quadest.grid.AxsymGrid')
                g = gridOrGeom;
            else
                % Geometry supplied — build a temporary grid
                g = quadest.grid.AxsymGrid(gridOrGeom, ...
                    'nth', p.Results.nth, 'nph', p.Results.nph);
            end

            nth = g.nth;
            nph = g.nph;

            % Reshape into (nth × nph) matrices
            X = reshape(g.x(:,1), nth, nph);
            Y = reshape(g.x(:,2), nth, nph);
            Z = reshape(g.x(:,3), nth, nph);

            % Wrap azimuthally so the surface closes
            X = [X, X(:,1)];
            Y = [Y, Y(:,1)];
            Z = [Z, Z(:,1)];

            surf(X, Y, Z, ...
                'FaceColor', p.Results.FaceColor, ...
                'EdgeColor', p.Results.EdgeColor, ...
                'FaceAlpha', p.Results.FaceAlpha);
            axis equal;
        end

        function drawGeometryOutline(geom, varargin)
            % DRAWGEOMETRYOUTLINE Draw the geometry cross-section at y=0
            %
            %   drawGeometryOutline(geom)
            %   drawGeometryOutline(geom, 'Color', 'w', 'LineWidth', 2)
            
            p = inputParser;
            addParameter(p, 'Color', [0.8, 0.8, 0.8]);
            addParameter(p, 'LineWidth', 2);
            parse(p, varargin{:});
            
            theta_plot = linspace(0, pi, 100)';
            x_sph = geom.at(theta_plot);
            z_sph = geom.ct(theta_plot);
            
            % Draw both sides
            plot(x_sph, z_sph, '-', 'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth);
            plot(-x_sph, z_sph, '-', 'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth);
        end
        
        function map = divergingColormap(m)
            % DIVERGINGCOLORMAP Red-blue diverging colormap
            %
            %   map = divergingColormap(m) returns an m×3 colormap
            %
            % Based on Kenneth Moreland's diverging color map algorithm.
            % Uses perceptually uniform interpolation through CIELAB space.
            
            if nargin < 1
                m = size(get(gcf, 'colormap'), 1);
            end
            
            % RGB endpoints (blue to red)
            rgb1 = [0.230, 0.299, 0.754];  % Blue
            rgb2 = [0.706, 0.016, 0.150];  % Red
            
            s = linspace(0, 1, m);
            map = zeros(m, 3);
            for i = 1:m
                map(i, :) = quadest.util.Plotting.divergingMapValue(s(i), rgb1, rgb2);
            end
        end

        function alignfigs(varargin)
        %alignfigs   Arrange all open figure windows nicely.
        %   ALIGNFIGS() arranges all the open current figure windows to fit nicely on the
        %   screen in a grid. If they won't fit, then the size is reduced uniformly.
        %
        %   ALIGNFIGS(N) specifies that there should be N columns.
        %
        %   Note that if not all windows have the same size, they will be adjusted to be
        %   so (to the smallest open window).
        %
        %   Examples:
        %       close all, for k = 1:8, figure, end, alignfigs()
        %       close all, for k = 1:15, figure, end, alignfigs(5)

        % Nick Hale - Dec 2014.

        h = findobj('Type', 'figure');
        nh = numel(h);

        % Windows are taller than they claim to be?
        fudge1 = 182; % Also the top on is different...
        fudgeN = 158;

        % Use these values if removing the toolbar.
        % fudge1 = 130;
        % fudgeN = 100;

        if ( nh < 2 )
            return
        end

        % Sort into the right order:
        [~, idx] = sort(cell2mat(get(h, 'Number')));
        h = h(idx);

        % Get the sizes:
        for k = 1:numel(h)
            p(k,:) = get(h(k), 'position');
        %     set(h(k), 'toolbar', 'none')
        end
        % Adjust to be the same size:
        p(:,3:4) = repmat(min(p(:,3:4)),nh, 1);

        % Get the screen size:
        ss = get(0, 'ScreenSize');
        ss = ss(3:4);

        % Determine the grid size:
        nxMax = floor(ss(1)/p(1,3));
        if ( nargin == 0 )
            nx = nxMax;
        else
            nx = varargin{1};
        end

        % Screen isn't large enough! Shrink the figures:
        ny = floor( (ss(2)+fudgeN) / (p(1,4)+fudgeN) );
        if ( nh > nx*ny || nx > nxMax )
            for k = 1:nh
                set(h(k), 'position', [p(k,1:2), .9*p(k,3:4)]);
            end
            alignfigs(varargin{:});
            return
        end

        % Add a margin at the left to centre:
        leftMargin = (ss(1) - nx*p(1,3))/2;

        % Place the first figure in the top lefT:
        p(1,:) = [leftMargin ss(2) - p(1,4) p(1,3:4)];
        set(h(1), 'position', p(1,:));
        p(1,:) = get(h(1),'position');

        % Initialise px and py:
        if ( nx > 1 )
            px = p(1,1) + p(1,3);
            py = p(1,2);
        else
            % Need to treat single column as a special case.
            px = leftMargin;
            py = p(1,2) - p(1,4) - fudge1;
        end

        % Loop over the remaining windows:
        for k = 2:numel(h)
            % Set to the new px and py values
            p(k,1:2) = [px, py];
            set(h(k), 'position', p(k,:));
            if ( mod(k, nx) )
                % Move to the next column:
                px = px + p(k-1,3);
            else
                % Move to the next row. 
                px = leftMargin;
                % Update py:
                if ( k == nx )
                    % First row is a special case (not sure why!)
                    py = py - p(k-1,4) - fudge1;
                else
                    py = py - p(k-1,4) - .5*fudgeN;
                end
            end        
        end

        % Bring all the windows to the front:
        for k = 1:nh
            figure(h(k))
        end
        end
    end
    
    methods (Static, Access = private)
        function plotSingleFactor_(xv, zv, error, estimates, mask_ext, geom, levels, FS, estuni)
            % PLOTSINGLEFACTOR_ Core single-factor error contour plot (private helper)
            %
            % Plots a filled contour of the quadrature error overlaid with line
            % contours of the error estimate in the current axes.
            %
            % Inputs:
            %   xv, zv    - Coordinate vectors
            %   error     - [M×1] error values at exterior points
            %   estimates - [M×1] error estimates at exterior points
            %   mask_ext  - Logical mask for exterior points
            %   geom      - AxsymGeometry object
            %   levels    - Contour levels
            %   FS        - Font size
            %   estuni    - (optional) [M×1] uniform estimates for green dashed overlay
            
            if nargin < 9
                estuni = [];
            end
            
            M = numel(xv);
            N = numel(zv);
            
            % Reshape error to grid (log10 scale)
            err_full = NaN(numel(mask_ext), 1);
            err_full(mask_ext) = log10(abs(error) + eps);
            Err = reshape(err_full, M, N);
            
            % Reshape estimates to grid
            est_full = NaN(numel(mask_ext), 1);
            est_full(mask_ext) = log10(abs(estimates) + eps);
            Est = reshape(est_full, M, N);
            
            % Clamp to minimum level
            minLevel = levels(end);
            Err(Err < minLevel) = minLevel;
            Est(Est < minLevel) = minLevel;
            
            % Handle Inf values
            Err(isinf(Err)) = -16;
            Est(isinf(Est)) = -16;
            
            % Meshgrid for plotting
            [X, Z] = meshgrid(xv, zv);
            
            % Filled contour of error
            contourf(X, Z, Err, levels, 'LineColor', 'none');
            hold on;
            
            % Line contour of estimates (black solid)
            [C, h] = contour(X, Z, Est, levels(2:end), '-k', 'LineWidth', 2);
            clabel(C, h, 'LabelSpacing', 1150, 'FontSize', 14);
            
            legendEntries = {'Quadrature Error', 'Error Estimate'};
            
            % Optional: uniform estimates overlay (green dashed)
            if ~isempty(estuni)
                estuni_full = NaN(numel(mask_ext), 1);
                estuni_full(mask_ext) = log10(abs(estuni) + eps);
                Estuni = reshape(estuni_full, M, N);
                Estuni(Estuni < minLevel) = minLevel;
                Estuni(isinf(Estuni)) = -16;
                
                [C2, h2] = contour(X, Z, Estuni, levels, '--g', 'LineWidth', 2);
                clabel(C2, h2, 'LabelSpacing', 300, 'Color', 'green', 'FontSize', 14);
                legendEntries{3} = 'Uniform Estimate';
            end
            
            % Draw geometry outline
            quadest.util.Plotting.drawGeometryOutline(geom);
            
            view(0, 90);
            axis equal;
            clim([levels(end), levels(1)]);
            xlim([min(xv), max(xv)]);
            ylim([min(zv), max(zv)]);
            
            % Colorbar
            cbar = colorbar('Ticks', fliplr(levels), 'TickLabelInterpreter', 'latex', 'FontSize', FS);
            xlabel(cbar, '$\log_{10}(\textrm{Absolute error})$', 'FontSize', FS, 'Interpreter', 'latex');
            
            % Colormap
            colormap(gca, quadest.util.Plotting.divergingColormap(length(levels) - 1));
            
            % Legend
            lgd = legend(legendEntries, 'Interpreter', 'latex', 'Location', 'southwest');
            lgd.FontSize = 15;
            
            xlabel('$x$', 'Interpreter', 'latex', 'FontSize', FS);
            ylabel('$z$', 'Interpreter', 'latex', 'FontSize', FS);
            set(gca, 'FontSize', FS);
            hold off;
        end
        function rgb = divergingMapValue(s, rgb1, rgb2)
            % Interpolate single value in diverging colormap
            
            lab1 = quadest.util.Plotting.rgbToLab(rgb1);
            lab2 = quadest.util.Plotting.rgbToLab(rgb2);
            
            msh1 = quadest.util.Plotting.labToMsh(lab1);
            msh2 = quadest.util.Plotting.labToMsh(lab2);
            
            % Insert white midpoint for saturated distinct colors
            if msh1(2) > 0.05 && msh2(2) > 0.05 && ...
               quadest.util.Plotting.angleDiff(msh1(3), msh2(3)) > 0.33*pi
                Mmid = max([msh1(1), msh2(1), 88.0]);
                if s < 0.5
                    msh2 = [Mmid, 0.0, 0.0];
                    s = 2.0 * s;
                else
                    msh1 = [Mmid, 0.0, 0.0];
                    s = 2.0 * s - 1.0;
                end
            end
            
            % Adjust hue for unsaturated colors
            if msh1(2) < 0.05 && msh2(2) > 0.05
                msh1(3) = quadest.util.Plotting.adjustHue(msh2, msh1(1));
            elseif msh2(2) < 0.05 && msh1(2) > 0.05
                msh2(3) = quadest.util.Plotting.adjustHue(msh1, msh2(1));
            end
            
            % Interpolate
            mshTmp = (1 - s) * msh1 + s * msh2;
            
            % Convert back
            labTmp = quadest.util.Plotting.mshToLab(mshTmp);
            rgb = quadest.util.Plotting.labToRgb(labTmp);
        end
        
        function msh = labToMsh(lab)
            L = lab(1); a = lab(2); b = lab(3);
            M = sqrt(L^2 + a^2 + b^2);
            s_val = (M > 0.001) * acos(L / max(M, eps));
            h = (s_val > 0.001) * atan2(b, a);
            msh = [M, s_val, h];
        end
        
        function lab = mshToLab(msh)
            M = msh(1); s = msh(2); h = msh(3);
            L = M * cos(s);
            a = M * sin(s) * cos(h);
            b = M * sin(s) * sin(h);
            lab = [L, a, b];
        end
        
        function adiff = angleDiff(a1, a2)
            v1 = [cos(a1), sin(a1)];
            v2 = [cos(a2), sin(a2)];
            adiff = acos(max(-1, min(1, dot(v1, v2))));
        end
        
        function h = adjustHue(msh, unsatM)
            if msh(1) >= unsatM - 0.1
                h = msh(3);
            else
                hueSpin = msh(2) * sqrt(unsatM^2 - msh(1)^2) / (msh(1) * sin(msh(2)) + eps);
                if msh(3) > -0.3 * pi
                    h = msh(3) + hueSpin;
                else
                    h = msh(3) - hueSpin;
                end
            end
        end
        
        function lab = rgbToLab(rgb)
            xyz = quadest.util.Plotting.rgbToXyz(rgb);
            lab = quadest.util.Plotting.xyzToLab(xyz);
        end
        
        function rgb = labToRgb(lab)
            xyz = quadest.util.Plotting.labToXyz(lab);
            rgb = quadest.util.Plotting.xyzToRgb(xyz);
        end
        
        function xyz = rgbToXyz(rgb)
            r = rgb(1); g = rgb(2); b = rgb(3);
            
            % sRGB gamma correction
            if r > 0.04045, r = ((r + 0.055) / 1.055)^2.4;
            else, r = r / 12.92; end
            if g > 0.04045, g = ((g + 0.055) / 1.055)^2.4;
            else, g = g / 12.92; end
            if b > 0.04045, b = ((b + 0.055) / 1.055)^2.4;
            else, b = b / 12.92; end
            
            x = r * 0.4124 + g * 0.3576 + b * 0.1805;
            y = r * 0.2126 + g * 0.7152 + b * 0.0722;
            z = r * 0.0193 + g * 0.1192 + b * 0.9505;
            xyz = [x, y, z];
        end
        
        function lab = xyzToLab(xyz)
            ref_X = 0.9505; ref_Y = 1.000; ref_Z = 1.089;
            var_X = xyz(1) / ref_X;
            var_Y = xyz(2) / ref_Y;
            var_Z = xyz(3) / ref_Z;
            
            if var_X > 0.008856, var_X = var_X^(1/3);
            else, var_X = 7.787 * var_X + 16.0/116.0; end
            if var_Y > 0.008856, var_Y = var_Y^(1/3);
            else, var_Y = 7.787 * var_Y + 16.0/116.0; end
            if var_Z > 0.008856, var_Z = var_Z^(1/3);
            else, var_Z = 7.787 * var_Z + 16.0/116.0; end
            
            L = 116 * var_Y - 16;
            a = 500 * (var_X - var_Y);
            b = 200 * (var_Y - var_Z);
            lab = [L, a, b];
        end
        
        function xyz = labToXyz(lab)
            L = lab(1); a = lab(2); b = lab(3);
            
            var_Y = (L + 16) / 116;
            var_X = a / 500 + var_Y;
            var_Z = var_Y - b / 200;
            
            if var_Y^3 > 0.008856, var_Y = var_Y^3;
            else, var_Y = (var_Y - 16.0/116.0) / 7.787; end
            if var_X^3 > 0.008856, var_X = var_X^3;
            else, var_X = (var_X - 16.0/116.0) / 7.787; end
            if var_Z^3 > 0.008856, var_Z = var_Z^3;
            else, var_Z = (var_Z - 16.0/116.0) / 7.787; end
            
            ref_X = 0.9505; ref_Y = 1.000; ref_Z = 1.089;
            xyz = [ref_X * var_X, ref_Y * var_Y, ref_Z * var_Z];
        end
        
        function rgb = xyzToRgb(xyz)
            x = xyz(1); y = xyz(2); z = xyz(3);
            
            r = x *  3.2406 + y * -1.5372 + z * -0.4986;
            g = x * -0.9689 + y *  1.8758 + z *  0.0415;
            b = x *  0.0557 + y * -0.2040 + z *  1.0570;
            
            % sRGB gamma correction
            if r > 0.0031308, r = 1.055 * r^(1/2.4) - 0.055;
            else, r = 12.92 * r; end
            if g > 0.0031308, g = 1.055 * g^(1/2.4) - 0.055;
            else, g = 12.92 * g; end
            if b > 0.0031308, b = 1.055 * b^(1/2.4) - 0.055;
            else, b = 12.92 * b; end
            
            % Clip to valid range
            maxVal = max([r, g, b]);
            if maxVal > 1.0
                r = r / maxVal;
                g = g / maxVal;
                b = b / maxVal;
            end
            rgb = [max(0, r), max(0, g), max(0, b)];
        end
    end
end
