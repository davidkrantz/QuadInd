% paper_quadrature_regimes_spheroid.m
% Conceptual top-down figure: quadrature regimes near two prolate spheroids.

clear; close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
addpath(projectDir);

savefig = 0;

% Prolate spheroids viewed from above. The major axis lies in the target
% plane; rotations make the two-body configuration less symmetric.
a = 1.0;      % minor semi-axis
c = 2.25;     % major semi-axis
nth = 40;
nph = 40;
centers = [-1.90, -0.70; ...
            1.80,  0.66];
anglesDeg = [15; -24];

% Conceptual regime boundaries. At overlaps, the target plane is colored by
% the maximum required quadrature: S3Q > upsampled > standard.
sS3Q = 1.22;
sUpsampled = 1.72;

bg = 0.9*[0.930, 0.935, 0.925];              % warm light gray
colStandard = [0.250, 0.300, 0.340];
colUpsampled = [0.000, 0.450, 0.700];
colS3Q = [0.780, 0.320, 0.120];
fillUpsampled = [0.730, 0.875, 0.930];
fillS3Q = [0.955, 0.760, 0.650];

fig = figure( ...
    'Color', 'w', ...
    'InvertHardcopy', 'off', ...
    'Renderer', 'painters');
ax = axes(fig, 'Position', [0.02, 0.02, 0.96, 0.96]);
hold(ax, 'on');
set(ax, 'Color', bg);

% Draw the target plane as one classified field. This avoids darkening or
% ambiguity where regions from the two spheroids overlap.
if savefig
    xv = linspace(-7.2, 7.2, 5000);
    yv = linspace(-4.7, 4.7, 5000);
else
    xv = linspace(-7.2, 7.2, 500);
    yv = linspace(-4.7, 4.7, 500);
end
[X, Y] = meshgrid(xv, yv);
regime = classifyTargets(X, Y, a, c, centers, anglesDeg, sS3Q, sUpsampled);
imagesc(ax, xv, yv, regime);
set(ax, 'YDir', 'normal');
colormap(ax, [bg; fillUpsampled; fillS3Q]);
clim(ax, [1, 3]);

% Region outlines are drawn as unions, not individual spheroid ellipses.
drawUnionOutline(ax, X, Y, regime >= 2, colUpsampled, 1.45);
drawUnionOutline(ax, X, Y, regime >= 3, colS3Q, 1.45);

% Projected surface discretization on each spheroid.
for ibody = 1:size(centers, 1)
    drawProjectedSpheroidMesh(ax, a, c, centers(ibody,:), ...
        anglesDeg(ibody), nth, nph);
end

axis(ax, 'equal');
axis(ax, 'off');
xlim(ax, [min(xv), max(xv)]);
ylim(ax, [min(yv), max(yv)]);

% Labels and leader lines, following the Bagge--Tornberg illustration style.
annotation(fig, 'textbox', [0.12, 0.800, 0.55, 0.08], ...
    'String', 'standard quadrature', ...
    'Interpreter', 'latex', ...
    'Color', colStandard, ...
    'FontSize', 23, ...
    'FontWeight', 'bold', ...
    'EdgeColor', 'none', ...
    'BackgroundColor', 'none', ...
    'FitBoxToText', 'on');

annotation(fig, 'arrow', [0.750, 0.600], [0.310, 0.460], ...
    'Color', colS3Q, ...
    'LineWidth', 1.35, ...
    'HeadLength', 8, ...
    'HeadWidth', 8);
annotation(fig, 'textbox', [0.750, 0.260, 0.18, 0.08], ...
    'String', 'S3Q', ...
    'Interpreter', 'latex', ...
    'Color', colS3Q, ...
    'FontSize', 22, ...
    'FontWeight', 'bold', ...
    'EdgeColor', 'none', ...
    'BackgroundColor', 'none', ...
    'FitBoxToText', 'on');

annotation(fig, 'arrow', [0.220, 0.180], [0.220, 0.330], ...
    'Color', colUpsampled, ...
    'LineWidth', 1.35, ...
    'HeadLength', 8, ...
    'HeadWidth', 8);
annotation(fig, 'textbox', [0.145, 0.150, 0.62, 0.08], ...
    'String', 'upsampled quadrature', ...
    'Interpreter', 'latex', ...
    'Color', colUpsampled, ...
    'FontSize', 21, ...
    'FontWeight', 'bold', ...
    'EdgeColor', 'none', ...
    'BackgroundColor', 'none', ...
    'FitBoxToText', 'on');

if savefig
    figDir = fullfile(projectDir, 'figs');
    if ~exist(figDir, 'dir'); mkdir(figDir); end
    outFile = fullfile(figDir, 'quadrature_regimes_spheroid.pdf');
    exportgraphics(fig, outFile, 'Resolution', 600, 'BackgroundColor', bg);
    fprintf('Saved %s\n', outFile);
end

function regime = classifyTargets(X, Y, a, c, centers, anglesDeg, sS3Q, sUpsampled)
    regime = ones(size(X));  % 1 = standard, 2 = upsampled, 3 = S3Q
    for ibody = 1:size(centers, 1)
        [Xloc, Yloc] = toLocalPlane(X, Y, centers(ibody,:), anglesDeg(ibody));
        inUpsampled = (Xloc./(sUpsampled*c)).^2 + (Yloc./(sUpsampled*a)).^2 <= 1;
        inS3Q = (Xloc./(sS3Q*c)).^2 + (Yloc./(sS3Q*a)).^2 <= 1;
        regime(inUpsampled) = max(regime(inUpsampled), 2);
        regime(inS3Q) = max(regime(inS3Q), 3);
    end
end

function drawUnionOutline(ax, X, Y, mask, color, lineWidth)
    contour(ax, X, Y, double(mask), [0.5, 0.5], ...
        'Color', color, ...
        'LineWidth', lineWidth);
end

function drawEllipseOutline(ax, a, c, scale, center, angleDeg, color, lineWidth)
    t = linspace(0, 2*pi, 700);
    xLocal = scale * c * cos(t);
    yLocal = scale * a * sin(t);
    [x, y] = toGlobalPlane(xLocal, yLocal, center, angleDeg);
    plot(ax, x, y, '-', 'Color', color, 'LineWidth', lineWidth);
end

function drawProjectedSpheroidMesh(ax, a, c, center, angleDeg, nth, nph)
    faceColor = 0.8 * [1, 1, 1];
    edgeColor = [0.3, 0.3, 0.3];

    theta = linspace(0, pi, nth);
    phi = linspace(0, 2*pi, nph);
    [Theta, Phi] = ndgrid(theta, phi);

    xLocal = c * cos(Theta);
    yLocal = a * sin(Theta) .* cos(Phi);
    [x, y] = toGlobalPlane(xLocal, yLocal, center, angleDeg);

    % Light fill under the wire mesh.
    drawEllipseOutline(ax, a, c, 1.0, center, angleDeg, edgeColor, 1.15);
    t = linspace(0, 2*pi, 500);
    [xf, yf] = toGlobalPlane(c*cos(t), a*sin(t), center, angleDeg);
    patch(ax, xf, yf, faceColor, ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.96);

    % Projected latitude and longitude curves.
    for itheta = 2:2:(nth-1)
        plot(ax, x(itheta,:), y(itheta,:), '-', ...
            'Color', edgeColor, ...
            'LineWidth', 0.42);
    end
    for iphi = 1:2:nph
        plot(ax, x(:,iphi), y(:,iphi), '-', ...
            'Color', edgeColor, ...
            'LineWidth', 0.42);
    end
    drawEllipseOutline(ax, a, c, 1.0, center, angleDeg, edgeColor, 1.15);
end

function [xGlobal, yGlobal] = toGlobalPlane(xLocal, yLocal, center, angleDeg)
    ca = cosd(angleDeg);
    sa = sind(angleDeg);
    xGlobal = center(1) + ca .* xLocal - sa .* yLocal;
    yGlobal = center(2) + sa .* xLocal + ca .* yLocal;
end

function [xLocal, yLocal] = toLocalPlane(xGlobal, yGlobal, center, angleDeg)
    ca = cosd(angleDeg);
    sa = sind(angleDeg);
    dx = xGlobal - center(1);
    dy = yGlobal - center(2);
    xLocal = ca .* dx + sa .* dy;
    yLocal = -sa .* dx + ca .* dy;
end
