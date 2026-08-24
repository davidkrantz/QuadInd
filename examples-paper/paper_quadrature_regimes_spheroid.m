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

% Conceptual regime boundaries for one active spheroid at a time.
sS3Q = 1.22;
sUpsampled = 1.72;

bg = 0.9*[0.930, 0.935, 0.925];              % warm light gray
colUpsampled = [0.000, 0.450, 0.700];
colS3Q = [0.780, 0.320, 0.120];
fillUpsampled = [0.730, 0.875, 0.930];
fillS3Q = [0.955, 0.760, 0.650];

% Identical, tightly cropped panel geometry for side-by-side placement.
xLimits = [-6.05, 5.80];
yLimits = [-3.05, 3.15];
figureSizeCm = [8.6, 4.5];

if savefig
    xv = linspace(xLimits(1), xLimits(2), 5000);
    yv = linspace(yLimits(1), yLimits(2), 5000);
else
    xv = linspace(xLimits(1), xLimits(2), 500);
    yv = linspace(yLimits(1), yLimits(2), 500);
end
[X, Y] = meshgrid(xv, yv);

for activeBody = 1:size(centers, 1)
    fig = figure( ...
        'Units', 'centimeters', ...
        'Position', [2, 2, figureSizeCm], ...
        'Color', 'w', ...
        'InvertHardcopy', 'off', ...
        'Renderer', 'painters');
    set(fig, ...
        'PaperUnits', 'centimeters', ...
        'PaperSize', figureSizeCm, ...
        'PaperPosition', [0, 0, figureSizeCm]);
    ax = axes(fig, 'Position', [0, 0, 1, 1]);
    hold(ax, 'on');
    set(ax, 'Color', bg);

    % Classify targets using only the active spheroid. The other spheroid is
    % retained in the drawing solely as geometric context.
    regime = classifyTargetsForBody(X, Y, a, c, centers(activeBody,:), ...
        anglesDeg(activeBody), sS3Q, sUpsampled);
    imagesc(ax, xv, yv, regime);
    set(ax, 'YDir', 'normal');
    colormap(ax, [bg; fillUpsampled; fillS3Q]);
    clim(ax, [1, 3]);

    drawUnionOutline(ax, X, Y, regime >= 2, colUpsampled, 1.45);
    drawUnionOutline(ax, X, Y, regime >= 3, colS3Q, 1.45);

    for ibody = 1:size(centers, 1)
        drawProjectedSpheroidMesh(ax, a, c, centers(ibody,:), ...
            anglesDeg(ibody), nth, nph);
    end

    if activeBody == 1
        drawRegimeLabels(ax, colUpsampled, colS3Q);
    end

    axis(ax, 'equal');
    axis(ax, 'off');
    xlim(ax, xLimits);
    ylim(ax, yLimits);

    if savefig
        figDir = fullfile(projectDir, 'figs');
        if ~exist(figDir, 'dir'); mkdir(figDir); end
        bodyNames = {'left', 'right'};
        outFile = fullfile(figDir, sprintf( ...
            'quadrature_regimes_%s_spheroid.pdf', bodyNames{activeBody}));
        print(fig, outFile, '-dpdf', '-painters');
        fprintf('Saved %s\n', outFile);
    end
end

function regime = classifyTargetsForBody(X, Y, a, c, center, angleDeg, sS3Q, sUpsampled)
    regime = ones(size(X));  % 1 = standard, 2 = upsampled, 3 = S3Q
    [Xloc, Yloc] = toLocalPlane(X, Y, center, angleDeg);
    inUpsampled = (Xloc./(sUpsampled*c)).^2 + (Yloc./(sUpsampled*a)).^2 <= 1;
    inS3Q = (Xloc./(sS3Q*c)).^2 + (Yloc./(sS3Q*a)).^2 <= 1;
    regime(inUpsampled) = 2;
    regime(inS3Q) = 3;
end

function drawUnionOutline(ax, X, Y, mask, color, lineWidth)
    contour(ax, X, Y, double(mask), [0.5, 0.5], ...
        'Color', color, ...
        'LineWidth', lineWidth);
end

function drawRegimeLabels(ax, colUpsampled, colS3Q)
    colStandard = [0.34, 0.34, 0.34];
    fontSize = 9;

    xLabel = -5.72;
    yLabel = [2.65, 2.65-0.5, 2.65-1.0];
    labels = {'standard quadrature', 'upsampled quadrature', 'special quadrature (S3Q)'};
    colors = [colStandard; colUpsampled; colS3Q];

    for ilabel = 1:numel(labels)
        text(ax, xLabel, yLabel(ilabel), labels{ilabel}, ...
            'Color', colors(ilabel,:), ...
            'FontSize', fontSize, ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle');
    end
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
