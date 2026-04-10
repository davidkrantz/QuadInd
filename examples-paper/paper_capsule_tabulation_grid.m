% paper_capsule_tabulation_grid.m - Paper figure for capsule geometry

clear; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
addpath(projectDir);

savefig = 0;
FS = 16;

nth = 40;
nph = 60;

geom = quadest.geometry.Capsule('R', 1, 'L', 6, 'kappa', 3);

% Reconstruct the tabulation grid coordinates 
% (matches UniformEstimateBuilder.buildTabGrid)
cfg = quadest.util.Config();
maxExtent = 2 * (geom.maxRadius() + geom.maxHeight());

% First quarter plane only: rxy >= 0, z >= 0
rxy_tab = linspace(0, maxExtent, cfg.ntab - 1);
z_tab   = linspace(0, maxExtent, cfg.nztab - 1);
rxy_tab = [-rxy_tab(2), rxy_tab];
z_tab = [-z_tab(2), z_tab];

[RXY, ZABS] = ndgrid(rxy_tab, z_tab);
evalPoints = [RXY(:), zeros(numel(RXY), 1), ZABS(:)];

extMask = geom.isExterior(evalPoints);

figure('DefaultAxesFontSize',FS);
plot3(zeros(size(RXY(extMask))), RXY(extMask), ZABS(extMask), '.', ...
    'Color', 'b', 'MarkerSize', 6);
hold on;
quadest.util.Plotting.plotSurfaceMesh(geom);
axis equal;
ylim([rxy_tab(1), maxExtent/2]);
zlim([z_tab(1), maxExtent/2]);
ylabel('$\rho$', 'Interpreter', 'latex', 'FontSize', FS);
zlabel('$\zeta$', 'Interpreter', 'latex', 'FontSize', FS);
xticklabels([]); yticklabels([]); zticklabels([]);
set(gca, 'FontSize', FS);
view(90,0);

if savefig
    if ~exist('../figs', 'dir'); mkdir('../figs'); end
    disp('saving figures...');
    exportgraphics(figure(1),'../figs/capsule_tabulation_grid.pdf','Resolution',400);
    disp('sucessfully saved figures');
end
