% Plot errors, categorized by different quadratute methods used
function plot_quaderrs(precomp, body, xpts, u, uref, FS)
% Distance plot
[maskS3Q, distS3Q, X2surf] = quadregions(precomp, body, xpts.xout, xpts.mask);
xdist = zeros(size(xpts.xt, 1));
xdist(xpts.mask) = X2surf;

% Extract masks
Smask = maskS3Q.SQ_mask;
Dirmask= maskS3Q.direct_mask;
upmask = maskS3Q.upsamp_mask;

UerrMask = sqrt( sum( (u(Smask, :) - uref(Smask, :) ).^2, 2) );
UerrMaskDir = sqrt( sum( (u(Dirmask, :) - uref(Dirmask, :) ).^2, 2) );
UerrMaskUp = cell(1,numel(maskS3Q.upfac)-1);
for k = 1:numel(maskS3Q.upfac)
    UerrMaskUp{k} = sqrt( sum( (u(upmask{k}, :) - uref(upmask{k}, :) ).^2, 2) );
end

% Extract number of target points in each quadrature region
nbr_targets_direct = sum(maskS3Q.direct_mask,1);
nbr_targets_upsamp = zeros(1,numel(maskS3Q.upfac));
for k = 1:numel(maskS3Q.upfac)
    nbr_targets_upsamp(k) = sum(maskS3Q.upsamp_mask{k},1);
end
nbr_targets_sq = sum(Smask);
fprintf('   # Direct points: %d \n', nbr_targets_direct);
fprintf('   # Upsamp points: %d \n', sum(nbr_targets_upsamp));
for k = 1:numel(maskS3Q.upfac)
fprintf('   #   Factor %i:      %d \n', maskS3Q.upfac(k),nbr_targets_upsamp(k));
end
fprintf('   # SQ points:     %d \n', nbr_targets_sq);

% Legend for upcoming plots
legents = cell(1,numel(maskS3Q.upfac)+1);
for k = 1:numel(maskS3Q.upfac)
    legents{k} = sprintf('Upsampling = %i',maskS3Q.upfac(k));
end
legents{1} = sprintf('Direct Quadrature'); 
legents{numel(maskS3Q.upfac)+1} = sprintf('Special Quadrature'); 

%% Distance to particle vs S3Q SOL error, with the different quadrature 
% regions marked
figure();clf;
ax1 = gca;
hold on
plot(xdist(distS3Q.direct_mask), UerrMaskDir, 'k.')
ax1.ColorOrderIndex = 2;
for k = 2:numel(maskS3Q.upfac)
    plot(xdist(distS3Q.upsamp_mask{k}), UerrMaskUp{k}, '.')
end
ax1.ColorOrderIndex = 1;
plot(xdist(distS3Q.SQ_mask), UerrMask, '.')
yline(precomp.tol,'k--','Tolerance')
xlabel('Distance from particle')
ylabel('Error')
set(gca, 'yscale', 'log')
set(gca, 'xscale', 'log')
title('Solution Error','fontsize',12)
grid on;
legend(legents, 'Interpreter', 'LaTex', 'location', 'best');

figure();clf;
ax2 = gca;
hold on
rquad = zeros(numel(xpts.xv), numel(xpts.zv));

rmask = reshape(distS3Q.direct_mask, numel(xpts.xv), numel(xpts.zv));
rquad(rmask) = numel(maskS3Q.upfac);
for k = 2:numel(maskS3Q.upfac)
  rumask = reshape(distS3Q.upsamp_mask{k}, numel(xpts.xv), numel(xpts.zv));
  rquad(rumask) = numel(maskS3Q.upfac) - k + 1;
end
rmask = reshape(distS3Q.SQ_mask, numel(xpts.xv), numel(xpts.zv));
rquad(rmask) = 0;
rquad = -abs(rquad)-0.5;

levels = [-0.5:-1:-numel(maskS3Q.upfac)-0.5];


rx    = reshape(xpts.xt, numel(xpts.xv), numel(xpts.zv), 3);
contourf(rx(:, :, 1), rx(:, :, 3), rquad, levels,'LineColor', 'none')
hold on
meshing(body) 
view(0, 90)

clim([-length(levels) 0])
axis([xpts.xv(1) xpts.xv(end) xpts.zv(1) xpts.zv(end)])
cbar = colorbar;
cbar = colorbar('Ticks',fliplr(levels),'TickLabelInterpreter','latex','FontSize',FS);
axis equal 
cbar.TickLabels = legents;

colormap(dcmap_redblue(numel(legents)));    %(8)

% title(sprintf('levels = 10^{%s}', num2str(levels)))
xlabel('$x$','interpreter','latex','fontsize',FS);
ylabel('$z$','interpreter','latex','fontsize',FS);
zlabel('$y$','interpreter','latex','fontsize',FS);
set(gca,'FontSize',FS);
set(gca,'FontSize',FS)
xlim([min(xpts.xv) max(xpts.xv)])
ylim([min(xpts.zv) max(xpts.zv)])


end