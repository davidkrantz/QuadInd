% Plot contour plots of error and error estimates
function plot_estcontour(xpts, body, density, estplt, upfac, u_ref, FS)
uni_errest = estplt.Errest_uni;
estuni = NaN*zeros(size(xpts.xt, 1), 1);
estuni(xpts.mask) = log10(abs(sum(uni_errest, 2)));
dir_errest = estplt.estup_stor{upfac};
estdir = NaN*zeros(size(xpts.xt, 1), 1);
estdir(xpts.mask) = log10(abs(sum(dir_errest, 2)));


fgrid = upsample(body, upfac);
fdensity = fgrid.interpolate(density);
fnorm = bsxfun(@times, fgrid.n, fgrid.w);
u_dir = stresslet(xpts.xt, fgrid.x, fnorm, fdensity);

% Error
UdirerrSOL = sqrt( sum( (u_dir - u_ref).^2, 2) );

% Plot error estimates
plotest(body, xpts.xt, xpts.xv, xpts.zv, xpts.X, xpts.Z, UdirerrSOL, estdir, xpts.mask, FS)

% Show uniform (unmodified) estimates too
plotest(body, xpts.xt, xpts.xv, xpts.zv, xpts.X, xpts.Z, UdirerrSOL, estdir, xpts.mask, FS, estuni)

end

%% plotest: Error estimates plot
function plotest(bodies, x, xv, zv, X, Z, Uerr, est, mask, FS, estuni)
xx =reshape(x,numel(xv),numel(zv),3);



% est = NaN*zeros(size(x, 1), 1);

err           = log10(abs(Uerr));
% est(mask_ext) = log10(abs(sum(errest, 2)));
Err           = reshape(err.*mask,numel(xv),numel(zv));
Est           = reshape(est.*mask,numel(xv),numel(zv));

if nargin >10
Estuni     = reshape(estuni.*mask,numel(xv),numel(zv));
end

figure()
clf
levels = [0:-2:-10];
[nn1,nn2]=size(Est);
for i=1:nn1
    for j=1:nn2
        if Est(i,j)<levels(end)
            Est(i,j)=levels(end);
        end
        if nargin >10
            if Estuni(i,j)<levels(end)
                Estuni(i,j)=levels(end);
            end
        end
        if Err(i,j)<levels(end)
            Err(i,j)=levels(end);
        end
    end
end
TF =find(isinf(Err));
for i=1:length(TF)
    [ind1,ind2]=ind2sub(size(Err),TF(i));
    Err(ind1,ind2)=-16;
end
contourf(xx(:,:,1),xx(:,:,3),Err, levels, 'LineColor','none');
Legend_ent{1} = sprintf('Direct Quadrature Error');
hold on

[C, h] = contour(X,Z,Est,levels(2:end), '-k', 'LineWidth',2);

if nargin >10
[C2, h2] = contour(X,Z,Estuni,levels, '--g', 'LineWidth',2);
end

meshing(bodies) 
view(0, 90)

Legend_ent{2} = sprintf('Error Estimate');
clim([levels(end) levels(1)])
clabel(C,h,'LabelSpacing',1150, 'FontSize',14)
if nargin >10
clabel(C2,h2,'LabelSpacing',300, 'Color','green', 'FontSize',14)
% Legend_ent{3} = sprintf('Uniform Estimate');
Legend_ent{3} = sprintf('SI Estimate');
end
axis([xv(1) xv(end) zv(1) zv(end)])
cbar = colorbar;
cbar = colorbar('Ticks',fliplr(levels),'TickLabelInterpreter','latex','FontSize',FS);
xlabel(cbar,'$\log_{10}(\textrm{Absolute error})$','fontsize',FS,'Interpreter','latex');
axis equal 
lgd2 = legend(Legend_ent, 'Interpreter', 'LaTex', 'location', 'southwest');
%lgd2 = legend(Legend_ent, 'Interpreter', 'LaTex', 'location', 'best');
lgd2.FontSize = 15;
colormap(dcmap_redblue(length(levels)-1));    %(8)

% title(sprintf('levels = 10^{%s}', num2str(levels)))
xlabel('$x$','interpreter','latex','fontsize',FS);
ylabel('$z$','interpreter','latex','fontsize',FS);
zlabel('$y$','interpreter','latex','fontsize',FS);
set(gca,'FontSize',FS);
set(gca,'FontSize',FS)
xlim([min(xv) max(xv)])
ylim([min(zv) max(zv)])

end