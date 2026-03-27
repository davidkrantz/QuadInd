% Determine different quadrature regions used
function [masks, dist, xt2surf] = quadregions(precomp, body_grid, x_ext, mask_ext)

% Extract quadrature types for each particles
quadregchk = fieldnames(precomp);
% Remove upsampfac from list if it exists
quadregs = {};
for iq = 1:length(quadregchk)
  quadchk = quadregchk{iq};
  if ~strcmp(quadchk, 'upsampfac') && ~strcmp(quadchk, 'tol')
    quadregs{end+1} = quadchk;
    % quadregs(iq) = [];
  end
end

% Determine the distance between the target point and closest
% discretization point (locally)
[~, xt2surf] = dsearchn( body_grid.x,x_ext ) ;

masks = [];
for iq = 1:length(quadregs)
  % Upsampling must be split
  if ~strcmp(quadregs{iq}, 'upsamp_mask')
    % Create "masks" to contain the values (logical)
    % eval(sprintf('masks.%s = logical(zeros(size(mask_ext, 1), 1));', quadregs{iq}));
    eval(sprintf('masks.%s(mask_ext, 1) = precomp.%s;', quadregs{iq}, quadregs{iq}));
  else % Different upsampling regions
    upfac = precomp.upsampfac;
    masks.upfac = upfac;
    masks.upsamp_mask = cell(1,numel(upfac));
    for iu = 1:length(upfac)
      % Create "masks" to contain the values (logical)
      %eval(sprintf('masks.%s%d = logical(zeros(size(mask_ext, 1), bodies.num_bodies));', quadregs{iq}, upfac(iu)));
      tmp_mask = false(size(mask_ext, 1), 1);
      %eval(sprintf('masks.%s%d(mask_ext, ii) = bodies.bodies{1, ii}.quadmask.upsamp_mask{%d};', quadregs{iq}, upfac(iu), iu));
      tmp_mask(mask_ext, 1) = precomp.upsamp_mask{iu};
      masks.upsamp_mask{iu} = tmp_mask;
    end
  end
end

dist = [];
quadmasks = quadregs(:);
for iq = 1:length(quadmasks)
  % Create mask, to determine if current quadrature method was used on either point
  if ~strcmp(quadregs{iq}, 'upsamp_mask')
    eval( sprintf('dist.%s = (mask_ext + masks.%s )>1;', quadmasks{iq}, quadmasks{iq}) );
  else
    upfac = precomp.upsampfac;
    dist.upsamp_mask = cell(1,numel(upfac));
    for iu = 1:length(upfac)
      dist.upsamp_mask{iu} = (mask_ext + masks.upsamp_mask{iu} )>1;
    end
  end
end
end
