% Color map for contour plots
function map = dcmap_redblue(m)
%DCMAP_REDBLUE    Red and blue diverging colormap
%   DCMAP_REDBLUE(M) returns an M-by-3 matrix containing an HSV colormap.
%   DCMAP_REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
if nargin < 1, m = size(get(gcf,'colormap'),1); end
map = diverging_map(linspace(0,1,m), [0.230, 0.299, 0.754], [0.706, 0.016, 0.150]);
end
