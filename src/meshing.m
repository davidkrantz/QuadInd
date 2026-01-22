% Function to plot particle mesh
function meshing(body,varargin)
bdy = body;
mesh_grid(bdy.x,bdy.nth,bdy.nph,varargin{:})
end
%% mesh_grid: Creates plot of particle grid
function mesh_grid(xvec, m, n, varargin)
% mesh_grid(xvec, m, n)
% mesh_grid(xvec, m, n, q)
% mesh_grid(xvec, m, n, 'PropertyName', PropertyValue, ...)
%
% Mesh plot of grid
X = reshape(xvec(:,1),m,n);
Y = reshape(xvec(:,2),m,n);
Z = reshape(xvec(:,3),m,n);
% Wrap around in theta direction
X = [X, X(:,1)];
Y = [Y, Y(:,1)];
Z = [Z, Z(:,1)];
% Draw surface
% if numel(varargin) == 0
    surf(X, Z, 0.0*Y, 'FaceColor', 0.8*[1 1 1]);
% elseif numel(varargin) == 1
%     C = reshape(varargin{1},m,n);
%     C = [C, C(:,1)];
%     surf(X, Y, Z, C);
% else
%     surf(X, Y, Z, varargin{:});
% end
axis equal
end
