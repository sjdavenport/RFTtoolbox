function int = integrateTriangulation(trep, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function computing the exact integral of the linear interpolant of a
% function given on a triangulation.
% This is done by evaluating the areas of all the triangles and
% then using these areas as weights for the midpoint(mean)-value
% on each triangle.
% The code is due to:
% https://stackoverflow.com/questions/23688669/how-to-integrate-over-a-discrete-2d-surface-in-matlab
%
% Input: 
%   trep (triangulation): e.g., output of delaunayTriangulation.m can be
%                         used
%   z: values of the function on the nodes
%
% Output:
%   int (numeric): integral of the function defined on the triangulation
%
% Example:
% % Create a set of points in 2D space
% P = [0 0; 1 0; 1 1; 0 1; 0.5 0.5];
% 
% % Delaunay triangulate the points
% trep = delaunayTriangulation(P);
% 
% % Define values of a function on the nodes
% z = [1; 2; 3; 4; 5];
% 
% % Compute the integral of the function defined on the triangulation
% int = integrateTriangulation(trep, z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = trep.Points; T = trep.ConnectivityList;
d21 = P(T(:,2),:)-P(T(:,1),:);
d31 = P(T(:,3),:)-P(T(:,1),:);
areas = abs(1/2*(d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1)));
int = areas'*mean(z(T),2);

end