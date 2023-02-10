function [ points, xlevellocs, ylevellocs ] = sample_boundary( f, xvals, nintersections, f_grid )
% sample_boundary( f, xvals, nitersections ) evaluates the function on a 
% grid defined by the x and y values, and then uses the fmincon optimization 
% method with the sqp algorithm to find the x and y coordinates of the 
% sampled points on the function boundary. The results are returned in the 
% form of matrices of x and y coordinates.
%--------------------------------------------------------------------------
% ARGUMENTS
% f: a function handle of the form f(x,y)
% xvals: a cell array of length 2, containing vectors of the x and y values
%       over which to sample the function nitersections: number of sections 
%       in which the sampling is divided
%--------------------------------------------------------------------------
% OUTPUT
% points: a D by npoints array 
% xlevellocs: a matrix of size (length(xvals{1}), nitersections) containing
%         the x-coordinates of the sampled points on the function boundary
% ylevellocs: a matrix of size (length(xvals{2}), nitersections) containing
%          the y-coordinates of the sampled points on the function boundary
%--------------------------------------------------------------------------
% EXAMPLES
% f = @(x,y) x.^2 + y.^2 - 1;
% xvals = {-1.5:0.1:1.5, -1.5:0.1:1.5};
% [ xlevellocs, ylevellocs, points ] = sample_boundary( f, xvals, 4 )
% 
% subplot(2,1,1)
% for J = 1:4
%    plot(xvals{1}, ylevellocs(:,J), '*')
%    hold on
% end
% subplot(2,1,2)
% for J = 1:4
%    plot(xlevellocs(:,J),xvals{2}, '*')
%    hold on
% end
%
% plot(points(1,:), points(2,:), '*')
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
[x1, x2] = meshgrid(xvals{1}, xvals{2});
if ~exist('f_grid', 'var')
    try
        f_grid = f(x1,x2);
    catch
        f_grid = reshape(f(x1(:)', x2(:)'), [length(xvals{2}), length(xvals{1})]);
    end
end

% Set the options for the fmincon optimization
options = optimoptions(@fmincon,'Display','off'); % Ensures that no output
options.Algorithm = 'sqp';
options = optimset('Display','off');

ylevellocs = zeros([length(xvals{1}), nintersections]);
D = 2;
ypoints = zeros([D, length(xvals{1})*nintersections]);

for I = 1:length(xvals{1})
    loader(I, length(xvals{1}), 'Points in x direction, percent done:')
    sample_grid_f = f_grid(:,I);
    [~, init_est_index] = sort(abs(sample_grid_f));
    init_zero_ests = xvals{2}(init_est_index(1:nintersections));
    for J = 1:nintersections
        ylevellocs(I, J) = fzero(@(y) f(xvals{1}(I), y), init_zero_ests(J), options);
    end
end
for J = 1:nintersections
    ypoints(1,(1:length(xvals{1})) + (J-1)*length(xvals{1})) = xvals{1};
    ypoints(2,(1:length(xvals{1})) + (J-1)*length(xvals{1})) = ylevellocs(:,J)';
end

xlevellocs = zeros([length(xvals{2}), nintersections]);
xpoints = zeros([D, length(xvals{2})*nintersections]);

for I = 1:length(xvals{2})
    loader(I, length(xvals{2}), 'Points in y direction, percent done:')
    sample_grid_f = f_grid(I,:);
    [~, init_est_index] = sort(abs(sample_grid_f));
    init_zero_ests = xvals{1}(init_est_index(1:nintersections));
    for J = 1:nintersections
        xlevellocs(I, J) = fzero(@(x) f(x, xvals{2}(I)), init_zero_ests(J), options);
    end
end

for J = 1:nintersections
    xpoints(1,(1:length(xvals{2})) + (J-1)*length(xvals{2})) = xlevellocs(:,J)';
    xpoints(2,(1:length(xvals{2})) + (J-1)*length(xvals{2})) = xvals{2};
end

points = [ypoints, xpoints];

plot(points(1,:), points(2,:), '*')

end

