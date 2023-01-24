function [ out ] = ballker( x, rad, leq )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% x         a D by nvalues matrix where D is the number of dimensions and
%           nvalues is the number of points at which to evaluate the kernel
% rad       half the size of the sides of the box
% leq       0/1 determines whether or not to include the boundary in the
%           support of the kernel or not. 0 (default) excludes it.
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 2
    rad = 1;
end
if nargin < 3
    leq = 0;
end


%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
x_sos = sum(x.^2, 1);

%%  Main Function Loop
%--------------------------------------------------------------------------
if leq
    out = x_sos <= rad^2;
else
    out = x_sos < rad^2;
end

end

