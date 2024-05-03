function out = boxker( x, rad, leq )
% boxker( x, rad )
%--------------------------------------------------------------------------
% ARGUMENTS
% x         a D by nvalues matrix where D is the number of dimensions and
%           nvalues is the number of points at which to evaluate the kernel
% rad       half the size of the sides of the box
% leq       0/1 determines whether or not to include the boundary in the
%           support of the kernel or not. 0 (default) excludes it.
%--------------------------------------------------------------------------
% OUTPUT
% out       vector of 0s and 1s each element of which indicates whether or 
%           not the elements are within the box of radius rad
%--------------------------------------------------------------------------
% EXAMPLES
% boxker([1/2,1/4])
% boxker([1/2,1/4]')
% boxker([1/2,1/4, 2/3])
% boxker([1/2,1/4;2/3,0.2])
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 2
    rad = 1/2;
end
if nargin < 3
    leq = 0;
end

if leq
    out = all(x <= rad,1).*all(x >= -rad,1);
else
    out = all(x < rad,1).*all(x > -rad,1);
end

% if leq
%     out = all(x <= rad).*all(x >= -rad);
% else
%     out = all(x < rad).*all(x > -rad);
% end

end

% lessthan = sum(x <= rad);
% greaterthan = sum(x >= -rad);
% 
% test = (lessthan+greaterthan) == 2*D
