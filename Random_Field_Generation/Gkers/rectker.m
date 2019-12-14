function out = rectker( x, rad, leq )
% rectker( x, rad )
%--------------------------------------------------------------------------
% ARGUMENTS
% x
% rad
%--------------------------------------------------------------------------
% OUTPUT
% out
%--------------------------------------------------------------------------
% EXAMPLES
% rectker([1/2,1/4])
% rectker([1/2,1/4, 2/3])
% rectker([1/2,1/4;2/3,0.2])
%--------------------------------------------------------------------------
% AUTHOR: Samuel J. Davenport
if nargin < 2
    rad = 1/2;
end
if nargin < 3
    leq = 0;
end

if leq
    out = all(x <= rad).*all(x >= -rad);
else
    out = all(x < rad).*all(x > -rad);
end

end

% lessthan = sum(x <= rad);
% greaterthan = sum(x >= -rad);
% 
% test = (lessthan+greaterthan) == 2*D
