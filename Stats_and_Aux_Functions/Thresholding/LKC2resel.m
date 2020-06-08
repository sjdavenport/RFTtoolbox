function resel_vec = LKC2resel( LKCs, L_0 )
% LKC2RESEL( LKCs ) takes in a vector of LKCs and returns the corresponding
% vector of resels.
%--------------------------------------------------------------------------
% ARGUMENTS
% LKCs      a vector of the LKCs such that L(i) is L_i i.e. the ith LKC for
%           i = 1:D. Note that L_0 is not included!
% L_0       the 0th LKC i.e. L_0. If this is not provided as an input it is
%           taken to be 1
%--------------------------------------------------------------------------
% OUTPUT
% resel_vec     the vector of resels
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
resel_vec = zeros(1,4);

LKCs = LKCs(LKCs > 0);
if size(LKCs,1) > 1
    LKCs = LKCs';
end

D = length(LKCs);

scaling_vec = repmat(sqrt(4*log(2)), 1, D).^(1:D);

resel_vec(2:(D+1)) = LKCs./scaling_vec;

if nargin >= 2
    resel_vec(1) = L_0;
else
    resel_vec(1) = 1;
end

end

