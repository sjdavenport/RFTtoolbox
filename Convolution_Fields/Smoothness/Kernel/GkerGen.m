function val = GkerGen( x, Sigma )
% GKERGEN( x, Sigma ) calculates the Gaussian Kernel
%--------------------------------------------------------------------------
% ARGUMENTS
% x        a D by nevals matrix where each column is a D-dimensional 
%          vector at which to evaluate the kernel.
% Sigma    the covariance matrix
%--------------------------------------------------------------------------
% OUTPUT
% val      a number giving the value of the kernel
%--------------------------------------------------------------------------
% EXAMPLES
% x = [0,0]'; GkerGen( x )
%
% x = [0, 1, -1, 2; 0, 1, -1, 2];
% val = GkerGen(x);
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
D = size(x, 1);

if nargin < 2
    Sigma = eye(D);
end

Sigmainv = inv(Sigma);
detSigmahalf = det(Sigma)^(1/2);
val = zeros([1,size(x,2)]);
for I = 1:size(x,2)
    val(I) = exp(-x(:,I)'*Sigmainv*x(:,I)/2)/(sqrt(2*pi)^D)/detSigmahalf;
end

end

