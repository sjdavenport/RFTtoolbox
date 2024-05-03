function out = Gratio( x, mu1, mu2, Sigma1, Sigma2 )
% Gratio calculates the ratio between two normal pdfs.
%--------------------------------------------------------------------------
% ARGUMENTS
% 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport

% num = exp(-(1/2)*(x-mu1)'*inv(Sigma1)*(x-mu1));
% denom = exp(-(1/2)*(x-mu2)'*inv(Sigma2)*(x-mu2));
% 
% out = num./denom

out = mvnpdf(x, mu1,Sigma1)./mvnpdf(x, mu2,Sigma2);

% normpdf(x, mu1,Sigma1)/normpdf(x, mu2,Sigma2)

end

