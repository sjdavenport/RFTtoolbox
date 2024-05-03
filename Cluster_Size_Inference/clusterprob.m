function prob = clusterprob( k, FWHM, u, D )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% k     the cluster extent
% FWHM    the smoothness
% u       the cluster defining threshold
% D       the number of dimensions
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% k = 150; FWHM = 5; u = 2.5; D = 3;
% clusterprob( k, FWHM, u, D )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
if D == 1
    omega = 2;
elseif D == 2
    omega = pi;
elseif D == 3
    omega = (4/3)*pi;
end

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
Lambdahalf = ((4*log(2))^(D/2))/FWHM^D;
prob = exp(-(k*u^D*2^(-D/2)*omega^(-1)*Lambdahalf)^(2/D));

end

