function [ Delta, Lambda, Omega ] = pointderivcov( fields, point )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % default option of opt1
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
field_data = fields.field;
N = size(field_data,2);

spacing = fields.xvals{1}(2) - fields.xvals{1}(1);
deriv = diff(field_data,1)/spacing;
deriv2 = diff(deriv,1)/spacing;

deriv = deriv - mean(deriv,2);

index = fields.xvals{1} == point;
Lambda = mean(deriv(index,:).*deriv(index, :), 2)*(N-1)/N;
Delta = mean(deriv(index,:).*deriv2(index, :), 2)*(N-1)/N;
Omega = mean(deriv2(index,:).*deriv2(index, :), 2)*(N-1)/N;

end

