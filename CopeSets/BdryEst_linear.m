function bdryValues = BdryEst_linear( field, c, delta )
% BdryEst_erodedilate( field, c, delta, outfield ) estimates the boundary
% of the excursion set by eroding and delating the mask to get a small
% neighbourhood of the boundary.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   field  
%   c         threshold value for excursion
%   delta     
%--------------------------------------------------------------------------
% OUTPUT
%  field   
%
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR:  Fabian Telschow
%--------------------------------------------------------------------------
%% Check mandatory input and get important constants
%--------------------------------------------------------------------------


%% Add/check optional values
%--------------------------------------------------------------------------
switch nargin
    case 2
        delta     = 666;
end


%% Main function  
%--------------------------------------------------------------------------
% Get useful constants
sf = field.fieldsize;
N  = sf( end );
index = repmat( {':'}, 1, field.D );

% Compute the estimated signal if no function provided
if sum( delta(:) ) == 666
    delta = mean( field );
    delta = delta.field;
end

% ensure there are values above the threshold!
if sum( delta(:) >= c ) == 0
    error( "The threshold c is to high. There are now values exceeding it!" )
end

% Get the values for the linearly interpolated boundary
bdry_param  = getBdryparams( delta, c );
bdryValues = zeros( bdry_param.length, N );

% Get the linearly interpolated boundary values of the field
for n = 1:N
    bdryValues( :, n ) = getBdryvalues( field.field( index{:}, n ),...
                                        bdry_param );
end

end