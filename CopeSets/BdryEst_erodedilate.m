function [ field, mask ] = BdryEst_erodedilate( field, c, delta, outfield )
% BdryEst_erodedilate( field, c, delta, outfield ) estimates the boundary
% of the excursion set by eroding and dilating the mask to get a small
% neighbourhood of the boundary.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   field  
%   c         threshold value for excursion
%   delta     
%   outfield  boolean. If true the output field is an Object of class
%             field. Default 1.
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
        outfield  = 1;
    case 3
        outfield  = 1;
end


%% Main function  
%--------------------------------------------------------------------------
% Get the dimension of the field
sf = field.fieldsize;

% Compute the excursion set estimated from the provided signal
if sum( delta(:) ) == 666
    % Compute the excursion set estimated from the sample mean
    mfield = mean( field );
    A_c = ( mfield.field >= c);
    
else
    % Compute the excursion set estimated from the provided signal
    A_c = ( delta  >= c);   
end

% ensure there are values above the threshold!
if sum( A_c ) == 0
    error( "The threshold c is to high. There are now values exceeding it!" )
end

mask = dilate_mask( A_c, 1 ) & dilate_mask( ~A_c, 1 );

field.mask = mask;

if ~outfield
    field = reshape( field.field( ...
                        repmat( mask, [ ones( [ 1, field.D ] ), sf( end ) ] )...
                        ),...
                     [ sum( mask(:) ) sf( end ) ] );
end

return
