function xvals_hr = xvals_highres( xvals, resadd, enlarge )
%  xvals_highres( xvals, resadd, enlarge ) computes a high resolution
% version of a xvals.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   xvals   an 1 x D cell array containing vectors
%   resadd  the amount of equidistant voxels introduced inbetween the
%           voxels  
% Optional
%   enlarge numeric denoting the amount of voxels added by dilating the
%           high resolution mask. Default ceil(resadd/2), which if resadd
%           is odd means that the boundary voxels are exactly half the
%           distance between the voxels shifted with respect to the
%           original mask voxels.
%--------------------------------------------------------------------------
% OUTPUT
%   xvals_hr  a 1 x D cell array containing the resolution
%             increased version of the xvals vector 
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
% -------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow
%--------------------------------------------------------------------------


%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

% Check whether the mask is logical
if ~iscell( xvals )
    error( "xvals must be a 1 x D cell!" );
end

% Check whether resadd is numeric
if ~isnumeric( resadd )
    error( "resadd must be a positive natural number!" );
else
    % Ensure that resadd is a positive integer
    resadd = ceil( abs( resadd ) );
end


%% Add/check optional input
%--------------------------------------------------------------------------
if ~exist( 'enlarge', 'var' )
   % Default option of 'enlarge'
    enlarge = ceil( resadd / 2 );
end

if ~isnumeric( enlarge )
    error( "enlarge must be a postive integer!" )
else
    if enlarge ~= ceil( enlarge ) || enlarge < 0
       error( "enlarge must be a postive integer!" )
    end
end

%% Main function
%--------------------------------------------------------------------------
D = length( xvals );

% Get the difference between voxels with resolution increase
dx = NaN * ones( [ 1 D ] );
for d = 1:D
    dx(d) = xvals{d}(2)-xvals{d}(1);
end
dx_hr = dx ./ ( resadd + 1 );

% Get the xvals_hr cell
xvals_hr = cell( [ 1, D ] );
for d = 1:D
    xvals_hr{d} = ( xvals{d}(1) - enlarge * dx_hr(d) ) : dx_hr(d) : ...
                                    ( xvals{d}(end) + enlarge * dx_hr(d) );
end

return