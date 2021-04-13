function obj = cnfield( varmask, FWHMcor, voxmap, FWHM, fibersize, xvals )
% wnfield( masksize, fibersize, mask ) constructs a Fields object having
% white noise in the fiber.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  varmask     Possible values are
%              - a 1 x D vector containing the size of the mask. The field
%                'mask' is then set to true( varmask ). 
%              - a T_1 x ... x T_D logical array containing the mask.
%  FWHMcor     a 1 x D vector containing the FWHM applied to correlate the
%              white noise
%  voxmap      a size varmask array representing a bijection of the domain
%
% Optional
%  FWHM        a 1 x D vector containing the FWHM applied after applying
%              the voxelmap
%  fibersize   a vector containing the size of the fiber. Default is 0,
%              i.e. the field is scalar.
%  xvals       a 1 x D cell array containing the xvals.
%              Default { 1:masksize(1), ..., 1:masksize(D) }.
%
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class Fields representing white noise, which is not
%      masked. 
%
%--------------------------------------------------------------------------
% EXAMPLES
% %% % Simple example with whole domain mask
% %% scalar field
% lat_data = wnfield( [4 2 3] )
%
% %% many subjects field
% lat_data = wnfield( [4 2 3], 100 )
%
% %% Simple example with mask
% mask = true( [ 4, 12 ] )
% mask = logical( pad_vals( mask ) )
% lat_data = wnfield( mask, 1 )
% figure, clf,
% imagesc( lat_data ), colorbar
% title( 'not masked field' )
% % Generate masked data
% lat_data = Mask( wnfield( mask, 1 ) )
% figure, clf,
% imagesc( lat_data )
% title( 'masked field' )
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------

% Check input argument and make it a vector
if  ~isnumeric( varmask ) && ~islogical( varmask ) 
    error( 'varmask must be a 1xD or Dx1 numerical vector or a logical array.' )
else
    svarmask = size( varmask );
    if isnumeric( varmask ) && ~( length( svarmask ) == 2 ...
                                                && any( svarmask == 1 ) )
        error( 'The masksize must be a 1xD or Dx1 numerical vector.' )
    elseif isnumeric( varmask ) && svarmask( 2 ) == 1
        varmask = varmask';
    end
end

%% Check optional input
%--------------------------------------------------------------------------

% Fill the defaults of the optional parameters if necessary
if ~exist( 'FWHM', 'var')
    FWHM = 0;
end

% Fill the defaults of the optional parameters if necessary
if ~exist( 'fibersize', 'var')
    fibersize = 1;
end

% Check the optional inputs
if  ~isnumeric( fibersize )
    error( 'The fibersize must be a 1xK or Kx1 numerical vector.' )
else
    tmp = size( fibersize );
    if ~( length( tmp ) == 2 && any( tmp == 1 ) )
        error( 'The fibersize must be a 1xK or Kx1 numerical vector.' )
    elseif tmp( 2 ) == 1
        fibersize = fibersize';
    end
end

%% Main function
%--------------------------------------------------------------------------

% Get a white noise field
obj = wnfield( varmask, fibersize );
if exist( 'xvals', 'var' )
    obj.xvals = xvals;
end

% Smooth the white noise field to introduce correlation
params = ConvFieldParams( FWHMcor * ones( [ 1 obj.D ] ), 0, 0, false );
obj  = convfield( obj, params );

% Generate voxelmap for subjects
N = prod( obj.masksize );
voxmap = voxmap(:)';
voxmap1 = voxmap;
for n = 2:prod( obj.fibersize )
    voxmap = [ voxmap, ( (n-1) * N ) + voxmap1 ];
end

% Subref the field
obj.field = reshape( obj.field( voxmap ),...
                     obj.fieldsize );
   
if FWHM ~= 0
    params = ConvFieldParams( FWHM * ones( [ 1 obj.D ] ), 0, 0, false );
    obj  = convfield( obj, params );
end

return