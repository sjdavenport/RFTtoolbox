function obj = wnfield( varmask, fibersize, xvals )
% wnfield( masksize, fibersize, mask ) constructs a Fields object having
% white noise in the fiber.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  varmask     Possible values are
%              - a 1 x D vector containing the size of the mask. The field
%                'mask' is then set to true( varmask ). 
%              - a T_1 x ... x T_D logical array containing the mask.
%
% Optional
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
%
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------

% Check input argument and make it a vector
if  ~isnumeric( varmask ) && ~islogical( varmask ) 
    error( "varmask must be a 1xD or Dx1 numerical vector or a logical array." )
else
    svarmask = size( varmask );
    if isnumeric( varmask ) && ~( length( svarmask ) == 2 ...
                                                && any( svarmask == 1 ) )
        error( "The masksize must be a 1xD or Dx1 numerical vector." )
    elseif isnumeric( varmask ) && svarmask( 2 ) == 1
        varmask = varmask';
    end
end

%% Check optional input
%--------------------------------------------------------------------------

% Fill the defaults of the optional parameters if necessary
if ~exist( 'fibersize', 'var')
    fibersize = 1;
end

% Check the optional inputs
if  ~isnumeric( fibersize )
    error( "The fibersize must be a 1xK or Kx1 numerical vector." )
else
    tmp = size( fibersize );
    if ~( length( tmp ) == 2 && any( tmp == 1 ) )
        error( "The fibersize must be a 1xK or Kx1 numerical vector." )
    elseif tmp( 2 ) == 1
        fibersize = fibersize';
    end
end

%% Main function
%--------------------------------------------------------------------------

if islogical( varmask)
    obj = Field( varmask );
    obj.field = randn( [ size( varmask ) fibersize ] );    
else
    obj = Field( varmask );
    obj.field = randn( [ varmask fibersize ] );
end

if exist( 'xvals', 'var' )
    obj.xvals = xvals;
end

return