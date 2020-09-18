function obj = wfield( varmask, fibersize, field_type, field_params, xvals )
% WFIELD( varmask, fibersize, xvals ) constructs a Fields object having
% white noise in the fiber.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  varmask     Possible values are
%              - a 1 x D vector containing the size of the mask. The field
%                'mask' is then set to true( varmask ).
%              - a T_1 x ... x T_D logical array containing the mask.
% Optional
%  fibersize   a vector containing the size of the fiber. Default is 0,
%              i.e. the field is scalar.
%  field_type  a string giving the type of field. 'N': normal, 'T': tfield
%              'L': Laplacian field. Default is 'N' i.e. a white Gaussian
%              field.
%  field_params   if field_type is not 'N', then if it is 'T' field_params
%                 is the degrees of freedom of the t-statistic, if it is 'L' 
%                 then field_params is the scale of the Laplacian
%  xvals       a 1 x D cell array containing the xvals.
%              Default { 1:masksize(1), ..., 1:masksize(D) }.
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class Fields representing white noise, which is not
%      masked. 
%--------------------------------------------------------------------------
% EXAMPLES
% %% % Simple example with whole domain mask
% %% scalar field
% lat_data = wfield( [4 2 3] )
%
% %% many subjects field
% lat_data = wfield( [4 2 3], 100 )
%
% %% Simple example with mask
% mask = true( [ 4, 12 ] )
% mask = logical( pad_vals( mask ) )
% lat_data = wfield( mask, 1 )
% figure, clf,
% imagesc( lat_data ), colorbar
% title( 'not masked field' )
% % Generate masked data
% lat_data = Mask( wfield( mask, 1 ) )
% figure, clf,
% imagesc( lat_data )
% title( 'masked field' )
%
% % Degree 3 t-field
% wfield( [5,5], 10, 'T', 3)
%
% % Scale 1 Laplacian field
% wfield( [5,5], 10, 'L', 1)
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
if ~exist( 'fibersize', 'var')
    fibersize = 1;
end

if ~exist( 'field_type', 'var')
    field_type = 'normal';
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

% @Fabian this loop is unecessary right?
obj = Field( varmask );
if strcmp(field_type, 'normal') || strcmp(field_type, 'N')
    obj.field = randn( [ obj.masksize(1:obj.D) fibersize ] );
elseif strcmp(field_type, 't') || strcmp(field_type, 'T')
    obj.field = trnd( field_params, [ obj.masksize(1:obj.D) fibersize ] );    
elseif strcmp(field_type, 'l') || strcmp(field_type, 'L')
    obj.field = rlap( field_params, [ obj.masksize(1:obj.D) fibersize ] );
elseif strcmp(field_type, 'skew') || strcmp(field_type, 'S') || strcmp(field_type, 's')
    obj.field = (randn( [ obj.masksize(1:obj.D) fibersize ] ).^2-1)/sqrt(2);
elseif strcmp(field_type, 's2') || strcmp(field_type, 'S2')
    obj.field = (randn( [ obj.masksize(1:obj.D) fibersize ] ).^2-1)/sqrt(2) - (randn( [ obj.masksize(1:obj.D) fibersize ] ).^2-1)/sqrt(2);
elseif strcmp(field_type, 'uniform') || strcmp(field_type, 'U') || strcmp(field_type, 'u')
    obj.field = rand( [ obj.masksize(1:obj.D) fibersize ] ) - 1/2;
elseif strcmp(field_type, 'p') || strcmp(field_type, 'pearson') || strcmp(field_type, 'P')
    obj.field = pearsrnd(0,1,1,4, obj.masksize(1:obj.D),fibersize);
else
    error('This field type has not been implemented')
end

if exist( 'xvals', 'var' )
    obj.xvals = xvals;
end

return