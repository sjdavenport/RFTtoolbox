%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the Field class object and its functions
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%% Basic constructor
%% %% Generation of empty field and filling the properties: nargin = 1
F = Field()
% Fill the mask field
F.mask = true( [ 3 4 5 ] )
% Fill the fields field. Note that the first dimensions have to agree with
% the dimensions of the mask
F.field = randn( [ 3 4 5 1000] )
% Fill the xvals field. Note that it must be compatible with mask
% dimensions
F.xvals = { 1:3, [ 0.5, 3, 6, 6.5 ], 15:19 }
%% %% nargin = 1
%% % Input can be a vector defining the size of the mask.
% mask field is defaulted to be true( varargin{1} ) and xvals is
% defaulted by a sequence. Hence only a field need to be added.
F = Field( [ 4 5 ] )
F.field = randn( [ 4 5 10 18 1000] )
% mask or xvals can be changed, if needed. Note that masked property will
% change to 0, since the field is nonzero/NaN/-Inf outside the mask.
mask = true( [ 2, 3 ] )
mask = logical( pad_vals( mask ) )
F.mask = mask
figure, clf
imagesc( F.field( :, :, 2 ) ), colorbar
% To mask the field simply apply the Mask function
F = Mask(F)
figure, clf
imagesc( F.field( :, :, 2 ) ), colorbar
%% % Input can be a 1xD cell array containing vectors.
% mask field is defaulted to be of compatible dimension with the xvals.
% Hence only a field need to be added.
F = Field({ 1:3, [ 0.5, 3, 6, 6.5 ], 15:19 })
% Fill field
F.field = randn( [ F.masksize 100 ] )

%% %% nargin = 2
%% % Input can be a field and a mask.
% xvals is defaulted, yet can be changed
mask = true( [ 4, 12 ] )
mask = logical( pad_vals( mask ) )

F = Field( randn( size( mask ) ), mask )
%% % Input can be a field and or xvals.
% mask is defaulted, yet can be changed.
F = Field( randn( [ 3 4 5 19 ] ), { 1:3, [ 0.5, 3, 6, 6.5 ] } )
%% % Input can be a field and a dimension of the mask.
% mask and xvals are defaulted, yet can be changed.
F = Field( randn( [ 3 4 5 19 ] ), 1 )
F = Field( randn( [ 3 4 5 19 ] ), 3 )
%% %% nargin = 3
%% % Input can specify a field, a mask and xvals.
F = Field( randn( [ 3 4 5 19 ] ), true( [ 3 1 ] ), { [ 0.1, 0.3, 0.7 ] } )