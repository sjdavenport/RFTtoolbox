%% Introduction to Fields class
% Author: Fabian Telschow
%% Section 1: Introduction to Fields
% This tutorial provides a short introduction to the Fields class from the.
% RFTtoolbox. Please install the <http://github.com/sjdavenport/RFTtoolbox 
% RFTtoolbox> in order to run this correctly.
% See <http://sjdavenport.github.io> for more details.

% First you'll need to add the RFTtoolbox to your path
% by changing the following line:
folder_where_RFTtoolbox_is_saved = 'C:/Users/12SDa/davenpor/davenpor/Toolboxes/RFTtoolbox';
addpath( genpath( folder_where_RFTtoolbox_is_saved ) )

%% Basics on classes
% A class in Matlab is essentially a structure, which additionally has 
% methods solely working on them. The entries of the list are called
% properties.
% The main benefit is that a project can be structured better by coding
% important object as classes. Moreover, classes allow to overload basic 
% operators, which makes the code more readable.
%
% The basic object for a class is the Constructor, which constructs an
% object of the class. The BasicConstructor in matlab is ClassName().
MyField = Field()
% Note that all properties are set to a Default value. It is not always
% possible to construct an empty class as for example with the SepKernel
% class from this package, where the develeoper does not allowed for an
% empty object:
SepKernel()
%%
% A Field object has 3 important properties:
%
% - the field property an array of size T_1 x ... x T_D x F_1 x ... x F_K
% representing the actual field
%
% - the mask property an logical T_1 x ... x T_D array indicating which
% voxels belong to the domain of the field
%
% - xvals property indicating the voxelsizes in each dimension
%
% Since empty fields are not interesting, we can fill the properties as if
% it is a Matlab structure. 
MyField.mask    = true( [10 10 3] );
MyField.xvals   = {1:10, (2:11)/2, 2:4};
A = randn( [ 10 10 3 8 7] );
MyField.field = A;
MyField
%%
% Note that the field/mask and xvals need to be compatible, e.g. the
% following causes an error, since the xvals and field property have
% different dimensions than the mask:
MyField.mask    = true( [10 10] )
%% Other Basic Constructors
% The code above Field can be constructed in one line by recognising that
% the Fields() method allows for several different input variations, e.g.
MyField2 = Field( A, true( [10 10 3] ), {1:10, (2:11)/2, 2:4} )
MyField3 = Field( A, true( [10 10 3] ) )
MyField3.xvals = {1:10, (2:11)/2, 2:4}
% the function isequal checks whether all properties of a Field are the
% same:
isequal( MyField, MyField2 )
isequal( MyField, MyField3 )
%% Subreferencing a Field object
% A field object can be subreferenced as a usual array, which shows the
% power of overloading methods. Note that you either need to index the mask
% dimensions only or all dimensions.
MyField_sub = MyField(:, 1:3, 2)
% Note that the dimension D gets adapted automatically if you slice out a
% dimension
%% Plotting 1D and 2D Fields
% Also the functions plot and imagesc are overloaded such that you can use
% them directly on scalar fields, i.e. fiberD = 1 and D=1/2.
plot( MyField( 1, :, 1, 1, 1 ) )
imagesc( MyField( 1, :, :, 1, 1 ) )
title("Dimensions get scaled by the xvals values")
%% Special Fields
% Two special fields can be constructed similarily easily, i.e. a white noise
% field we basically constructed before
mask = true( [10 10 3] );
mask(1:5,:,1) = false;
fibersize   = 1;
xvals       = {1:10, (2:11)/2, 2:4};
wn = wfield( mask, fibersize, 'N', 0, xvals );
imagesc( wn(:,1,:) )
% A constant fiber field, i.e. each mask value contains the same array A
A  = ones( [ 2 2 3 ] )
cf = constfield( A, mask, xvals )
imagesc( cf(:,1,:, 2,2,2) ), colorbar
% Note that in most cases xvals and mask are optional and if not provided
% are set to standard values.
%% Masking a field
% If you looked carefully for wn we changed the mask to be not everywhere
% true and the masked-property changed its value to 0. this is because the
% field property is not equal to all 0, -Inf or NaN. In order to apply the
% mask to your field you can simply use the Mask() method and you will see
% that the masked property has changed:
wn_masked = Mask( wn )
figure, clf,
subplot(1,2,1)
imagesc( wn(:,:,1) )
title("not masked field")
colorbar
subplot(1,2,2)
imagesc( wn_masked(:,:,1) )
title("masked field")
colorbar
% In principle the Mask( field, val, mask ) function allows for an
% arbitrary value val as input and a different mask, than stored in the
% Field object. By default val=0,
%% Algebra with Fields
% Fields cannot only represent scalar fields but also vector fields or
% other objects from fiber bundles over the mask. Therefore basic algebra
% operations especially for matrix multiplication, transposing etc. are
% available for fields overloading the basic operators. You can play around
% with some of the examples.
%
% Vector scalar multiplication
v = constfield( [1 2 3], mask )
tmp = 3*v
figure,clf,
subplot(1,3,1)
imagesc( tmp(:,1,:,1,1) ), colorbar
subplot(1,3,2)
imagesc( tmp(:,1,:,1,2) ), colorbar
subplot(1,3,3)
imagesc( tmp(:,1,:,1,3) ), colorbar
% Vector field scalar field multiplication
v = constfield( [1 2 3], mask )
w = wfield( mask )
tmp = w .* v
figure, clf
imagesc( tmp(:,1,:,1,1) ), colorbar
% Matrix multiplication I: Vector field times Vector field
v = constfield( [1 2 3]', mask )
w = wnfield( mask, 3 )
tmp = v.' * w
figure, clf
imagesc( tmp(:,1,:,1) ), colorbar
% Matrix multiplication II: Matrix times Vector field
v = constfield( [1 2 3]', mask )
A = wfield( mask, [3 3] )
tmp = A * v
figure, clf
imagesc( tmp(:,1,:,1) ), colorbar
% Matrix multiplication III: Matrix times Matrix field
A = constfield( ones( [ 4, 3 ] ), mask )
B = wfield( mask, [3 3] )
tmp = A * B
figure, clf
imagesc( tmp(:,1,:,1,1) ), colorbar
% Note that there are many other possible operations which are implmented
% or might be implemented at some point. If something urgent is missing
% please contact one of the developers.
%% Collapse
% Another sometimes useful function is collapse, which collapses either the
% domain or the fiber to a one dimensional space.
A = wfield( mask, [3 3] )
collapse( A )
collapse(A, 'fiber')
%% ConvField Class
% This package also includes a ConvField class, which is a subclass of the
% Fields class. A subclass contains all properties of the parent class and
% all methods implemented for the parent class can be used on the subclass
% object. However, special methods and properties can be introduced for the
% subclass. In this toolbox the ConvField class describes a field contained
% from smoothing with a kernel. It is generated as follows:

% Get lattice data
wn = wfield( mask, 1 );
% Input data
% Kernel: SepKernel object. Numeric creates a seperable isotropic Gaussian
% kernel with FWHM kernel
kernel_FWHM = 3;
% Derivative type: 0/1/2 supported giving the actual field or its first or
% second derivative
derivtype = 0; 
% resolution increase of the field, i.e. number of voxels added inbetween
% each voxel
resadd = 3;
% lat_mask: mask the lat_data before smoothing or not
lat_masked = true;
% enlarge the mask by voxres in high resolution. 
enlarge = ceil( resadd / 2 );
% Get the params object for a convfield
params = ConvFieldParams(kernel_FWHM*ones([1 3]),...
                         resadd,...
                         enlarge,...
                         true );
% Construct ConvField object.
cfield = convfield( wn, params, derivtype )
imagesc( cfield(:,:,3) )

% Note that a ConvField basically tracks the Kernel, resadd and enlarge
% addtionally to the other field properties
%