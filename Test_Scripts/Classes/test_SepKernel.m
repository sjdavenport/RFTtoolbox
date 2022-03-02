git%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests and illustrates the SepKernel class and
%%%     its functions
%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% Basic constructor
% Define function handles as kernels
f  = @(x)(x.^2 + 2);
f1 = @(x)(2 * x);
f2 = @(x) 2;

g  = @(x)(x.^3 + 3);
g1 = @(x)(3 * x);
g2 = @(x) 3;

% define truncations for a 3D SepKernel
truncation = [3 3 3];
dtruncation = [4 3 3];
d2truncation = [4 3 5];
truncMatrix = [truncation; dtruncation; d2truncation];

% Generate a SepKernel object with D-dimensional domain and fill
D = 3;
NewKernel = SepKernel(D)

% Fill all the entries for the kernel and the derivative kernel
kern = cell([1 D]);
for d = 1:D
    kern{d} = f;
end
NewKernel.kernel     = kern;
NewKernel.truncation = truncation;

for d = 1:D
    kern{d} = f1;
end
NewKernel.dkernel  = kern;
NewKernel.dtruncation = dtruncation;

for d = 1:D
    kern{d} = f1;
end
NewKernel.d2kernel = kern;
NewKernel.d2truncation = d2truncation;

NewKernel.adjust = [1 1 3];

% Corresponding fields are filled!
NewKernel

% The above is complicated. The SepKernel constructor can also be used
% in the way that you provide the informations piece by peice or
% simultaneously.
% First argument:  Dimension of the domain
% Second argument: Function handle or 1 x D cell containing function handles
%                  for the kernel
% Third argument:  Function handle or 1 x D cell containing function handles
%                  for the derivative of the kernel
% Fourth argument: Function handle or 1 x D cell containing function handles
%                  for the second derivative of the kernel
% Fifth argument:  numeric, 1 x D vector or 3 x D matrix containing the
%                  truncation for the kernels
% Sixth argument:  numeric or 1 x D vector containing the adjust value for
%                  the kernel
NewKernel2 = SepKernel(D, f, f1, f2, truncMatrix, [1 1 3])  % same result as above using

%% % Generating Gaussian kernels depending on FWHM
D = 2;
GaussKern = SepKernel(D, 5)
% Equivalently you can use GaussKernel:
GaussKern2 = GaussKernel(D, 5)

% Adding truncation
GaussKern = SepKernel(D, 5, [1 3])
% additionally adding the truncation and adjust

GaussKern = SepKernel(D, 5, [1 3], [2 1]) % truncation and adjust

%% % Non-isotropic Gaussian kernels
GaussKern = GaussKernel(D, 2);
% Equivalently you can use GaussKernel:
GaussKern2 = GaussKernel(D, 10);

non_isotropic_kernel = GaussKern;
non_isotropic_kernel.truncation(2) = GaussKern2.dtruncation(1);
non_isotropic_kernel.kernel{2} = GaussKern2.kernel{1};
non_isotropic_kernel.dtruncation(2) = GaussKern2.dtruncation(1);
non_isotropic_kernel.dkernel{2} = GaussKern2.dkernel{1};
non_isotropic_kernel.d2truncation(2) = GaussKern2.dtruncation(1);
non_isotropic_kernel.d2kernel{2} = GaussKern2.d2kernel{1};

lat_data = wfield([40,40]);
params = ConvFieldParams(non_isotropic_kernel, 5);
f = convfield(lat_data, params);
figure(1)
imagesc(f)

% Equal code as above, but simpler
non_isotropic_kernel = GaussKernel(2, [2 10]);
params = ConvFieldParams(non_isotropic_kernel, 5)
f = convfield(lat_data, params);
figure(2)
imagesc(f)


%% %-----------------------------------------------------------------------
% You can generate also Sep kernels with different kernels in different
% directions. Here f,g are 1D functions and f1,g1 its derivatives and f2,g2
% the second derivatives.
KernWithDifferentSmoothings = SepKernel(2, {f, g}, {f1, g1}, {f2, g2}, [2 3], [1 3])

%% %-----------------------------------------------------------------------
% Generating general kernels: in the Sep Kernel class it is possible to add
% the numerical derivatives instead of providing all. Just leave the
% derivatives as they are in the generation or provide @(x) NaN handles if
% you use mixtures.
D = 3;
KernWithoutDeriv = SepKernel(D, f, truncMatrix, [1 1 3])

% Fill the derivatives by numerical ones:
h = 0.01;
KernWithDerivatives  = NumericDerivatives( KernWithoutDeriv, h )
% or with default h
KernWithDerivatives2 = NumericDerivatives( KernWithoutDeriv )

% Test whether first and second derivatives are almost equal
x = 0.5;
trueDeriv = f1(0.5)
KernWithDerivatives.dkernel{1}(x)
KernWithDerivatives2.dkernel{1}(x)

trueDeriv2 = f2(0.5)
KernWithDerivatives.d2kernel{1}(x)
KernWithDerivatives2.d2kernel{1}(x)

% You can also fill in missing derivatives
Ex2 = SepKernel(2, {f, g}, {@(x) NaN, g1}, {f2, @(x) NaN}, [2 3], [1 3])
Ex2derivatives = NumericDerivatives( Ex2 )

