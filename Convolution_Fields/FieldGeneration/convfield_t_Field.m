function [tcfield, cfields] = convfield_t_Field( lat_data, Kernel, resadd, enlarge )
% CONVFIELD_T( lat_data, FWHM, resadd ) computes a convolution t field
% with specified FWHM with resadd additional voxels added between the
% original points of the lattice.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  a Dim by nsubj array of data
%  Kernel    either an object of class SepKernel or a numeric.
%            If class SepKernel:
%              if derivtype = 0: the fields 'kernel' and 'truncation' must
%                                be specified.
%              if derivtype = 1: the fields 'dkernel' and 'dtruncation'
%                                must be specified.
%              if derivtype = 2: the fields 'd2kernel' and 'd2truncation'
%                                must be specified.
%
%            If Kernel is numeric, the convolution field is generated by 
%            smoothing with an isotropic Gaussian kernel with FWHM = Kernel.
%            Truncation and adjust_kernel are set to be default values.
% Optional
%  resadd        the amount of voxels added equidistantly inbetween the
%                existing voxels. Default is 1.
%  xvals_vecs    a D-dimensional cell array whose entries are vectors giving 
%               the xvalues at each each dimension along the lattice. It 
%               assumes a regular, rectangular lattice (though within a given
%               dimension the voxels can be spaced irregularly).
%               I.e suppose that your initial lattice grid is a
%               4by5 2D grid with 4 voxels in the x direction and 5 in
%               the y direction. And that the x-values take the values:
%               [1,2,3,4] and the y-values take the values: [0,2,4,6,8].
%               Then you would take xvals_vecs = {[1,2,3,4], [0,2,4,6,8]}.
%               The default is to assume that the spacing between the
%               voxels is 1. If only one xval_vec direction is set the
%               others are taken to range up from 1 with increment given by
%               the set direction.
%--------------------------------------------------------------------------
% OUTPUT
% tcfield         the convolution tfield
%--------------------------------------------------------------------------
% EXAMPLES
% %% 1D convolution t field
% nvox = 10; nsubj = 20; resadd = 20; FWHM = 2;
% lat_data = normrnd(0,1,[nvox,nsubj]);
% [tconvfield, xvals_vecs] = convfield_t(lat_data, FWHM, resadd);
% lattice_tfield = convfield_t(lat_data, FWHM, 0);
% plot(1:nvox, lattice_tfield, 'o-')
% hold on
% plot(xvals_vecs{1}, tconvfield)
% title('1D convolution t fields')
% legend('Convolution field', 'Lattice Evaluation')
% xlabel('voxels')
% ylabel('t field')
% 
% %% 2D convolution field

%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Get the dimensions of the data
Dim = lat_data.masksize;

% Obtain the number of dimensions
D = lat_data.D;

%%  add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'resadd', 'var' )
   % default option of resadd
   resadd = 1;
end

if ~exist( 'enlarge', 'var' )
   % default option of resadd
   enlarge = ceil(resadd/2);
end

%%  main function
%--------------------------------------------------------------------------
% Obtain the convolution fields for each subject
cfields  = convfield_Field( lat_data, Kernel, 0, resadd, 1, enlarge );

% Calculate the t-statistic
if D > 1
    tcfield = mvtstat(cfields.field, spacep(Dim,resadd)+2*enlarge);
else
    tcfield = mvtstat(cfields.field);
end

end
