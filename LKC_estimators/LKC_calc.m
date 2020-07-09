function [ L, L0 ] = LKC_calc( lat_data, Kernel, resadd, lat_masked, version )
% LKC_calc( lat_data, Kernel, resadd, lat_masked, version ) calculates the
% LKCS.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  an object of class Field with fiberD == 1.
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
%  resadd     the amount of voxels added equidistantly in between the
%             existing voxels. Default is 1.
%  lat_masked a logical, if true lat_data is masked by the provided mask in
%             the lat_data Field object. Default is true.
% Optional
%  version      not included yet
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
cfield  = convfield_Field( lat_data, Kernel, 0, resadd, lat_masked, enlarge );
dcfield = convfield_Field( lat_data, Kernel, 1, resadd, lat_masked, enlarge );
    
[L,L0] = LKC_voxmfd_est( cfield, dcfield );

end
