function [tcfield, cfields] = convfield_t( lat_data, params )
% CONVFIELD_T( lat_data, FWHM, resadd ) computes a convolution t field
% with specified FWHM with resadd additional voxels added between the
% original points of the lattice.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  a Dim by nsubj array of data
%  params    an object of class ConvFieldParams.
% Optional
%
%--------------------------------------------------------------------------
% OUTPUT
% tcfield  an object of class Field representing the convolution tfield
% cfields  an object of class Field representing the convolution fields of
%          the lat_data
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
D = lat_data.D;

%%  main function
%--------------------------------------------------------------------------
% Obtain the convolution fields for each subject
cfields  = convfield( lat_data, params, 0 );

% Generate an empty Field with the given mask
tcfield = Field( cfields.mask );

% Calculate the t-statistic
if D > 1
    tcfield.field = mvtstat( cfields.field, spacep( Dim, params.resadd ) + 2* params.enlarge );
else
    tcfield.field = mvtstat( cfields.field );
end

tcfield.xvals = cfields.xvals;
tcfield = Mask(tcfield);
end

