function [tcout, xvals_vecs] = convfield_t( lat_data, FWHM, resadd)
% CONVFIELD_T( lat_data, FWHM, spacing) computes a convolution t field
% with specified FWHM at points given by spacing
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      a Dim by nsubj array of data
% FWHM          the FWHM of the kernel with which to do smoothing
% resAdd        the amount of voxels added equidistantly inbetween the
%               existing voxels. Default is 1.
%--------------------------------------------------------------------------
% OUTPUT
% tcout         the convolution tfield
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
s_lat_data = size(lat_data);
Dim = s_lat_data(1:end-1);

% Obtain the number of dimensions
D = length(s_lat_data) - 1;

%%  main function
%--------------------------------------------------------------------------
% Obtain the convolution fields for each subject
[setofconvfields, xvals_vecs] = convfield_struct( lat_data, FWHM, resadd, D );

% Calculate the t-statistic
if D > 1
    tcout = mvtstat(setofconvfields, spacep(Dim,resadd));
else
    tcout = mvtstat(setofconvfields);
end

end
% 
% Dim_lat = size(lat_data);
% D = length(Dim_lat) - 1;
% nsubj = Dim_lat(end);
% 
% setofconvfields = zeros(Dim_lat);
% 
% if D == 1
%     for I = 1:nsubj
%         setofconvfields( :, I ) = convfield( lat_data, FWHM, spacing, D );
%     end
% elseif D == 2
%     for I = 1:nsubj
%         setofconvfields( :, :, I ) = convfield( lat_data, FWHM, spacing, D );
%     end
% elseif D == 3
%     for I = 1:nsubj
%         setofconvfields( :, :, :, I ) = convfield( lat_data, FWHM, spacing, D );
%     end
% end


