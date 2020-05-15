function [tcout, xvals_vecs] = convfield_t( lat_data, FWHM, spacing)
% CONVFIELD_T( lat_data, FWHM, spacing) computes a convolution t field
% with specified FWHM at points given by spacing
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      a Dim by nsubj array of data
% FWHM          the FWHM of the kernel with which to do smoothing
% spacing       the interval at which to compute the convolution field.
%               I.e. if size(lat_data) = [10,20] and spacing = 0.1 the field 
%               will be evaluated at the points 1:0.1:10
%--------------------------------------------------------------------------
% OUTPUT
% tcout         the convolution tfield
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
D = length(size(lat_data)) - 1;

[setofconvfields, xvals_vecs] = convfield( lat_data, FWHM, spacing, D );

tcout = mvtstat(setofconvfields);

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


