function zerolocs = varonederivzeros( lat_data, FWHM, std_field, initial_zero_locs )
% CONVDERIVZEROS( lat_data, Kernel, initial_zero_locs ) find the zeros of
% the derivative of a convolution field
%--------------------------------------------------------------------------
% ARGUMENTS
% 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% % 1D deriv root finding
% FWHM = 3;
% nvox = 10;
% lat_data = normrnd(0,1,1,nvox);
% maxestloc = lmindices(lat_data, 1)
% convderivzeros( lat_data, FWHM, maxestloc )
% findconvpeaks( lat_data, FWHM,[NaN,maxestloc] )
% ninter = 0.05;
% xvals_fine = 1:ninter:nvox;
% plot(xvals_fine, convfield(lat_data, FWHM, 0.05, 1))
%
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

D = size(initial_locs, 1);

field = @(x) applyconvfield( x, lat_data, FWHM );

% varonefield = @(x) field(x)./sqrt(var_field(x)); 
varonefield = @(x) field(x)./std_field(x); 

h = 0.00001;
derivfield = getderivs(varonefield, D, h);

zerolocs = fzero(derivfield, initial_zero_locs);

end

