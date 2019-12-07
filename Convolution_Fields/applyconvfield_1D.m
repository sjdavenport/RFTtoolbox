function field = applyconvfield(tval, xvalues_at_voxels, Kernel, data_mean)
%varx = @(x) applyconvfield(x, xvalues_at_voxels, FWHM, ones(1,100));
% varx_deriv = @(x) applyconvfield(x, xvalues_at_voxels, @(x)Gkerderiv(x, FWHM), ones(1,100));
%(varx(10+h) - varx(10))/h
% varx_deriv(10)

if isnumeric(Kernel)
    if Kernel < 1
        warning('Are you sure the FWHM and increm have been written the right way around?')
    end
    Kernel = @(x) Gker(x,Kernel);
end

if size(tval, 2) == 1
    tval = tval';
end
if size(xvalues_at_voxels, 1) == 1
    % Need xvalues_at_voxels to be a column vector
    xvalues_at_voxels = xvalues_at_voxels';
end
L = length(xvalues_at_voxels);
eval_mate = repmat(tval,L, 1) - repmat(xvalues_at_voxels,1, length(tval)); %repmat(tval, matrix(rep(tval, each = L), L) - rep(xvalues_at_voxels, length(tval));
field = data_mean*Kernel(eval_mate);
end