function ss = smooth1(P,Q,s,dtype)
% SMOOTH1 
%--------------------------------------------------------------------------
% ARGUMENTS
% P         image(s) to be smoothed (or 3D array).
% Q         filename for smoothed image (or 3D array).
% s         [sx sy sz] Gaussian filter width {FWHM} in mm (or edges).
% dtype     datatype [Default: 0 == same datatype as P].
%--------------------------------------------------------------------------
% OUTPUT
% ss        Sums of squares of kernel.
%--------------------------------------------------------------------------
% SEE ALSO
% spm_type, spm_fileparts, spm_slice_vol, spm_matrix (used with 1 arg)

if ischar(Q) && isstruct(P)
    [pth,nam,ext,num] = spm_fileparts(Q);
    Q         = P;
    Q.fname   = fullfile(pth,[nam,ext]);
    if ~isempty(num)
        Q.n   = str2num(num);
    end
    if ~isfield(Q,'descrip'), Q.descrip = sprintf('SPM compatible'); end
    Q.descrip = sprintf('%s - conv(%g,%g,%g)',Q.descrip, s);

    if dtype~=0 % Need to figure out some rescaling
        Q.dt(1) = dtype;
        if ~isfinite(spm_type(Q.dt(1),'maxval')),
            Q.pinfo = [1 0 0]'; % float or double, so scalefactor of 1
        else
            % Need to determine the range of intensities
            if isfinite(spm_type(P.dt(1),'maxval')),
                % Integer types have a defined maximum value
                maxv = spm_type(P.dt(1),'maxval')*P.pinfo(1) + P.pinfo(2);
            else
                % Need to find the max and min values in original image
                mx = -Inf;
                mn =  Inf;
                for pl=1:P.dim(3)
                    tmp = spm_slice_vol(P,spm_matrix([0 0 pl]),P.dim(1:2),0);
                    tmp = tmp(isfinite(tmp));
                    mx  = max(max(tmp),mx);
                    mn  = min(min(tmp),mn);
                end
                maxv = max(mx,-mn);
            end
            sf      = maxv/spm_type(Q.dt(1),'maxval');
            Q.pinfo = [sf 0 0]';
        end
    end
end
