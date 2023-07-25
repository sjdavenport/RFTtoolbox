function [tfced] = tfce(image,H,E,connectivity,dh)
% tfce(image,H,E,connectivity,dh) performs TFCE
%--------------------------------------------------------------------------
% ARGUMENTS
%  image: a 3D matlab array
%  H: height exponent (default is 2)
%  E: extent exponent (default is 0.5)
%  connectivity: connectivity used to compute the connected components
%  dh: size of steps for cluster formation
%--------------------------------------------------------------------------
% OUTPUT

%--------------------------------------------------------------------------
% EXAMPLES
% dim = [50,50]; nsubj = 50; FWHM = 2;
% Sig = 0.1*peakgen(1, 10, 8, dim);
% data = wfield(dim, nsubj);
% data.field = data.field + Sig;
% tstat = convfield_t(data, FWHM);
% tstat_tfce = tfce(tstat.field,2,0.5,8,0.05)
% subplot(1,2,1)
% surf(tstat.field)
% subplot(1,2,2)
% surf(tstat_tfce)
%--------------------------------------------------------------------------
% AUTHOR: Mark Allen Thornton and Samuel Davenport
%--------------------------------------------------------------------------
D = length(size(image));
if ~exist( 'connectivity', 'var' )
   % Default value
   if D == 2
       connectivity = 8;
   elseif D == 3
       connectivity = 26;
   end
end

if ~exist( 'H', 'var' )
   % Default value
   H = 2;
end

if ~exist( 'E', 'var' )
   % Default value
   E = 0.5;
end

if ~exist( 'dh', 'var' )
   % Default value
   dh = 0.1;
end

% set cluster thresholds
threshs = 0:dh:max(image(:));
threshs = threshs(2:end);
nthreshs = length(threshs);

% find positive voxels (greater than first threshold)
nvox = length(image(:));

% find connected components
vals = zeros(nvox,1);
cc = arrayfun(@(x) bwconncomp(bsxfun(@ge,image,x),connectivity), threshs);
for h = 1:nthreshs
    clustsize = zeros(nvox,1);
    ccc = cc(h);
    voxpercc = cellfun(@numel,ccc.PixelIdxList);
    for c = 1:ccc.NumObjects
        clustsize(ccc.PixelIdxList{c}) = voxpercc(c);
    end
    % calculate transform
    curvals = (clustsize.^E).*(threshs(h)^H);
    vals = vals + curvals;
end
tfced = NaN(size(image));
tfced(:) = vals.*dh;

end

