function [ indices, npeaks ] = lmindices( map, top, mask )
% LMINDICES( map, top, mask ) finds the indices of the top local maxima of
% a given map. It uses a connectivity criterion of 18.
%--------------------------------------------------------------------------
% ARGUMENTS
% map   a 3D map which could be vectorized (ie map = map(:)) but not
%       necessary I guess.
% top  specifies the number of the top local maxima to find if its
%       intially set as inf it then counts the number of peaks within the mask!
% mask  an image of 1s and 0s that masks the data. (Needs to be 91 by 109
%       by 91).
%--------------------------------------------------------------------------
% OUTPUT
% indices   a vector of indices of local maxima.
% npeaks    the number of peaks returned = length(indices) this is only
%           really relevant when top = Inf as then the number of peaks that
%           is returned is variable.
%--------------------------------------------------------------------------
% EXAMPLES
% Sig = imgload('oldmean');
% indices = lmindices( Sig(:), 1, ones(size(Sig)) )
% Sig(indices)
%
% oldmean = imgload('oldmean');
% top_lm_indices = lmindices(oldmean, 20)
%
% a = zeros([91,109,91]);
% a(16,100,40) = 5;
% max_index = lmindices(a, 1, 'all') %need to use 'all' here as (16,100,40)
% % doesn't lie within the MNI mask of the brain.
% a(max_index)
%
% a(10,50,35) = 3;
% top_2_lms = lmindices(a, 2, 'all');
% a(top_2_lms(2))
%
% To investigate:
% Why does the following give an odd answer:
% Sig = gensig(3, 6, 2, 10);
% lmindices(Sig,3) %It finds maxima that are not local maxima!
% %Ie: convind(340934) = [48,41,35]. But Sig(48,41,35) < Sig(48,41,36)!!
% %Something to do with the boundary I think.
%
% noise = noisegen(stdsize, 1, 30);
% noise(30,40,50) = 100;
% lmindices(noise, Inf, 'all')
%
% noise = noisegen(stdsize, 1, 30);
% lmindices(noise, Inf, 'all')
%--------------------------------------------------------------------------
%Put more error checking in this function!!
if nargin < 2
    top = 1;
end
if nargin < 3
    mask = 'all';
%     mask = imgload('MNImask');
end

global stdsize
if strcmp(mask, 'all')
    mask = ones(stdsize);
end

if ~isequal(numel(map), 902629)
    error('The map must be of size [91,109,91] or a vectorized version of this.')
end


if ~isequal(size(mask), stdsize)
    error('The mask must be of size [91,109,91]')
end

if any(mask(:))
    map = nan2zero(map);
    mask = nan2zero(mask);
    
    % if mask_out_nans
    %    map = nan2zero(reshape(map,stdsize).*mask);
    %    mask = ones(stdsize);
    % end
    % nanmask = reshape(isnan(map),stdsize).*mask; %Finds the nans that lie within your imposed mask. Need to make them zeros.
    % if any(nanmask(:))
    %     if mask_out_nans
    %         notnanmask = ~nanmask; %Gives 0s where there are nans within your imposed mask and 1s everywhere else.
    %         mask = mask.*notnanmask;
    %     else
    %         error('There are nans in your images which need to be masked out')
    %     end
    % end
    
    s = size(mask);
    [x,y,z] = ind2sub(s,find(mask));
    XYZ     = [x y z]';
    
    %Old incorrect version:
    % map = map.*mask;
    % [~,Z,M] = spm_max(map(:),XYZ);
    
    map = map(:);
    mask = mask(:);
    
    maskedmap = map(mask>0);
    
    [~,Z,M] = spm_max(maskedmap,XYZ);
    [~, sidx] = sort(Z, 'descend');
    %Is Z a vector of the lm values? Check this!
    
    if isinf(top)
        top = length(sidx);
    elseif top > length(sidx)
        top = length(sidx); %Return the number of local maxima found.
        warning('The number of desired maxima is less than the number within the mask')
    end
    
    topMlocs = M(:,sidx(1:top));
    indices = sub2ind(s,topMlocs(1,:),topMlocs(2,:),topMlocs(3,:));
    indices = indices(map(indices) ~= 0);
    %If you want a way arround this so that you're able to identify maxima that
    %are zero you could instead assign the NaNs in a linear way so that they
    %take values below the minimum of the image and are arranged so that
    %neighbouring ones are never location maxima.
    
    npeaks = length(indices);
else
    indices = NaN;
    npeaks = 0;
end

end

