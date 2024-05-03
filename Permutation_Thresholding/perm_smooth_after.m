function [ threshold, vec_of_maxima, permuted_tstat_store ] = ... 
    perm_smooth_after( data, mask, FWHM, alpha, nperm, show_loader, store_perms )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% data = randn(91,109,10);
% MNImask = imgload('MNImask');
% MNImask2D = MNImask(:,:,50);
% [ threshold, vec_of_maxima ] = perm_smooth_after( data, MNImask2D>0, FWHM )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% s_data = size(data);
dim = size(mask);
D = length(dim);
data_vectorized = vec_data(data, mask);

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'nperm', 'var' )
   % Default value
   nperm = 1000;
end

if ~exist( 'show_loader', 'var' )
   % Default value
   show_loader = 0;
end

if ~exist( 'alpha', 'var' )
   % Default value
   alpha = 0.05;
end

if ~exist('store_perms', 'var')
    store_perms = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
nsubj = size(data_vectorized,2);

mask = logical(mask);

tstat = unwrap(nan2zero(mvtstat(data_vectorized)), mask);
[tstat_smooth, ss] = fconv(tstat, FWHM, D);
tstat_smooth = tstat_smooth.*mask./sqrt(ss);
    
vec_of_maxima = zeros(1,nperm);
vec_of_maxima(1) = max(tstat_smooth(:));

% Compute bernoulli random variables for the sign flipping
random_berns = 2*(binornd(1,0.5, nsubj, nperm )-1/2);

if store_perms == 1
    permuted_tstat_store = zeros([sum(mask(:)), nperm]);
    permuted_tstat_store(:,1) = tstat_smooth(mask);
else
    permuted_tstat_store = NaN;
end

for I = 2:nperm
    if show_loader == 1
        loader(I-1, nperm-1, 'perm progress:');
    end
    
    random_berns_for_iter = random_berns(:, I);
    random_sample_negative = find(random_berns_for_iter < 0);

    data_perm = data_vectorized;
    data_perm(:,random_sample_negative) = -data_vectorized(:,random_sample_negative);
    
    tstat_perm = unwrap(nan2zero(mvtstat(data_perm)), mask);
    [tstat_smooth_perm, ss] = fconv(tstat_perm, FWHM, D);
    tstat_smooth_perm = tstat_smooth_perm.*mask./sqrt(ss);
    
    if store_perms == 1
        permuted_tstat_store(:,I) = tstat_smooth_perm(mask);
    end
    
%     vec_of_maxima(I) = max(tstat_smooth_perm(:));
    vec_of_maxima(I) = max(tstat_smooth_perm(mask));
end

threshold = prctile(vec_of_maxima, 100*(1-alpha) );

end

