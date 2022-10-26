%% 2D ramp
dim = [100 100]; 
mu = repmat(linspace(1, 3), dim(2), 1);
nsubj = 240;
FWHM = 3;
c = 2;
lvls = 0.90;
nsim = 1000;

thresh    = zeros([dim, 1, nsim, 2]);
quantiles = zeros([1 nsim]);
hatmu     = zeros([dim, nsim]);
hatsigma  = zeros([dim, nsim]);
len_bdry  = zeros([1 nsim]);

for m = 1:nsim
    m
    noise = noisegen( dim, nsubj, FWHM );
    data = noise + mu;
    data = Field(data, true(dim)); % convert to field
    [ thresh_tmp, quantiles_tmp, hatmu_tmp, hatsigma_tmp, len_bdry_tmp ] = ...
        CopeSets( data, c, lvls );
    
    thresh(:,:,1,m,:) = thresh_tmp;
    quantiles(m)    = quantiles_tmp;
    hatmu(:,:,m)    = hatmu_tmp.field;
    hatsigma(:,:,m) = hatsigma_tmp.field;
    len_bdry(m) = len_bdry_tmp;
end

%%
% only voxels
coveringRateNaive = CovRateLvlSets( mu, hatmu, thresh, c, 0 )

% linear interpolated voxels
coveringRateInterpol = CovRateLvlSets( mu, hatmu, thresh, c, 1 )