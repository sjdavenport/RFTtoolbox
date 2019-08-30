function data = noisegen( Dim, nSubj, FWHM, shape_of_array )
% NOISEGEN( Dim, nSubj, FWHM, shape_of_array ) generates an array of N(0,1) 
% noise smoothes this with a gaussian kernel with a certain FWHM and scales 
% so that the resulting field has variance 1. By default this array is Dim 
% by nSubj, but this can be changed using the shape_of_array parameter.
%--------------------------------------------------------------------------
% ARGUMENTS
% Dim       The dimensions of the image.
% nSubj     The number of subjects.
% FWHM      The FWHM of the kernel
% shape_of_array   0/1/2, default is 0. Determines the shape of the array.
%           If 0 data is Dim by nSubj,
%           If 1 data is nSubj by Dim,
%           If 2 data is prod(Dim) by nSubj instead of
%           nSubj by Dim. 
%           If 3 data is nSubj by prod(Dim) instead of
%           Dim by nSubj.
%--------------------------------------------------------------------------
% OUTPUT
% data      an array of size Dim by nSubj where the first
%           index runs over the number of subjects: for each subject giving
%           an image, points of which are identified by the rest of the
%           indices. (Note that the shape of the array can be changed using
%           the shape of array parameter.)
%--------------------------------------------------------------------------
% EXAMPLES
% noise = noisegen(160,20,6);
% plot(mean(noise,2), 'linewidth', 2)
%
% Dim = [100,100]
% noise = noisegen(Dim, 20, 6);
% noise_mean = mean(noise,3);
% surf(noise_mean)
%
% global stdsize
% noise = noisegen(stdsize, 20, 6, 1);
% noise_mean = mean(noise,1);
%--------------------------------------------------------------------------
% Authors:
% Sam Davenport and Thomas Nichols
if nargin < 2
    nSubj  = 20; 
end
if nargin < 3
    FWHM = 6;
end
if nargin < 4
    shape_of_array = 0;
end

nDim    = length(Dim); %The number of dimensions.

rimFWHM = 1.7; %The number of standard deviations or something like that in either direction.
wDim    = Dim + 2*ceil(rimFWHM*FWHM)*ones(1,nDim);  %The increased dimension, needed to deal with the edge effect.

%The Trunc variables describe the subset that corresponds to the image 
%rather than the surrounding voxels that we have added on.
Trunc_x = {(ceil(rimFWHM*FWHM)+1):(ceil(rimFWHM*FWHM)+Dim(1))}; %The size of the image in the x direction.

if nDim == 1
    TrnInd = cat(2, Trunc_x);
elseif nDim==2
    %Concatenates Trunc_x and Trunc_y into one array. Why is this necessary?
    Trunc_y = {(ceil(rimFWHM*FWHM)+1):(ceil(rimFWHM*FWHM)+Dim(2))}; %The size of the image in the y direction.
    TrnInd = cat(2, Trunc_x, Trunc_y); %Note cat(2,A,B) == [A,B]
else
    Trunc_y = {(ceil(rimFWHM*FWHM)+1):(ceil(rimFWHM*FWHM)+Dim(2))}; %The size of the image in the y direction. 
    Trunc_z = {(ceil(rimFWHM*FWHM)+1):(ceil(rimFWHM*FWHM)+Dim(3))}; %The size of the image in the z direction.
    TrnInd  = cat(2, Trunc_x, Trunc_y, Trunc_z);
end

%Initializes an array to store the noise
RawNoise = zeros([wDim nSubj]);

%Initialize Data Matrix. This allows for some special cases.
if shape_of_array == 2
    data = zeros( prod(Dim), nSubj );
elseif shape_of_array == 3
    data = zeros(  nSubj, prod(Dim) );
elseif shape_of_array == 1
    data   = zeros( [nSubj Dim] );
else
    data   = zeros( [Dim nSubj] );
end

%Loop over subjects and generate noise for each one.
for subj = 1:nSubj
    %Set the noise to be iid N(0,1).
    if nDim == 1
        RawNoise(:, subj)       = randn([1,wDim]);
    elseif nDim == 2
        RawNoise(:,:,subj)      = randn(wDim);
    else
        RawNoise(:,:,:, subj)   = randn(wDim);
    end
    
    % Smooth noise
    if FWHM == 0 %Ie if no smoothing is to be applied.
        if (nDim==1)
            Noises    = RawNoise(:,subj);
        elseif (nDim==2)
            Noises    = RawNoise(:,:,subj);
        else
            Noises    = RawNoise(:,:,:,subj);
        end
    else
        Noises    = zeros(wDim);
        %Smooths the noise and divides by tt which is the sum of squares of
        %the kernel. This ensures that noise is unit variance.
        if nDim == 1
            [Noises,tt] = spm_conv_mod(RawNoise(:,subj),FWHM,FWHM);
        elseif nDim == 2
            [Noises,tt] = spm_conv_mod(RawNoise(:,:,subj),FWHM,FWHM);
        else
            tt       = spm_smooth_mod(RawNoise(:,:,:,subj),Noises,FWHM);
        end
        Noises    = Noises/sqrt(tt); %Done to standardize.
    end
  
    % Truncate to avoid edge effects
    if nDim == 1 
        trunNoises    = Noises(TrnInd{1});
    elseif nDim == 2
        trunNoises    = Noises(TrnInd{1},TrnInd{2});
    else
        trunNoises    = Noises(TrnInd{1},TrnInd{2},TrnInd{3});
    end
    
    % Store smoothed truncated signal image in data array.
    if shape_of_array == 2
        data(:,subj) = trunNoises(:);
    elseif shape_of_array == 3
        data(subj,:) = trunNoises(:);
    elseif shape_of_array == 0
        if nDim == 1
            data(:,subj) = trunNoises;
        elseif nDim == 2
            data(:,:,subj) = trunNoises;
        else
            data(:,:,:,subj) = trunNoises;
        end
    elseif shape_of_array == 1
        if nDim == 1
            data(subj,:) = trunNoises;
        elseif nDim == 2
            data(subj,:,:) = trunNoises;
        else
            data(subj,:,:,:) = trunNoises;
        end
    end
end

end

