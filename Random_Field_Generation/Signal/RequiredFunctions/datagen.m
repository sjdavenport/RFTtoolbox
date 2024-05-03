function [ data ] = datagen( Dim, nSubj, FWHM, Mag, Rad )
% DATAGEN( Dim, nSubj, Smo, Mag, Rad, rimFWHM ) generates 2d smoothed
% images with signal and correlated gaussian noise that comes from smoothed
% N(0,1) random variables.
%--------------------------------------------------------------------------
% ARGUMENTS
% nSubj     is the number of subjects which corresponds to the number of
%           separate images to generate.
% Dim       Gives different options that correspond to the dimensions.
%           5: 2D noise. 6: 3D noise.
% Rad       A positive number that is the equatorial radius of Spheroid.
% Mag       the magnitude of the signal.
% Smo       a real number that equals FWHM_x = FWHM_y = FWHM_z ie we
%           assume that the smoothing is equal in all directions. This can
%           also be a vector in which case it will specify different levels
%           of smoothing.
% rimFWHM   the width of padding around the image to be truncated to reduce
%           the edge effect.
%--------------------------------------------------------------------------
% OUTPUT
% data      an array of size Dim(1) by Dim(2) by nSubj where the third
%           index runs over the number of subject: for each subject giving
%           an image, points of which are identified by the first two
%           indices.
%--------------------------------------------------------------------------
% EXAMPLES
% data = datagen(2);
% surf(mean(data, 3));
%
% data = datagen(3, 20);
% mean_est = mean(data, 4);
% surf(mean_est(:,:,50));
% surf(mean_est(:,:,30));
%
% 2D random noise:
% data = datagen(5);
% surf(mean(data,3))
%
% % 3D noise with signal
% data = datagen(stdsize, 20, 6);
% m_data = mean(data, 4);
% surf(m_data(:,:,50))
%--------------------------------------------------------------------------
% SEE ALSO
% spm_conv, spm_smooth, SpheroidSignal,

%Set to sum(100*clock) to ensure that this is different each time.
randn('seed',sum(100*clock));   %-Random number generator initializaiton

%DEFAULT VARS
reshape = 0;
if (nargin < 1)
    Dim = [100,100];
end
if (Dim == 2)
    Dim = [100,100];
end
if (Dim == 3)
    Dim = [91, 109, 91];
end
if (Dim == 4)
    reshape = 1;
    Dim = [91, 109, 91];
end
use_sig = 1;
if (Dim == 5)
    %Generate random noise images. 2D.
    use_sig = 0;
    Dim = [100, 100];
end
if (Dim == 6)
    %Generate random noise images. 3D.
    use_sig = 0;
    Dim = [91, 109, 91];
end

if (nargin < 2)
    nSubj  = 20;  % Number of subjects
end
if (nargin < 3)
    FWHM = 6;
end
if (nargin < 4)
    Mag = 2;
end
if (nargin < 5)
    Rad = 10;
end

data = noisegen( Dim, nSubj, FWHM, reshape); %Generate noise.

if use_sig == 1
    Sig = gensig( Dim, 6, Mag, Rad ); %Note we have preset the signal. It may have a different smootheness to the noise!
    if reshape
        Sig = Sig(:);
        for subj = 1:nSubj
            data(:,subj) = data(:,subj) + Sig;
        end
    else
        for subj = 1:nSubj
            data(:,:,:,subj) = data(:,:,:,subj) + Sig;
        end
    end
end

