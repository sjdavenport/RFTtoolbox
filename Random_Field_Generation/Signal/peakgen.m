function Sig = peakgen( Mag, Rad, Smo, Dim, centre_locs )
% peakgen( Mag, Rad, Smo, Dim, centre_locs ) generates signal with peaks 
% that is defined everywhere (via convolution fields). WORK IN PROGRESS!!
%--------------------------------------------------------------------------
% ARGUMENTS
% Mag       is a vector of length npeaks giving the magnitude of the 
%           resulting signal at each peak. Note that if Mag is just a real
%           number rather than a vector then all peaks are taken to have
%           magnitude Mag. npeaks is then determined by the length of
%           centre_locs.
% Rad       is a vector of length npeaks giving the radius of the spheroid 
%           signal at each peak prior to being smoothed. Note that if Rad 
%           is just a real number rather than a vector then all peaks are 
%           taken to have radius Rad.
% Smo       is a vector of length npeaks that gives the smoothing applied
%           to the signal at each peak. The smoothing applied is Gaussian 
%           with the same FWHM in each x,y and z directions. Note that if Smo 
%           is just a real number rather than a vector then all peaks are
%           smoothed with FWHM Smo.
% Dim       a vector of length D giving the dimensions of the output image. 
%           For example Dim =[20,50] means that the output image is 20 x 50.
% centre_locs   is a cell array of length npeaks such that the nth entry is
%               a length D vector giving the coordinates of the centre location
%               of the nth peak.
%--------------------------------------------------------------------------
% OUTPUT
% Sig      an image of dimension Dim with peaks centred at the locations
%           specified by centre_locs.
%--------------------------------------------------------------------------
% EXAMPLES
% % 1D signal (NEED TO IMPLEMENT FOR 1D)!
% Sig = peakgen( 1, 3, 6, 100 )
% 
% % 2D signal
% Sig = peakgen([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
% surf(Sig)%
%
% % 2D Single peak plus noise
% Sig = peakgen(1, 10, 8, [90,90]); nsubj = 20;
% surf(Sig)
% lat_data = randn([90,90,nsubj]) + Sig;
% surf(mean(lat_data,3))
%--------------------------------------------------------------------------

%Set to sum(100*clock) to ensure that this is different each time.
randn('seed',sum(100*clock));   %-Random number generator initializaiton

%DEFAULT VARS

if nargin < 1
    Mag = 2;
end
if nargin < 2
    Rad = 10;
end
if nargin < 3
    Smo = 6;
end
if nargin < 4
    Dim = [91, 109, 91];
end
if nargin < 5
    centre_locs = {Dim/2 + 1/2};
end
    
npeaks = length(centre_locs);
if length(Mag) == 1
    Mag = repmat(Mag, 1, npeaks);
elseif length(Mag) ~= npeaks
    error('The number of peaks in Mag is not the same as in centre_peaks')
end

if length(Rad) == 1
    Rad = repmat(Rad, 1, npeaks);
elseif length(Rad) ~= npeaks
    error('The number of peaks in Rad is not the same as in centre_peaks')
end

if length(Smo) == 1
    Smo = repmat(Smo, 1, npeaks);
elseif length(Smo) ~= npeaks
    error('The number of peaks in Smo is not the same as in centre_peaks')
end

rimFWHM = 1.7;
%-----------Initialization of Some Variables
nDim    = length(Dim);

boundary2add = ceil(rimFWHM*max(Smo(:)))*ones(1,nDim);
wDim    = Dim + 2*boundary2add;  % Working image dimension

for peak = 1:npeaks
    centre_locs{peak} = centre_locs{peak} + boundary2add; % Need to do this to ensure the centre locs are in the right place even after boundary stuff.
end

%Below describes the actual bits of the image other than the extra
%truncation stuff that we have added on!
Trunc_x = {(ceil(rimFWHM*max(Smo(:)))+1):(ceil(rimFWHM*max(Smo(:)))+Dim(1))};
Trunc_y = {(ceil(rimFWHM*max(Smo(:)))+1):(ceil(rimFWHM*max(Smo(:)))+Dim(2))};

if nDim==2
    %Concatenates Trunc_x and Trunc_y into one array. Why is this
    %necessary? 
    TrnInd = cat(2, Trunc_x, Trunc_y); %Note cat(2,A,B) == [A,B]
else
    Trunc_z = {(ceil(rimFWHM*max(Smo(:)))+1):(ceil(rimFWHM*max(Smo(:)))+Dim(3))};
    TrnInd  = cat(2, Trunc_x, Trunc_y, Trunc_z);
end

Sig = zeros(wDim);
for peak = 1:npeaks
    Sig = Sig + SpheroidSignal(wDim, Rad(peak), Mag(peak), Smo(peak), centre_locs{peak}); %- Signal Should smooth here really!!
end

if nDim == 2
    Sig = Sig(TrnInd{1},TrnInd{2});
elseif nDim == 3
    Sig = Sig(TrnInd{1},TrnInd{2}, TrnInd{3});
end

end

