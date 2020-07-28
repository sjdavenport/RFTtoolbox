function EEC = EEC_spm( thresholds, L, L0, field_type, nsubj, use_old )
% EEC_spm( LKCs, field_type ) calculates the expected Euler characteristic 
% of the excursion set above the given thresholds
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  L    the postive LKCs
%  L0   the 0th LKC i.e. the Euler characteristic of the 
%  field_type  
%--------------------------------------------------------------------------
% OUTPUT
%  EEC   the expected Euler characteristic of the excursion set above the
%        given threshold
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = 6; resadd = 1; nsubj = 50; nvox = 100;
% lat_data = wnfield(nvox, nsubj);
% cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
% dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
% [L,L0] = LKC_voxmfd_est( cfield, dcfield );
% EEC_spm( [2,3], L, L0, 'T', nsubj )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

% Set the default number of subjects (i.e. if you're using a Gaussian
% field)
if ~exist('nsubj', 'var')
    nsubj = 1;
end

% Allow former use of SPM for comparison (to be removed!)
if ~exist('use_old', 'var')
    use_old = 0;
end


% Obtain the number of dimensions
D = length(L);

%%  Main Function Loop
%--------------------------------------------------------------------------
% Uses SPM to calculate the expected Euler Characteristic (note the SPM
% functions are very confusing and we need to replace this but good to
% check against)
resel_vec = LKC2resel( L, L0 );

if strcmp(field_type, 'Z') 
    if ~isequal(nsubj, 1)
        error('For Gaussian fields need df = 1')
    end
    df = 1;
elseif strcmp(field_type, 'T')
    df = [1,nsubj-1];
else
    error('fields other than Z and T are not coded')
end

% Calulate the expected Euler characteristic above each threshold using SPM
EEC = zeros(1, length(thresholds));
for I = 1:length(EEC)
    ECdensities  = spm_ECdensity(field_type,thresholds(I),df); 
    EEC(I) = resel_vec*ECdensities(1:D+1);
end
%Note that we could easily parallelize the above loop as you just need to
%calculate the EC densities to get vectors of length(EEC) for each 
% and then these vectors while scaling by the LKCs. (But it's fast anyways
% so doesn't really matter).

if use_old
    for I = 1:length(EEC)
        [~,~,EEC(I)] = spm_P_RF(1,0,thresholds(I),df,field_type,resel_vec,1);
    end
end

end

