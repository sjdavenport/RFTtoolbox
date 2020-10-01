function stat_field = statnoise( Dim, nsubj, FWHM, resadd )
% STATNOISE generates a stationary smooth noise field. (ATM only 1D.)
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% stat_field = statnoise( 100, 10, 3, 11 )
%
% global PIloc; load([PIloc, 'Variance/storevars']); FWHM = 20;
% stat_field = statnoise( 100, 50, FWHM, 11 )*(1/sqrt(allvars(FWHM)));
% plot(var(stat_field.field,0,2));
%
% nea = wfield(100,50); FWHM = 20; 
% params = ConvFieldParams(20, 11, 0);
% neac = convfield(nea, params)*(1/sqrt(allvars(FWHM)));
% plot(var(neac.field,0,2));
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
params = ConvFieldParams(FWHM, resadd, 0);
cutoff = ceil(4*FWHM2sigma(FWHM));

lat_data = wfield(Dim(1)+2*cutoff, nsubj);
smooth_field = convfield(lat_data, params);

stat_field = cut_field(smooth_field, cutoff);

end

