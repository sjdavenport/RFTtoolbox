function stat_field = statnoise( Dim, nsubj, params, truncmult )
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
% params = ConvFieldParams(
%
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
if ~exist('truncmult', 'var')
    truncmult = 1;
end
    
%%  Add/check optional values
%--------------------------------------------------------------------------
if params.enlarge ~= 0
    error('The params enlarge parameter must be set to 0')
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% params = ConvFieldParams(FWHM, resadd, 0);
params.kernel.truncation = truncmult*params.kernel.truncation;
cutoff = ceil(4*FWHM2sigma(FWHM));

lat_data = wfield(Dim(1)+2*cutoff, nsubj);
smooth_field = convfield(lat_data, params);
bigger_mask = smooth_field.mask;

number2cut = cutoff*(kernel.resadd+1);
if smooth_field.D == 1
    stat_field = smooth_field( (number2cut+1):bigger_mask(1) - number2cut);
elseif smooth_field.D == 2
    stat_field = smooth_field( (number2cut+1):bigger_mask(1) - number2cut,...
            (number2cut+1):bigger_mask(2) - number2cut);
elseif smooth_field.D == 3
    stat_field = smooth_field( (number2cut+1):bigger_mask(1) - number2cut,...
            (number2cut+1):bigger_mask(2) - number2cut, ...
            (number2cut+1):bigger_mask(3) - number2cut);
end
% stat_field = cut_field(smooth_field, cutoff);

end

