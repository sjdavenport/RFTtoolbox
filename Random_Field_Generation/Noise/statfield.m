function stat_field = statfield( Dim, nsubj, params, truncmult )
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
% % 1D 
% FWHM = 3; resadd = 3; params = ConvFieldParams(FWHM, resadd);
% Dim = 10; nsubj = 20; f = statnoise( Dim, nsubj, params )
%
% % 2D
% FWHM = [3,3]; resadd = 3; params = ConvFieldParams(FWHM, resadd);
% Dim = [10,10]; nsubj = 20; f = statnoise( Dim, nsubj, params )
% % Note this has size 37x37 as (resadd+1)*10 - resadd = 37, i.e. 4
% % associated with all except the final voxel point!
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
params.enlarge = 0;
if params.enlarge ~= 0
    error('The params enlarge parameter must be set to 0')
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% params = ConvFieldParams(FWHM, resadd, 0);

params.kernel.truncation = truncmult*params.kernel.truncation;
cutoff = params.kernel.truncation(1);

lat_data = wfield(Dim + 2*cutoff, nsubj);
smooth_field = convfield(lat_data, params);
bigger_masksize = smooth_field.masksize;

number2cut = cutoff*(params.resadd+1);
if smooth_field.D == 1
%     stat_field = cut_field(smooth_field, cutoff);
    stat_field = smooth_field( (number2cut+1):bigger_masksize(1) - number2cut);
elseif smooth_field.D == 2
    stat_field = smooth_field( (number2cut+1):bigger_masksize(1) - number2cut,...
            (number2cut+1):bigger_masksize(2) - number2cut);
elseif smooth_field.D == 3
    stat_field = smooth_field( (number2cut+1):bigger_masksize(1) - number2cut,...
            (number2cut+1):bigger_masksize(2) - number2cut, ...
            (number2cut+1):bigger_masksize(3) - number2cut);
end
stat_field.kernel = smooth_field.kernel;
for I = 1:smooth_field.D
    stat_field.xvals{I} = stat_field.xvals{I} - cutoff;
end
% stat_field = cut_field(smooth_field, cutoff);

end

