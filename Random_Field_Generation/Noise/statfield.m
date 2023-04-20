function [stat_field, deriv_field, deriv2_field] = statfield( Dim, nsubj, params, truncmult, do_derivs )
% STATFIELD generates a stationary smooth noise field. (ATM only 1D.)
% For resadd greater than 1 this is approximately stationary for FWHM >= 3.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%     truncmult     default set to 1
%     do_derivs     0/1 determining whether to evaluate the derivative and 
%                   second derivative of the field as well as the field
%                   itself. Default is 1 i.e. to do so! (Though this has
%                   only been implemented in 1D.)
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% % 1D 
% FWHM = 3; resadd = 3; params = ConvFieldParams(FWHM, resadd);
% Dim = 10; nsubj = 20; f = statfield( Dim, nsubj, params )
%
% % 2D
% FWHM = [3,3]; resadd = 3; params = ConvFieldParams(FWHM, resadd);
% Dim = [10,10]; nsubj = 20; f = statfield( Dim, nsubj, params )
% % Note this has size 37x37 as (resadd+1)*10 - resadd = 37, i.e. 4
% % associated with all except the final voxel point!
% 
% stat_field = statnoise( 100, 10, 3, 11 )
%
% global PIloc; load([PIloc, 'Variance/storevars']); FWHM = 20;
% params = ConvFieldParams(FWHM, 0)
% std_est = sqrt(allvars(FWHM))*sqrt(100/99);
% stat_field = statfield( 100, 10000, params)*(1/std_est);
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

if ~exist('do_derivs', 'var')
    do_derivs = 1;
end

% if ~exist('setunitvariance', 'var')
%     setunitvariance = 0;
% end

% Allow for default FWHM input
D = 1;
if isnumeric(params)
    FWHM = params
    params = ConvFieldParams( repmat(params,1,D), 0 );
else
    FWHM = 0;
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

% Set these to NaN so that in higher dimensions they are not generated!
deriv_field = NaN;
deriv2_field = NaN;
if smooth_field.D == 1
%     stat_field = cut_field(smooth_field, cutoff);
    stat_field = smooth_field( (number2cut+1):bigger_masksize(1) - number2cut);
    if do_derivs 
        smooth_field_derivs = convfield(lat_data, params, 1);
        smooth_field_derivs2 = convfield(lat_data, params, 2); 
        deriv_field = smooth_field_derivs( (number2cut+1):bigger_masksize(1) - number2cut);
        deriv2_field = smooth_field_derivs2( (number2cut+1):bigger_masksize(1) - number2cut);
    end
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

% if setunitvariance
%     if FWHM
%        warning('To set unit variance you need the file storecfieldvar.m')
%         global PIloc
%         load([PIloc,'Variance/storevars'], 'allvars')
%         stat_field = (stat_field./sqrt(allvars(FWHM)))*(100/99);
%     else
%         error("FWHM wasn't numeric so can't scale the variance")
%     end
% end

end

