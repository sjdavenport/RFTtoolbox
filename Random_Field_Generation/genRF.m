function rfs = genRF( nreal, df, FWHM, Dim, asvector )
% genRF( nreal, df, FWHM, Dim, asvector ) returns an nreal by prod(Dim) 
% set of random fields with degrees of freedom df (allowing generation of 
% Gaussian, t and F fields) which have a given FWHM. 
%--------------------------------------------------------------------------
% ARGUMENTS
% nreal   the number of realizations
% df    If df == 1, a GRF is generated. If df = [1, n] a t field with
%       n degrees of freedom is generated. If df = [m, n] an F field 
%       with m and n degrees of freedom is generated.
% FWHM  The smoothness.
% Dim   The dimensions of each field. Default is Dim = [91,109,91].
%       asvector determines whether an nreal by prod(Dim) or a nreal by Dim 
%       array is returned. Default takes as_vector = 1.
% asvector  0/1 if 1 then an nreal by prod(Dim) array is generated. If 0
%           then the array is nreal by Dim.
%--------------------------------------------------------------------------
% OUTPUT
% rfs   An nreal by prod(Dim) array where each row is a random field with
%       the specified degrees of freedom.
%--------------------------------------------------------------------------
% EXAMPLES
% rfs = genRF(10, [1,5], 6) %Generates 10 t-fields each with 5 degrees of
%                           %freedom and generated from FWHM 6 fields.
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport.
if nargin < 4
    Dim = [91,109,91];
end
if nargin < 5
    asvector = 1;
end
if asvector ~= 1
    error('Need to check through the t generation as haven''t set an option for this yet')
end

if asvector == 1
    rfs = zeros(nreal, prod(Dim));
else
    rfs = zeros(nreal, Dim);
end

if isequal(df, 1)
    if asvector == 1
        rfs = noisegen( Dim, nreal, FWHM, 3 );
    else
        rfs = noisegen( Dim, nreal, FWHM, 1 );
    end
else
    if length(df) ~= 2 
        error('For df not equal to 1 you need a vector of length 2 to specify the degrees of freedom in the numerator and the denominator.')
    end
end

if df(1) == 1 && length(df) > 1
    for real = 1:nreal
        grfs = noisegen( Dim, (df(2)+1), FWHM, 3 );
        [~,~,~,rfs(real, :)] = meanmos(grfs); %Generates the t-statistic.
    end
elseif df(1) > 1
    error('Haven''t designed generation of F fields yet, ie where the df in the numerator is not 1!')
    for real = 1:nreal
        grfs = noisegen( Dim, (df(2)+1), FWHM, 3 );
        [~,~,~,rfs(real, :)] = meanmos(grfs); %Generates the t-statistic.
    end
end

end

