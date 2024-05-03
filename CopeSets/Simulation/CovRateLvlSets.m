function [coveringRate] = CovRateLvlSets( truef, hatf, thresh, c, Bdry_test )

% Compute covering rate of level sets
% Input:
%  truef:  underlying true function on an D-dimensional array, it is saved as
%          an D-dimensional array, where the values correspond to the heights
%  hatf:   array of size [size(truef), N] containing N estimates of truef
%  thresh: array of size [size(truef), nquantiles, N, 2] corresponding to the thresholds for the
%          hatf, [ , , 1] are lower bounds, [ , , 2] are upper bounds
%  c:      targeted levelset
%  Bdry_test: if 0, applies old boundary test comparing just binarized sets.
%             if 1, also use linear interpolation boundary test. 
%Output:
% coveringRate is computed from all available hatf and computes the
% frequency that the true c-levelset of truef is completly contained in the
% thresholded hatfs.
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (fabian.telschow@hu-berlin.de), Alex Bowring
% Last changes: 10/25/2018
%__________________________________________________________________________

%%%%%% Compute variables of the field/simulation
sf    = size( truef );
D     = ndims( truef );
N     = size( hatf, D+1 );
index = repmat( {':'}, 1, D+1 );

truefm = reshape( repmat(reshape( truef, [prod(sf) 1] ), 1, N ) , [sf N] ) ;

% Finding the edges and weights for interpolating true field = c
if ( Bdry_test )
    bdry_params = getBdryparams(truef, c);
    bdry_length = bdry_params.length;
    low_thresh_bdry_values  = zeros(bdry_length, N);
    high_thresh_bdry_values = zeros(bdry_length, N);    
    hatf_bdry_values        = zeros(bdry_length, N);   
end

nquantiles = size(thresh,D+1);
coveringRate = zeros(1,nquantiles);
index2 = repmat( {':'}, 1, D );
index3 = repmat( {':'}, 1, 2 );

for i = 1:nquantiles
    thresh_tmp = squeeze(thresh(index2{:}, i, index3{:}));
    
    % ensure that for only one Sample the thresh_tmp is of the same form
    if N == 1
        thresh_tmp = reshape( thresh_tmp, [sf, 1, 2] );
    end

    if ( Bdry_test )
        % finding values of estimated field at locations where true
        % field = c 
        for n = 1:N          
            low_thresh_bdry_values(:,n)  = getBdryvalues(squeeze(thresh_tmp(index2{:},n,1)), bdry_params);
            high_thresh_bdry_values(:,n) = getBdryvalues(squeeze(thresh_tmp(index2{:},n,2)), bdry_params);
            hatf_bdry_values(:,n)        = getBdryvalues(hatf(index2{:},n), bdry_params);
        end
        violations = sum( ...
            any(reshape(...
                ( truefm < c & (hatf >= squeeze(thresh_tmp(index{:},2)))) |...
                ( truefm >= c & (hatf < squeeze(thresh_tmp(index{:},1))))...
            ,[prod(sf) N]), 1) | ...
            any(...% here might be a problem!
                hatf_bdry_values >= high_thresh_bdry_values | ...
                hatf_bdry_values < low_thresh_bdry_values...
            , 1 )) / N;
    else 
       violations = sum(any(reshape(...
        ( truefm < c & (hatf >= squeeze(thresh_tmp(index{:},2))) ) |...
        ( truefm >= c & (hatf < squeeze(thresh_tmp(index{:},1))))...
       ,[prod(sf) N]), 1 )) / N ; 
    end
    coveringRate(i) = 1-violations; 
end
