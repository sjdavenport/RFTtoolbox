function EC = EulerCharCrit( f, D, mask, version )
% EulerCharCrit(f, D, mask, version)
% Computes the Euler characteristic curve of random fields represented as a
% stepfunction.
%--------------------------------------------------------------------------
% ARGUMENTS
%   f    an array K_1, x ... x K_D x N of random fields for which the Euler
%        characteristic curves should be computed
%   D    an integer containing the dimension of the domain of the random
%        fields
%   mask an boolean array K_1, x ... x K_D having TRUE, for the values
%        voxels belonging to the mask. Default: true(size(f, 1:D)).
%   version a string. If "C" the fast C implementation is used (default)
%           if "matlab" a slow matlab only implementation is used.
%           If the string is "CN" with N an integer, the image is split into
%           N parts and EC is computed parallelized using parfor.
%           Note that this will require the parallel computing toolbox. 
%--------------------------------------------------------------------------
% OUTPUT
%   EC curve is computed for all N fields f. The 2 x N_crit array contains
%   in the first column the heights of the critical values and in the
%   second column the Euler characterisitc of the upper level set.
%   4-connectivity for 2D fields and 6-connectivity for 3D fields is used.
%--------------------------------------------------------------------------
% EXAMPLES  
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%%%%%%%%% Preliminary constants, input check etc.
if D > 3
    error( "The method is unfortunately until now only implemeted for D <= 3." )
end

switch D
    case 1
        cc = 2;
    case 2
        cc = 4;
    case 3
        cc = 6;
end

% Get constants from the input field
sf = size( f );

% introduce index for accessing fields of different dimensions
index  = repmat( {':'}, 1, D );

% Initialize the EC curve cells, one for each realisation of the field
if length( sf ) > D
    nEC = sf( end );
    EC  = cell( [ 1 nEC ] );
else
    nEC = 1;
    EC  = cell( 1 );
end

% Set default parameter, if not provided
if ~exist( 'mask', 'var' )
    % masks which voxels should be considered
    mask = true( sf(1:D) );
    L0 = 1;
elseif ~all( size( mask ) == sf( 1:D ) )
    error( "Please specify an input mask, which has the same dimension as the domain of the data" )
else
    L0 = EulerChar( mask, 0.5, D );
end

if ~exist( 'version', 'var' )
    version = "C";
end

% prepare for parallelisation, if assumed.
if ~isempty( str2num( version ) )
    Npar    = str2num( version );
    version = "C";
    
    % save the state of the CPU's are open already
    state_gcp = isempty(gcp('nocreate'));

    % open connection to CPUs, if not already established
    if( state_gcp && Npar > 1 ) 
        parpool( Npar );
        state_gcp = 42;
    end
else
    Npar    = 0;
    state_gcp = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set all values outside the mask to -Inf
if ~all(mask(:))
    for kk = 1:size( f, D+1 )
        tmp = f( index{:}, kk );
        tmp( ~mask ) = -Inf;
        f( index{:}, kk ) = tmp;
    end
end

% prepare the field for the EulerCharCrit_c input
if Npar == 0
    % Pad large negative values to the array, which is required for 
    % EulerCharCrit_c input
    f_tmp = f;
    if length( sf ) == D
        f = -Inf * ones( sf + repmat( 2, [ 1 D ] ) );
    else
        f = -Inf * ones( sf + [ repmat( 2, [ 1 D ] ) 0 ] );
    end

    switch D
        case 1
            f( 2:end-1, : ) = f_tmp;
        case 2
            f( 2:end-1, 2:end-1, : ) = f_tmp;
        case 3
            f( 2:end-1, 2:end-1, 2:end-1, : ) = f_tmp;
    end
    clear f_tmp
else
    % find Npar break points
    par_breaks = 1:floor( sf(1) / Npar ):sf(1);
    if length( par_breaks ) == Npar+1
        par_breaks(end) = sf(1);
    else
        par_breaks = [ par_breaks, sf(1)];
    end
    
    f_tmp = cell( [ 1 Npar ] );
    index  = repmat( {':'}, 1, D );
    
    for k = 1:Npar
        if length( sf ) == D
            if k < Npar
                f_tmp{k} = -Inf * ones( [ par_breaks(k+1)-par_breaks(k) sf(2:end) ]...
                                + repmat( 2, [ 1 D ] ) );
            else
                f_tmp{k} = -Inf * ones( [ par_breaks(k+1)-par_breaks(k)+1 sf(2:end) ]...
                                + repmat( 2, [ 1 D ] ) );
            end
        else
            if k < Npar
                f_tmp{k} = -Inf * ones( [ par_breaks(k+1)-par_breaks(k) sf(2:end) ]...
                                + [ repmat( 2, [ 1 D ] ) 0 ] );
            else
                f_tmp{k} = -Inf * ones( [ par_breaks(k+1)-par_breaks(k)+1 sf(2:end) ]...
                                + [ repmat( 2, [ 1 D ] ) 0 ] );
            end
        end
        
        if k == 1
            switch D
                case 1
                    f_tmp{k}(2:end, :) = f( ( 1:par_breaks(k+1) ),...
                                              index{:} );
                case 2
                    f_tmp{k}(2:end, 2:end-1, :) = f( ( 1:par_breaks(k+1) ),...
                                                       index{:} );
                case 3
                    f_tmp{k}(2:end, 2:end-1, 2:end-1, :) = f( ( 1:par_breaks(k+1) ),...
                                                                index{:} );
            end
        elseif k == Npar
            switch D
                case 1
                    f_tmp{k}(1:end-1, :) = f( ( (par_breaks(k)-1):par_breaks(k+1) ),...
                                             index{:} );
                case 2
                    f_tmp{k}(1:end-1, 2:end-1, :) = f( ( (par_breaks(k)-1):par_breaks(k+1) ),...
                                             index{:} );
                case 3
                    f_tmp{k}(1:end-1, 2:end-1, 2:end-1, :) = f( ( (par_breaks(k)-1):par_breaks(k+1) ),...
                                             index{:} );
            end
        else
            switch D
                case 1
                    f_tmp{k}(:, :) = f( ( ( par_breaks(k)-1):(par_breaks(k+1) ) ),...
                                             index{:} );
                case 2
                    f_tmp{k}(:, 2:end-1, :) = f( ( ( par_breaks(k)-1):(par_breaks(k+1) ) ),...
                                             index{:} );
                case 3
                    f_tmp{k}(:, 2:end-1, 2:end-1, :) = f( ( ( par_breaks(k)-1):(par_breaks(k+1) ) ),...
                                             index{:} );
            end
        end
    end
end

%%%% Compute the EC curves
if strcmp( version, "C" )
    % C based implementation
    if Npar < 2
        for n = 1:nEC
            ECn = EulerCharCrit_c( f( index{:}, n ), cc )';
            ECn = ECn( ECn( :, 2 ) ~= 0, : );
            ECn = ECn( ~isnan( ECn( :, 1 ) ), : );
            ECn = ECn( ECn( :, 1 )~=-Inf, : );

            [ ~, I ]   = sort( ECn( :, 1 ), 'ascend' );
            ECn     = ECn( I, : );
            EC{ n } = [ [ -Inf; ECn( :, 1 ); Inf ], [ L0; L0; L0 + ...
                        cumsum( ECn( :, 2 ) ) ] ];
        end
    else
        for n = 1:nEC
            % parallelize EC computation using parfor
            tmp = cell( [ 1 Npar ] );
            parfor k = 1:Npar
                index2  = repmat( {':'}, 1, D );
                tmp{k} = EulerCharCrit_c( f_tmp{k}( index2{:}, n ), cc )';
            end
            clear index2
            
            ECn = [];
            for k = 1:Npar
                ECn = [ECn; tmp{k}];
            end
            
            ECn = ECn( ECn( :, 2 ) ~= 0, : );
            ECn = ECn( ~isnan( ECn( :, 1 ) ), : );
            ECn = ECn( ECn( :, 1 )~=-Inf, : );

            [ ~, I ]   = sort( ECn( :, 1 ), 'ascend' );
            ECn     = ECn( I, : );
            EC{ n } = [ [ -Inf; ECn( :, 1 ); Inf ], [ L0; L0; L0 + ...
            cumsum( ECn( :, 2 ) ) ] ];
        end
    end
else
    % treat cases split by dimension
    switch D
        case 1
            % Compute signs of horizontal and vertical gradient
            df     = sign( f( 2:end, : ) - f( 1:end-1, : ) );
            minMax = df( 1:end-1, : ) - df( 2:end, : );
            f = f( 2:end-1, : );

            % Get locations of local minima and maxima and there numbe rin each
            % function
            mins = minMax == -2;
            maxs = minMax == 2;

            lmins = sum( mins );
            lmaxs = sum( maxs );

            % Compute the EC curves
            for k = 1:nEC
                dEC = [ [ f( mins ), ones( [ lmins( k ), 1 ] ) ];...
                        [ f( maxs ), -ones( [lmaxs( k ), 1 ] ) ] ];

                % sort crits according to its critical height
                [ ~, I ]    = sort( dEC( :, 1 ), 'ascend' );
                dEC         = dEC( I, : );
                dEC( 1, 2 ) = dEC( 1, 2 ) + L0;
                EC{ k }     = [ [ -Inf; dEC( :, 1 ); Inf ],...
                                [ L0; L0; L0 + cumsum( dEC( :, 2 ) ) ] ];
            end

        case 2

           for l = 1 : nEC
                %%%%%%% Find critical points
                CritMat = [ ];
                for i = 2:( sf( 1 ) + 1 )
                    for j = 2:( sf( 2 ) + 1 )
                            dEC = ECchange( f( i-1:i+1, j-1:j+1, l ), cc );
                            if dEC ~= 0
                                CritMat = [ CritMat; [ f(i, j, l) dEC ] ];
                            end
                    end
                end
                [ ~, I ]   = sort( CritMat( :, 1 ), 'ascend' );

                CritMat = CritMat( I, : );
                clear I
                EC{l}   = [ [ -Inf; CritMat( :, 1 ); Inf ], ...
                            [ L0; L0; L0 + cumsum( CritMat( :, 2 ) ) ] ];
           end
        case 3
            for l = 1 : nEC
                %%%%%%% Find critical points
                CritMat = [ ];
                for i = 2:( sf( 1 ) + 1 )
                    for j = 2:( sf( 2 ) + 1 )
                        for k = 2:( sf( 3 ) + 1 )
                            dEC = ECchange( f( i-1:i+1, j-1:j+1, k-1:k+1, l ), cc );
                            if dEC ~= 0
                                CritMat = [ CritMat; [ f(i, j, k, l) dEC ] ];
                            end
                        end
                    end
                end
                [ ~, I ]   = sort( CritMat( :, 1 ), 'ascend' );

                CritMat = CritMat( I, : );
                clear I
                EC{l}   = [ [ -Inf; CritMat( :, 1 ); Inf ], ...
                            [ L0; L0; L0 + cumsum( CritMat( :, 2 ) ) ] ];
            end
    end
end

% close connection to CPUs
if( state_gcp == 42 && Npar > 1 )   
    delete(gcp)
end

return