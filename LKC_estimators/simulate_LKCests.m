function results = simulate_LKCests( Msim, nsubj, methods, params,...
                                     data_gen )
% simulate_LKCests( Msim, methods, data_gen ) simulates different LKC
% estimators for random fields.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   Msim     an integer for the number of monte carlo simulations
%   nsubj    an integer for the number of samples in the estimates
%   methods  a structure possibly containing the fields
%            - convE filled with the 'version' variable from LKC_est() 
%            - HPE filled with a logical for normalize from LKC_HP_est()    
%            - bHPE filled with a 1 x 2 cell containing Mboot, normalize 
%              from LKC_HP_est()
%            - warpE filled with the normalize variable
%            - kiebelE filled with the k variable
%   params   an object of type ConvFieldParams
%
% Optional
%   data_gen  a function handle for the lattice data or a logical mask. If
%             it is a logical mask the lattice data is white noise with
%             mean zero and variance one.
%
%--------------------------------------------------------------------------
% OUTPUT
%   results  a structure containing the results of the simulation
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow  
%--------------------------------------------------------------------------


%% Check mandatory input and get important constants
%--------------------------------------------------------------------------
D = length( params.kernel.adjust );

if islogical( data_gen )
    data_gen = @(n) wfield( data_gen, n );
end

if isa( data_gen, 'function_handle' )
    tmp = data_gen(1);
    if tmp.D ~= D
        error( "The smoothing Kernel and the domain need to have the same dimension.")
    end
elseif isa( data_gen, 'Field' )
    index = repmat( {':'}, [ 1 D ] );
    if ~( data_gen.fiberD + data_gen.D ) == D + 1
        error( "data_gen must be an array with D+1" )
    end
else
    error( "data_gen must be either a Field, a function handle or a mask." )
end


%% Main function  
%--------------------------------------------------------------------------

% Preallocate values for the estimates and set default values, if not
% provided
if isfield( methods, "convE" )
    L_conv_ests            = NaN * ones( [ Msim D ] );
    L_conv_ests_nonstatInt = NaN * ones( [ Msim 1 ] );
    version = methods.convE;
end

if isfield( methods, "stationaryE" )
    L_stationary_ests = NaN * ones( [ Msim D ] );
    version2 = methods.stationaryE;
end

if isfield( methods, "HPE" )
    L_HP_ests    = NaN * ones( [ Msim D ] );
    normalizeHPE = methods.HPE(1);
end

if isfield( methods, "bHPE" )
    L_bHP_ests    = NaN * ones( [ Msim D ] );
    normalizebHPE = methods.bHPE(2);
    Mboot         = methods.bHPE(1);
end

if isfield( methods, "warpE" )
    L_warp_ests    = NaN * ones( [ Msim D ] );
    normalizewarp = methods.warpE;
end

if isfield( methods, "kiebelE" )
    L_kiebel_ests    = NaN * ones( [ Msim D ] );
    k = methods.kiebelE;
end

if isfield( methods, "formanE" )
    L_forman_ests  = NaN * ones( [ Msim D ] );
    L_kiebel2_ests = NaN * ones( [ Msim D ] );
    k = methods.formanE;
end

tic
for m = 1:Msim
    if isa( data_gen, 'function_handle' )
        lat_data = data_gen( nsubj );
    elseif isa( data_gen, 'Field' )
        lat_data = Field( data_gen.mask );
        lat_data.field = data_gen.field( index{:},...
                        randsample( data_gen.fibersize, nsubj ) );
    end

    % Generate convolution fields from lattice data
    [ cfield, ss ] = convfield( lat_data, params, 0 );
    cfield = ss .* cfield;

    % Compute derivatives if convE is used in this simulation
    if isfield( methods, "convE" )
        dcfield  = ss .* convfield( lat_data, params, 1 );
        if dcfield.D == 3
            if version(3) == 1
                d2cfield = ss .* convfield( lat_data, params, 2 );
            else
                d2cfield = Field();
            end
        else
            d2cfield = Field();
        end
        [ L_conv, ~, nonstatInt ] = LKC_voxmfd_est( cfield, dcfield,...
                                                     d2cfield,...
                                                     version );
        L_conv_ests( m, : ) = L_conv;
        L_conv_ests_nonstatInt( m ) = nonstatInt;
    end
    
    if isfield( methods, "stationaryE" )
        if ~exist( 'dcfield', 'var' )
            dcfield  = convfield( lat_data, params, 1 );
        end
        tmp = LKC_stationary_est( cfield, dcfield, version2 );
        L_stationary_ests(m,:) = tmp;
    end    

    if isfield( methods, "HPE" )
        tmp = LKC_HP_est( Mask(cfield), 1, normalizeHPE );
        L_HP_ests(m,:)   = tmp.hatL;
    end

    if isfield( methods, "bHPE" )
        tmp = LKC_HP_est( Mask(cfield), Mboot, normalizebHPE );
        L_bHP_ests(m,:)  = tmp.hatL;
    end
    
    if isfield( methods, "warpE" )
        tmp = Mask(cfield);
        s = sqrt( length(tmp.field(tmp.field~=0))/tmp.fibersize);
        cfield = Field( true( [ s s ] ) );
        cfield.field = reshape( tmp.field(tmp.field~=0), [s s tmp.fibersize] );
        
        L_warp_ests(m,:) = LKC_warp_est( Mask(cfield), normalizewarp);
    end
    
    if isfield( methods, "kiebelE" )
        L_kiebel_ests(m,:) = LKC_kiebel_est( cfield, k );
    end
    
    if isfield( methods, "formanE" )
       [L_f, L_k] = LKC_forman_est( cfield, k );
       L_forman_ests(m,:)  = L_f;
       L_kiebel2_ests(m,:) = L_k;
    end
end
simtime = toc;

% Prepare output
results = struct( 'simtime', simtime, 'nsubj', nsubj );

if isfield( methods, "convE" )
    results.convE = L_conv_ests;
    results.convE;
    results.nonstatInt = L_conv_ests_nonstatInt;
end

if isfield( methods, "stationaryE" )
    results.stationaryE = L_stationary_ests;
end

if isfield( methods, "HPE" )
    results.HPE = L_HP_ests;
end

if isfield( methods, "bHPE" )
    results.bHPE = L_bHP_ests;
end

if isfield( methods, "warpE" )
    results.warpE = L_warp_ests;
end

if isfield( methods, "kiebelE" )
    results.kiebelE = L_kiebel_ests;
end

if isfield( methods, "formanE" )
    results.formanE = L_forman_ests;
    results.kiebel2E = L_kiebel2_ests;
end

return