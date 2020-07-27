function results = simulate_LKCests( Msim, nsubj, methods, cfield_props,...
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
%   cfield_props  a structure containing input to cfield_Field()
%                 - lat_masked
%                 - kernel a SepKernel object
%                 - resadd
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

resadd     = cfield_props.resadd;
kernel     = cfield_props.kernel;
if isnumeric( kernel )
    kernel = SepKernel( cfield_props.D, kernel );
end
D          = length( kernel.adjust );
lat_masked = cfield_props.lat_masked;
enlarge    = cfield_props.enlarge;

%% Add/check optional values
%--------------------------------------------------------------------------

% This kind of code with exists is better than using nargin < xy, since
% Parameter can be easily permuted
if ~exist( 'data_gen', 'var' )
   error( "Please, enter a function handle or a mask for data_gen." )
end

if islogical( data_gen )
    data_gen = @(n) wnfield( data_gen, n );
end

tmp = data_gen(1);
if tmp.D ~= D
    error( "The smoothing Kernel and the domain need to have the same dimension.")
end


%% Main function  
%--------------------------------------------------------------------------

% Preallocate values for the estimates and set default values, if not
% provided
if isfield(methods, "convE")
    L_conv_ests            = NaN * ones( [ Msim D ] );
    L_conv_ests_nonstatInt = NaN * ones( [ Msim 1 ] );
    version = methods.convE;
end

if isfield(methods, "HPE")
    L_HP_ests    = NaN * ones( [ Msim D ] );
    normalizeHPE = methods.HPE(1);
end

if isfield(methods, "bHPE")
    L_bHP_ests    = NaN * ones( [ Msim D ] );
    normalizebHPE = methods.bHPE(2);
    Mboot         = methods.bHPE(1);
end

tic
for m = 1:Msim
    lat_data = data_gen( nsubj );

    % Generate convolution fields from lattice data
    cfield   = convfield_Field( lat_data, kernel, 0, resadd, lat_masked,...
                                enlarge );

    % Compute derivatives if convE is used in this simulation
    if isfield( methods, "convE" )
        dcfield  = convfield_Field( lat_data, kernel, 1, resadd, lat_masked,...
                                    enlarge );
        if dcfield.D == 3
            if version(3) == 1
                d2cfield = convfield_Field( lat_data, kernel, 2, resadd,...
                                            lat_masked, enlarge );
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

    if isfield( methods, "HPE" )
        tmp = LKC_HP_est( Mask(cfield), 1, normalizeHPE );
        L_HP_ests(m,:)   = tmp.hatL;
    end

    if isfield( methods, "bHPE" )
        tmp = LKC_HP_est( Mask(cfield), Mboot, normalizebHPE );
        L_bHP_ests(m,:)  = tmp.hatL;
    end
end
simtime = toc;

% Prepare output
results = struct( 'simtime', simtime, 'nsubj', nsubj );

if isfield( methods, "convE" )
    results.convE = L_conv_ests;
    results.convE;
    nonstatInt = L_conv_ests_nonstatInt;
end

if isfield( methods, "HPE" )
    results.HPE = L_HP_ests;
end

if isfield( methods, "bHPE" )
    results.bHPE = L_bHP_ests;
end

return