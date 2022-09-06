function out = mtimes( obj1, obj2 )
% mtimes( obj1, obj2 ) applies matrix multiplication within the fiber of a
% object of class Field and a constant vector or another object of class 
% Field with compatible dimension.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj1   object of type Field with fiberD < 3 or a vector.
%  obj2   object of type Field with fiberD < 3 or a vector.
%
%--------------------------------------------------------------------------
% OUTPUT
% out  an object of class Field with fiber obtained by the matrix
%      multiplication.
%
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input and get constants
%--------------------------------------------------------------------------

if isa( obj1, 'Field' ) || isa( obj2, 'Field' )
    if isa( obj1, 'Field' )
        sfiberA = obj1.fibersize;
        if length( sfiberA ) == 1
            sfiberA = [ 1 sfiberA ];
        end
        msizeA  = obj1.masksize;
        A       = obj1.field;
        out     = obj1;
        % Get index for the domain
        index  = repmat( {':'}, 1, obj1.D );
        Aisfield = 1;
        
    elseif isnumeric( obj1 )
        msizeA  = NaN;
        sfiberA = size( obj1 );
        A       = obj1;
        Aisfield = 0;
        
    else
        error( "obj1 must be either a numeric vector or matrix or a Field object." )
    end
    clear obj1;
    
    if isa( obj2, 'Field' )
        sfiberB = obj2.fibersize;
        if length( sfiberB ) == 1
            sfiberB = [ sfiberB 1 ];
        end
        msizeB  = obj2.masksize;
        B       = obj2.field;
        out     = obj2;
        % Get index for the domain
        index  = repmat( {':'}, 1, obj2.D );
        Bisfield = 1;
        
    elseif isnumeric( obj2 )
        sfiberB = size( obj2 );
        B       = obj2;
        msizeB  = NaN;
        Bisfield = 0;
    else
        error( "obj1 must be either a numeric vector or matrix or a Field object." )
    end
    clear obj2;
    
else
    error( "One of the inputs must be of class Field." )
end

% Check whether the fibers are valid for matrix multiplications
if length( sfiberA ) ~= 2 || length( sfiberB ) ~= 2
    error( strcat( "The fiber dimension of the input must be",...
                   " below 2 for carrying out matrix multiplications." ) )
end

% Check whether obj and v have fiber dimensions compatible for matrix
% multiplication
if sfiberA(2) ~= sfiberB(1) && ( ~all( sfiberA == 1 ) && ~all( sfiberB == 1 ) )
        error( "The dimensions are not suitable for matrix multiplications." )
end

% Check whether size of domains match
if all( ~isnan( msizeA ) ) && all( ~isnan( msizeB ) )
    if length( msizeA ) ~= length( msizeB )
        error( "The input must have the same dimension of the domain." )
    elseif ~all( msizeA == msizeB )
        error( "The input must have the same size of the domain." )
    end
end

%% Main function
%--------------------------------------------------------------------------

if all( ~isnan( msizeB ) )
    out.field = squeeze( zeros( [ msizeB, sfiberA(1), sfiberB(2) ] ) );
else
    out.field = squeeze( zeros( [ msizeA, sfiberA(1), sfiberB(2) ] ) );
end

% if (all( sfiberA == 1 ) && Aisfield == 0) || (all( sfiberB == 1 ) && Bisfield == 0)
if (all( sfiberA == 1 ) && Aisfield == 0) || (all( sfiberB == 1 ) && Bisfield == 0)
    out.field = A*B;

elseif ~all( isnan( msizeA ) ) && ~all( isnan( msizeB ) )
    if all( sfiberA ~= 1 ) && all( sfiberB ~= 1 ) % Matrix Matrix multiplication
        for d = 1:sfiberA(1)
            for dd = 1:sfiberB(2)
                for k = 1:sfiberB(1)
                    out.field( index{:}, d, dd ) = out.field( index{:}, d, dd ) +...
                                           A( index{:}, d, k ) .* B( index{:}, k, dd );
                end
            end
        end

    elseif all( sfiberA ~= 1 ) && sfiberB(2) == 1 % Matrix vector multiplication
        for d = 1:sfiberA(1)
            for k = 1:sfiberA(2)
                out.field( index{:}, d ) = out.field( index{:}, d ) +...
                                             A( index{:}, d, k ) .* B( index{:}, k );
            end
        end
    elseif sfiberA(1) == 1 && all( sfiberB ~= 1 )  % Vector Matrix multiplication
        for d = 1:sfiberB(2)
            for k = 1:sfiberA(2)
                out.field( index{:}, d ) = out.field( index{:}, d ) +...
                                             A( index{:}, k ) .* B( index{:}, k, d );
            end
        end
    elseif sfiberA(1) == 1 && sfiberB(2) == 1 % vector vector multiplication
        for k = 1:sfiberA(2)
            out.field( index{:} ) = out.field( index{:} ) +...
                                         A( index{:}, k ) .* B( index{:}, k );
        end        
    end

elseif all( isnan( msizeB ) )
    if all( sfiberA ~= 1 ) && all( sfiberB ~= 1 ) % Matrix Matrix multiplication
        for d = 1:sfiberA(1)
            for dd = 1:sfiberB(2)
                for k = 1:sfiberB(1)
                    out.field( index{:}, d, dd ) = out.field( index{:}, d, dd ) +...
                                           A( index{:}, d, k ) .* B( k, dd );
                end
            end
        end

    elseif all( sfiberA ~= 1 ) && sfiberB(2) == 1 % Matrix vector multiplication
        for d = 1:sfiberA(1)
            for k = 1:sfiberA(2)
                out.field( index{:}, d ) = out.field( index{:}, d ) +...
                                             A( index{:}, d, k ) .* B( k );
            end
        end
    elseif sfiberA(1) == 1 && all( sfiberB ~= 1 )  % Vector Matrix multiplication
        for d = 1:sfiberB(2)
            for k = 1:sfiberA(2)
                out.field( index{:}, d ) = out.field( index{:}, d ) +...
                                             A( index{:}, k ) .* B( k, d );
            end
        end
    elseif sfiberA(1) == 1 && sfiberB(2) == 1 % vector vector multiplication
        for k = 1:sfiberA(2)
            out.field( index{:} ) = out.field( index{:} ) +...
                                         A( index{:}, k ) .* B( k );
        end        
    end
    
elseif all( isnan( msizeA ) )
    if all( sfiberA ~= 1 ) && all( sfiberB ~= 1 ) % Matrix Matrix multiplication
        for d = 1:sfiberA(1)
            for dd = 1:sfiberB(2)
                for k = 1:sfiberB(1)
                    out.field( index{:}, d, dd ) = out.field( index{:}, d, dd ) +...
                                           A( d, k ) .* B( index{:}, k, dd );
                end
            end
        end

    elseif all( sfiberA ~= 1 ) && sfiberB(2) == 1 % Matrix vector multiplication
        for d = 1:sfiberA(1)
            for k = 1:sfiberA(2)
                out.field( index{:}, d ) = out.field( index{:}, d ) +...
                                             A( d, k ) .* B( index{:}, k );
            end
        end
    elseif sfiberA(1) == 1 && all( sfiberB ~= 1 )  % Vector Matrix multiplication
        for d = 1:sfiberB(2)
            for k = 1:sfiberA(2)
                out.field( index{:}, d ) = out.field( index{:}, d ) +...
                                             A( k ) .* B( index{:}, k, d );
            end
        end
    elseif sfiberA(1) == 1 && sfiberB(2) == 1 % vector vector multiplication
        for k = 1:sfiberA(2)
            out.field( index{:} ) = out.field( index{:} ) +...
                                         A( k ) .* B( index{:}, k );
        end        
    end
end

return