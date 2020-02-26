function EC = EulerChar(f, u, D, dir)
% EulerChar(f, D, mask, version)
% Computes the Euler characteristic of excursion sets from the vector u.
% Note that this implementation is extremely slow compared to
% EulerCharCrit.m so please use the latter, if there is no other reason. 
%--------------------------------------------------------------------------
% ARGUMENTS
%   f    an array K_1, x ... x K_D x N xnsim of random fields for which the
%        Euler characteristic curves should be computed
%   u    is a vector of thresholds
%   D    an integer containing the dimension of the domain of the random
%        fields
%   dir  1 or -1. If 1 EC of excursion sets above u are computed, if -1
%        EC of excursion sets below u are computed
%   version a string. If "C" the fast C implementation is used (default)
%           if "matlab" a slow matlab only implementation is used.
%--------------------------------------------------------------------------
% OUTPUT
%   EC of excursion sets given by u are computed for all available fields
%   over D dimensional domain. The function uses 4-connectivity for 2D
%   fields and 6-connectivity for 3D fields
%--------------------------------------------------------------------------
% EXAMPLES  
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow, Armin Schwartzman
%--------------------------------------------------------------------------

% rename input, since old code is used
N = D;
% Check input
if N >= 4
    error('N >= 4 not implemented')
end
sz = size(f);
if length( sz ) < N
    error('N is larger than array size f')
end

% Size of EC array
if length( sz ) == N
    sz = [ length(u) 1 ];
else
    sz = [ length(u) sz(N+1:end) ];
end

if ~exist('dir','var')
    dir = 1;
end

if dir ~= 1
    u = -u;
end

EC = zeros(length(u)+1, prod(sz(2:end)));   % Extra element deleted later
A1 = false(size(f));
for j = length(u):-1:1
    EC(j,:) = EC(j+1,:);
    if dir == 1
        A = (f > u(j));
    else
         A = (f < u(j));
    end
    D = A & ~A1;    % set added
    switch N
        case 1
            I = find(any(D, 1));
            vertices = sum(A(:,I), 1);
            edges = sum(A(1:end-1,I) & A(2:end,I), 1);
            EC(j,I) = vertices - edges;
        case 2
            I = find(any(any(D, 1), 2));
            % 4-connectivity
            vertices = sum(sum(A(:,:,I), 1), 2);
            edges    = sum(sum(A(1:end-1,:,I) & A(2:end,:,I), 1), 2) + ...
                       sum(sum(A(:,1:end-1,I) & A(:,2:end,I), 1), 2);
            faces    = sum(sum(A(1:end-1,1:end-1,I) & A(2:end,1:end-1,I) & ...
                            A(1:end-1, 2:end, I) & A(2:end, 2:end, I), 1), 2);
            EC(j,I) = vertices - edges + faces;
        case 3
            I = find(any(any(any(D, 1), 2), 3));
            % 6-connectivity
            vertices = sum(sum(sum(A(:,:,:,I), 1), 2), 3);
            edges = sum(sum(sum(A(1:end-1,:,:,I) & A(2:end,:,:,I), 1), 2), 3) + ...
                    sum(sum(sum(A(:,1:end-1,:,I) & A(:,2:end,:,I), 1), 2), 3) + ...
                    sum(sum(sum(A(:,:,1:end-1,I) & A(:,:,2:end,I), 1), 2), 3);
            faces = sum(sum(sum(A(1:end-1,1:end-1,:,I) & A(2:end,1:end-1,:,I) & ...
                                A(1:end-1, 2:end, :,I) & A(2:end, 2:end, :,I), 1), 2), 3) + ...
                    sum(sum(sum(A(1:end-1,:,1:end-1,I) & A(2:end,:,1:end-1,I) & ...
                                A(1:end-1,:, 2:end, I) & A(2:end,:, 2:end ,I), 1), 2), 3) + ...
                    sum(sum(sum(A(:,1:end-1,1:end-1,I) & A(:,2:end,1:end-1,I) & ...
                                A(:,1:end-1, 2:end, I) & A(:,2:end, 2:end ,I), 1), 2), 3);
            cells = sum(sum(sum(A(1:end-1,1:end-1,1:end-1,I) & A(2:end,1:end-1,1:end-1,I) & ...
                                A(1:end-1, 2:end, 1:end-1,I) & A(2:end, 2:end, 1:end-1,I) & ...
                                A(1:end-1,1:end-1, 2:end, I) & A(2:end,1:end-1, 2:end, I) & ...
                                A(1:end-1, 2:end,  2:end, I) & A(2:end, 2:end,  2:end, I), 1), 2), 3);
            EC(j,I) = vertices - edges + faces - cells;
    end
    A1 = A; % new set becomes old
end
EC = reshape(EC(1:end-1,:), sz);
