function LKC = LKC_warp_est( field, normalize )
% LKC_warp_EST( Y, mask, Mboot, normalize, version ) implements the warping
% estimator for LKCs proposed in Taylor et al 2007.
% Currently support is only for 1D and 2D rectangular mask, i.e. all true!
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory:
%  field  an object of class field containing observation of a random field.
%
% Optional:
%  normalize logical indicating whether Y needs to be standardized. 
%            Default 1, i.e., mean will be subtracted and data will be 
%            standardized to have empirical variance 1, if N>1 else 0.
%--------------------------------------------------------------------------
% OUTPUT
%  LKC     structure containing fields:
%          - hatL1: D x N or D x Mboot array containing the LKC estimates
%                   for each random field
%          - hatL: D x 1 vector of estimates of LKC for the sample Y. It
%                  is the average of hatL1.
%          - L0: integer containing the Euler characteristic of the mask
%                equivalently the zeroth LKC of the random fields.
%          - hatSIGMA: D x D estimate of the covariance matrix of the
%                      estimate hatL. It is computed using the empirical
%                      covariance of hatL1.
%          - hatSE:    D x 1 vector of standard errors for components of
%                      hatL based on the CLT. 
%          - confInt95: approximate 95% confidence intervals for hatL
%                       based on the standard CLT.
%--------------------------------------------------------------------------
% EXAMPLES
% -------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

% Get constants from the input
D  = field.D;
N  = field.fibersize;
%f  = Mask( field );
f  = field.field;

mask = field.mask;

% Check that method is implemented for dimension D
if D > 3
    error( 'D must be < 4. Higher dimensional domains have not been implemented')
end

%% Add and check optional values
%--------------------------------------------------------------------------

if ~exist( 'normalize', 'var' )
    % Default value of "normalize"
    if N == 1
       normalize = 0;
    else
       normalize = 1;
    end
end

%% Main function
%--------------------------------------------------------------------------
% Initiate output
LKC = zeros( 1, D );

% If necessary center and normalize the input fields
if( normalize )
    f = f - mean( f, D+1 );
    f = f ./ sqrt( sum( f.^2, D+1 ) );
else
    f = f / sqrt(N);
end

% Calculate the LKCs
switch D
    case 1
        % note that this implementation does not take disconnected domains
        % into account yet!
        f = f(mask, :);
        LKC(1) = sum( sqrt( sum( diff( f, 1, 1).^2, 2 ) ) );
    case 2
        % Compute length of horizontal/vertical edges, the
        % triangulation of the mask is given by the following pattern for
        % 3 x 3 square. This is a special choice and could be replaced by
        % any other triangulation
        %
        % o---o---o---o---o
        % | / | / | / | / |
        % o---o---o---o---o
        % | / | / | / | / |
        % o---o---o---o---o
        edges_vert = diff( f, 1, 1 );
        edges_horz = diff( f, 1, 2 );
        l_edges_vert = sqrt(sum( edges_vert.^2, 3 ));
        l_edges_horz = sqrt(sum( edges_horz.^2, 3 ));
        % edges_diag = sqrt(sum( (Z(1:(end-1),1:(end-1),:) - Z(2:end,2:end,:)).^2, 3 )) ;

        % Compute estimator L1 = sum_edges mu_1(edges) - sum_triangles mu_1(triangles)
        % Note that the interior cancels out!
        LKC(1) = 0.5*( sum( l_edges_horz(1,:) ) + sum( l_edges_horz(end,:) ) + ...
                       sum( l_edges_vert(:,1) ) + sum( l_edges_vert(:,end) ) ) ;

        % Compute estimator L2 = sum_triangles mu_2(triangles)
        % Taylor (2007) "Detecting Sparse Signals..." p920 eq. (21)
        LKC(2) = sum(sum(...
                    sqrt(...
                            l_edges_vert(:,1:(end-1)).^2 .* l_edges_horz(1:(end-1),:).^2 ...
                            - ( sum( edges_vert(:,1:(end-1),:) .* edges_horz(1:(end-1),:,:), 3 ) ).^2 ...
                         )...
                 )) / 2 + ...
                 sum(sum(...
                    sqrt(...
                            l_edges_vert(:,2:end).^2 .* l_edges_horz(2:end,:).^2 ...
                            - ( sum( edges_vert(:,2:end,:) .* edges_horz(2:end,:,:), 3 ) ).^2 ...
                         )...
                 )) / 2;
end

if D > 2
    error( "This estimator is not yet implemented for D>2." )
end

return