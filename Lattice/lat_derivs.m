function [ derivs, derivs2, fielditself ] = lat_derivs( data, points )
% LAT_DERIVS( lat_data, points ) takes a lattice of data and estimates the
% derivative at the specified points (Only 1D for now!)
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data     an object of class field
%  points   a D by npoints matrix where each column is a point in the
%           field at which to evaluate the derivative
%--------------------------------------------------------------------------
% OUTPUT
%  derivs   a D by npoints matrix
%  derivs2  a D by D by npoints matrix
%--------------------------------------------------------------------------
% EXAMPLES
% %% 1D Examples
% data = 1:10; y = Field(data', 2); y.xvals = {data};
% [derivs, deriv2, fielditself] = lat_derivs( y, 4 )
%
% x = -1:0.01:1; y = Field(x'.^2, 2); y.xvals = {x};
% [a,b,c] = lat_derivs(y, 0)
%
% %% 2D examples
% [x,y] = meshgrid(1:10,1:10); f = Field(x+y, 2);
% [derivs, deriv2] = lat_derivs( f, [4,4]' )
%
% [x,y] = meshgrid( -1:0.01:1, -1:0.01:1); f = Field(x.^2+y.^2, 2);
% f.xvals = {-1:0.01:1, -1:0.01:1}
% [derivs, deriv2] = lat_derivs( f, [0,0]' )
%
% [x,y] = meshgrid( -1:0.01:1, -1:0.01:1); f = Field(x.^2+ 3*x.*y + y.^2, 2);
% f.xvals = {-1:0.01:1, -1:0.01:1}
% [derivs, deriv2] = lat_derivs( f, [0,0]' )
% [derivs, deriv2] = lat_derivs( f, [-1,-1]' )
% [derivs, deriv2] = lat_derivs( f, [-1,1]' )
%
% f = statfield( 4, 1, 5, 1, 0 )
% [ derivs(I), derivs2(I), fielditself(I) ] = lat_derivs( f, 3 )
% latspacing = (f.xvals{1}(2) - f.xvals{1}(1));
% leftderiv = (f.field(13) - f.field(12))/latspacing;
% rightderiv = (f.field(14) - f.field(13))/latspacing;
% derivsnew(I) = (1/2)*(leftderiv + rightderiv);
% derivs2new(I) = (rightderiv - leftderiv)/latspacing;
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
npoints = size(points,2);
D = data.D;

if D ~= size(points, 1)
    error('The dimensions of the points must be the same as the data.')
end

% Initialize vectors to store the index for each point and the spacing in
% the dth dimension.
point_index = zeros(npoints, D);
lat_spacing = zeros(1, D);

for d = 1:D
    
    % Find the indices of the points that you want to find the deriviatives at
    for I = 1:npoints
        findpoint = find(data.xvals{d} == points(d,I));
        if ~isempty(findpoint)
            point_index(I,d) = findpoint;
        else
            error('Some points don''t exist.')
        end
    end
    
    % Determine the spacing in the dth direction
    lat_spacing(d) = data.xvals{d}(2) - data.xvals{d}(1);
end


%%  Main Function Loop
%--------------------------------------------------------------------------
derivs = zeros(D,npoints);
derivs2 = zeros(D,D,npoints);
fielditself = zeros(1, npoints);

for I = 1:npoints
    fielditself(I) = data.field(point_index(I,:));
    %PIE: point index entry, converted to a cell array to allow for
    %indexing
    PIE = point_index(I,:);
    
    %     if PIE(d) == 1 || PIE(d) == length(data.xvals{d})
    %         error('The off diagonal second derivative has not been implemented for edge voxels')
    %     end
    
    % Derivatives and diagonal 2nd derivatives
    for d = 1:D
        PIE2use = PIE;
        % Hacky fix to ignore the edge cases:
        if PIE(d) == 1
            PIE2use(d) = PIE(d) + 1;
        end
        if PIE(d) == length(data.xvals{d})
            PIE2use(d) = PIE(d) - 1;
        end
        
        sbvec = sbasis(d,D)';
        PIm2dex = num2cell(PIE2use-2*sbvec);
        PIm1dex = num2cell(PIE2use-sbvec);
        PIdex = num2cell(PIE2use);
        PI1dex = num2cell(PIE2use+sbvec);
        PI2dex = num2cell(PIE2use+2*sbvec);
        if PIE2use(d) > 1 && PIE2use(d) < length(data.xvals{d})
            derivs(d,I) = (1/2)*(data.field(PI1dex{:}) - data.field(PIm1dex{:}))/lat_spacing(d);
            derivs2(d,d,I) = (data.field(PI1dex{:}) - 2*data.field(PIdex{:}) + data.field(PIm1dex{:}))/lat_spacing(d)^2;
        elseif PIE(d) == 1
            derivs(d,I) = (data.field(PI1dex{:}) - data.field(PIdex{:}))/lat_spacing(d);
            derivs2(d,d,I) = (data.field(PI2dex{:}) - 2*data.field(PI1dex{:}) + data.field(PIdex{:}))/lat_spacing(d)^2;
        elseif PIE(d) == length(data.xvals{d})
            derivs(d,I) = (data.field(PIdex{:}) - data.field(PIm1dex{:}))/lat_spacing(d);
            derivs2(d,d,I) = (data.field(PIdex{:}) - 2*data.field(PIm1dex{:}) + data.field(PIm2dex{:}))/lat_spacing(d)^2;
        end
    end
    
    % Off diagonal hessian elements
    for d1 = 1:D
        PIE2use = PIE;
        % Hacky fix to ignore the edge cases:
        if PIE(d1) == 1
            PIE2use(d1) = PIE(d1) + 1;
        end
        if PIE(d1) == length(data.xvals{d1})
            PIE2use(d1) = PIE(d1) - 1;
        end
        
        sbvec1 = sbasis(d1,D)';
        for d2 = 1:(d1-1)
            PIE2usehere = PIE2use;
            if PIE2use(d2) == 1
                PIE2usehere(d2) = PIE2use(d2) + 1;
            end
            if PIE2use(d2) == length(data.xvals{d2})
                PIE2usehere(d2) = PIE2use(d2) - 1;
            end
            
            sbvec2 = sbasis(d2,D)';
            PIdex11 = num2cell(PIE2usehere+sbvec1+sbvec2);
            PIdex01 = num2cell(PIE2usehere+sbvec2);
            PIdex10 = num2cell(PIE2usehere+sbvec1);
            PIdex00 = num2cell(PIE2usehere);
            
            % Calculate ((f(x+1,y+1)-f(x,y+1))/h_1 - (f(x+1,y)-f(x,y))/h_1)/h_2
            % (or rather the generalization to D dimensions, the above
            % expression is the 2D one.
            derivs2(d1,d2,I) = ((data.field(PIdex11{:}) - data.field(PIdex01{:})) ...
                - (data.field(PIdex10{:}) - data.field(PIdex00{:})))/(lat_spacing(d1)*lat_spacing(d2));
            derivs2(d2,d1,I) = derivs2(d1,d2,I);
        end
    end
    
end

end

