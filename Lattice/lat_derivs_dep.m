function [ derivs, derivs2 ] = lat_derivs( data, points )
% LAT_DERIVS( lat_data, points ) takes a lattice of data and estimates the
% derivative at the specified points (Only 1D for now!)
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data     an object of class field
%  points       a cell array where each entry is a point of the data
%--------------------------------------------------------------------------
% OUTPUT
%  derivs
%--------------------------------------------------------------------------
% EXAMPLES
% %% 1D examples
% data = 1:10; y = Field(data', 2); y.xvals = {data};
% [derivs, deriv2] = lat_derivs( y, 4 ), derivs = lat_derivs( y, 4 )
% 
% x = -1:0.01:1; y = Field(x'.^2, 2); y.xvals = {x};
% [a,b] = lat_derivs(y, 0)
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
npoints = size(points,2);
D = data.D;

if D ~= length(points) && D ~= 1
    error('The dimensions of the points must be the same as the data.')
end

% Find the indicies of the points that you want to find the deriviatives at
point_index = zeros(1, npoints);
for I = 1:npoints
    findpoint = find(data.xvals{1} == points(I));
    if ~isempty(findpoint)
        point_index(I) = findpoint;
    else
        error('Some points don''t exist.')
    end
end

lat_spacing = data.xvals{1}(2) - data.xvals{1}(1);

%%  Main Function Loop
%--------------------------------------------------------------------------
derivs = zeros(1,npoints);
derivs2 = zeros(1,npoints);

if D == 1
    for I = 1:npoints
        PIE = point_index(I); %PIE: point index entry
        if PIE > 1 && PIE < length(data.xvals{1})
            derivs(I) = (1/2)*(data.field(PIE+1) - data.field(PIE-1))/lat_spacing;
            derivs2(I) = (data.field(PIE+1) - 2*data.field(PIE) + data.field(PIE-1))/lat_spacing^2;
        elseif PIE == 1
            derivs(I) = (data.field(PIE+1) - data.field(PIE))/lat_spacing;
            derivs2(I) = (data.field(PIE+2) - 2*data.field(PIE+1) + data.field(PIE))/lat_spacing^2;
        elseif PIE == length(data.xvals{1})
            derivs(I) = (data.field(PIE) - data.field(PIE-1))/lat_spacing;
            derivs2(I) = (data.field(PIE) - 2*data.field(PIE-1) + data.field(PIE-2))/lat_spacing^2;
        end
    end
else
    error('Need to code D > 1');
end

% if D == 1
%     weights = ones(1,3);
% elseif D > 3
%     error('Weights for D > 3 have not been coded');
% else
%     weights = ones(3*ones(1,D));
%     if D == 2
%        for I = [1,3]
%            for J = [1,3]
%                weights(I,J) = sqrt(2);
%            end
%        end
%        weights(2,2) = 1;
%     elseif D == 3
%         weights = weights*sqrt(2);
%         for I = [1,3]
%             for J = [1,3]
%                 for K = [1,3]
%                     weights(I,J,K) = sqrt(3);
%                 end
%             end
%         end
%         weights(1,2,2) = 1; weights(2,1,2) = 1; weights(3,2,2) = 1;
%         weights(2,3,2) = 1; weights(2,2,1) = 1; weights(2,2,3) = 1;
%         weights(2,2,2) = 1;
%     end
% end
%    
% for I = 1:npoints
%     point = points{I};
%     box = cell(1,D);
%     for d = 1:D
%         box{d} = (point(I) - 1):(point(I) + 1);
%     end
%     subset = data(box{:});
%     cell_point = num2cell(point);
%     differences = subset - data(cell_point{:});
%     derivs{I} = sum(differences./weights)/(3^D - 1)/lat_spacing; %minus one as don't count the middle point
% end

end

