function lat_data = Gmult( lat_data )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

nsubj = lat_data.fibersize;

Gvector = normrnd(0,1,1,nsubj);

indexD = repmat(':',1,lat_data.D);
for I = 1:nsubj
    index = {indexD, I};
    lat_data.field( index{:} ) = Gvector(I)*lat_data.field( index{:} );
end

end

