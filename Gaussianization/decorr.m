function [corrval, covmate] = decorr( lat_data )
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

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
basefield = mean(lat_data);
varfield = var(lat_data);
basefield.field = basefield.field*0;
nsubj = lat_data.fibersize;
D = lat_data.D;
if D == 1
    for I = 1:nsubj
        basefield.field = basefield.field + lat_data.field(:,I).*mean(lat_data.field(:,[1:(I-1),(I+1):nsubj]), 2);
    end
elseif D == 2
    for I = 1:nsubj
        basefield.field = basefield.field + lat_data.field(:,:,I).*mean(lat_data.field(:,:,[1:(I-1),(I+1):nsubj]),3);
    end 
end
basefield = Mask(basefield);
varfield = Mask(varfield);

corrval = sum(basefield.field(:))/sum(basefield.mask(:))/nsubj;
varval = sum(varfield.field(:))/sum(varfield.mask(:));

covmate = corrval*ones(nsubj);
covmate = covmate + (varval-corrval)*eye(nsubj);

end

