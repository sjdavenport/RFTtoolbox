function modul( iterand, niterand )
% MODUL( iterand, niterand ) allows you to easily check how a for loop is
% progressing by displaying iterand iff it is evenly divided by niterand
%--------------------------------------------------------------------------
% ARGUMENTS
% iterand   
% niterand 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('niterand', 'var')
    % Set the default option for niterand
    niterand = 100;
end

% Test to see whether the iterand is evenly divided by the niterand
if mod(iterand, niterand) == 0
    disp(iterand);
end

end

