function modul(iterand, niterand)
% MODUL( iterand, niterand ) allows you to easily check how a for loop is
% progressing by displaying iterand iff it is evenly divided by niterand
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  iterand    - an integer representing the current iteration number
% Optional
%  niterand   - an integer representing the interval at which to display
%              the iterand value (default: 100)
%--------------------------------------------------------------------------
% OUTPUT
% None (displays the value of iterand if it is evenly divided by niterand)
%--------------------------------------------------------------------------
% EXAMPLES
% for i = 1:1000
%     modul(i, 200);
% end
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

