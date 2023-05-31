function loader( I, totalI, message )
% loader( I, totalI ) is called with the current iteration and the total 
% number of iterations in the loop, resulting in a progress bar being 
% printed to the command window.
%--------------------------------------------------------------------------
% ARGUMENTS
% I: mandatory input, representing the current iteration of the loop.
% totalI: mandatory input, representing the total number of iterations in the loop.
%--------------------------------------------------------------------------
% EXAMPLES
% for I = 1:10^5
%   loader(I,10^5, 'Total progress:')
% end
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
percentDone = 100 * I / totalI;
msg = sprintf('%3.1f\n', percentDone); %Don't forget this semicolon
if I == 1
    fprintf('-------------------------------------------------------\n')
    if nargin < 3
        message = 'Percent done:   ';
    else
        message = [message,'   '];
    end
    reverseStr = '';
    fprintf(message)
else
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf([reverseStr, msg]);

end

