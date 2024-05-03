function [ out1, out2 ] = newfun( var1, var2, opt1 )
% NEWFUN( var1, var2, opt1 ) is described here [...]
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   var1   this is a mandatory variable.
%   var2   this is a mandatory variable.
%
% Optional
%   opt1   this is an optional parameter. Default 0.
%
%--------------------------------------------------------------------------
% OUTPUT
%   out1  this is the first output
%   out2  this is the second output
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow  
%--------------------------------------------------------------------------


%% Check mandatory input and get important constants
%--------------------------------------------------------------------------
%%%%%% HEADLINE 
%%% SUBHEADLINE
% Line or short paragraph comment

%% Add/check optional values
%--------------------------------------------------------------------------

% This kind of code with exists is better than using nargin < xy, since
% Parameter can be easily permuted
if ~exist( 'opt1', 'var' )
   % Default option of opt1
   opt1 = 0;
end


%% Main function  
%--------------------------------------------------------------------------

return