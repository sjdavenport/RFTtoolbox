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
% -------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow  
%--------------------------------------------------------------------------


%% ------------------------------------------------------------------------
%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
%%%%%% HEADLINE 
%%% SUBHEADLINE
% line or short paragraph comment

%% ------------------------------------------------------------------------
%  add/check optional input
%--------------------------------------------------------------------------

% this kind of code with exists is better than using nargin < xy, since
% parameter can be easily permuted
if ~exist( 'opt1', 'var' )
   % default option of opt1
   opt1 = 0;
end


%% ------------------------------------------------------------------------
%  main function
%--------------------------------------------------------------------------

return