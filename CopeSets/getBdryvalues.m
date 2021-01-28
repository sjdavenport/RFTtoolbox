function [bdry_values] = getBdryvalues(field, bdry_params)
% Linear interpolation of a field to a boundary.
% Input:
%  field:       random field over a domain in R^D, it is an (D+1)-dimensional array,
%               where the last dimension enumerates the samples
%  bdry_params: structure generated from getBdryparams

% Output:
%  bdry_values are the linear interpolated values of field along the
%  boundary.
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Alex Bowring
% Last changes: 10/25/2018
%__________________________________________________________________________

%%%%% Compute parameters of the field input
dim = size(field);
D   = length(dim);

%%%%% Treat cases of 1D and 2D domain of the field differently
switch D
    case 2
        %%%%% Linear Interpolation onto the boundary
        lshift_bdry_values = bdry_params.lshift.w1.*field(bdry_params.lshift.edges) + bdry_params.lshift.w2.*field(bdry_params.lshift.edges(:,[dim(2) 1:dim(2)-1]));
        rshift_bdry_values = bdry_params.rshift.w1.*field(bdry_params.rshift.edges) + bdry_params.rshift.w2.*field(bdry_params.rshift.edges(:,[2:dim(2) 1]));
        ushift_bdry_values = bdry_params.ushift.w1.*field(bdry_params.ushift.edges) + bdry_params.ushift.w2.*field(bdry_params.ushift.edges([dim(1) 1:dim(1)-1],:));
        dshift_bdry_values = bdry_params.dshift.w1.*field(bdry_params.dshift.edges) + bdry_params.dshift.w2.*field(bdry_params.dshift.edges([2:dim(1) 1],:));

        bdry_values = [lshift_bdry_values; rshift_bdry_values; ushift_bdry_values; dshift_bdry_values];

    case 3
        %%%%% Linear Interpolation onto the boundary
        lshift_bdry_values = bdry_params.lshift.w1.*field(bdry_params.lshift.edges) + bdry_params.lshift.w2.*field(bdry_params.lshift.edges(:,[dim(2) 1:dim(2)-1],:));
        rshift_bdry_values = bdry_params.rshift.w1.*field(bdry_params.rshift.edges) + bdry_params.rshift.w2.*field(bdry_params.rshift.edges(:,[2:dim(2) 1],:));
        ushift_bdry_values = bdry_params.ushift.w1.*field(bdry_params.ushift.edges) + bdry_params.ushift.w2.*field(bdry_params.ushift.edges([dim(1) 1:dim(1)-1],:,:));
        dshift_bdry_values = bdry_params.dshift.w1.*field(bdry_params.dshift.edges) + bdry_params.dshift.w2.*field(bdry_params.dshift.edges([2:dim(1) 1],:,:));
        bshift_bdry_values = bdry_params.bshift.w1.*field(bdry_params.bshift.edges) + bdry_params.bshift.w2.*field(bdry_params.bshift.edges(:,:,[dim(3) 1:dim(3)-1]));
        fshift_bdry_values = bdry_params.fshift.w1.*field(bdry_params.fshift.edges) + bdry_params.fshift.w2.*field(bdry_params.fshift.edges(:,:,[2:dim(3) 1]));
        
        bdry_values = [lshift_bdry_values; rshift_bdry_values; ushift_bdry_values; dshift_bdry_values, bshift_bdry_values, fshift_bdry_values];
end
end 
   
     