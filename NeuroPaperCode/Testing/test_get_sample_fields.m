%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the get_sample_fields function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
tot_nsubj = 30; nvox = 100; data = normrnd(0,1,nvox, tot_nsubj); 
mask = true(1,nvox)';
spfn = get_sample_fields( data, mask );

%% %% 2D Examples
%% Simple 2D example
tot_nsubj = 50; Dim = [20,30]; data = normrnd(0,1,[Dim,tot_nsubj]); 
mask = true(Dim);
spfn = get_sample_fields( data, mask );
spfn(20).lat_data

%% %% 3D Examples
%% Simple 3D example

%% Resting State Data example
mask = imgload('MNImask');
spfn = get_sample_fields( 'RS_2Block', mask );
a = spfn(20)
