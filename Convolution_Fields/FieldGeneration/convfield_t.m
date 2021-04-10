  function [tcfield, cfields] = convfield_t( lat_data, params )
% CONVFIELD_T( lat_data, FWHM, resadd ) computes a convolution t field
% with specified FWHM with resadd additional voxels added between the
% original points of the lattice.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  a Dim by nsubj array of data
%  params    an object of class ConvFieldParams. In D dimensions to specify
%            a given FWHM and resadd set 
%               params = ConvFieldParams(repmat(FWHM, 1, D), resadd)
%            If params is numeric, then FWHM = params and the t-statistic 
%            field, smoothing with this FWHM is obtained as a default, the
%            equivalent to setting 
%            params = ConvFieldParams(repmat(FWHM, 1, D), 0)
%--------------------------------------------------------------------------
% OUTPUT
% tcfield  an object of class Field representing the convolution tfield
% cfields  an object of class Field representing the convolution fields of
%          the lat_data
%--------------------------------------------------------------------------
% EXAMPLES
% %% 1D convolution t field
% nvox = 10; nsubj = 20; resadd = 20; FWHM = 2;
% lat_data = normrnd(0,1,[nvox,nsubj]);
% params = ConvFieldParams(FWHM, resadd);
% tcfield = convfield_t(lat_data, params);
% lattice_tfield = convfield_t(lat_data, FWHM);
% plot(1:nvox, lattice_tfield.field, 'o-')
% hold on
% plot(tcfield.xvals{1}, tcfield.field)
% title('1D convolution t field')
% legend('Convolution field', 'Lattice Evaluation', 'Location', 'Best')
% xlabel('voxels')
% ylabel('t field')
% 
% %% 2D convolution field
% Dim = [10,10]; nsubj = 20; resadd = 5; FWHM = 2;
% params = ConvFieldParams([FWHM, FWHM], resadd);
% lat_data = normrnd(0,1,[Dim,nsubj]);
% tcfield = convfield_t(lat_data, params);
% lattice_tfield = convfield_t(lat_data, FWHM);
% zmax = max(tcfield.field(:)); zmin = min(tcfield.field(:));
% subplot(1,2,1); imagesc(lattice_tfield); zlim([zmin,zmax])
% title('2D lattice t field')
% xlabel('x'); ylabel('y'); zlabel('t field')
% subplot(1,2,2); imagesc(tcfield); zlim([zmin,zmax])
% title('2D convolution t field')
% xlabel('x'); ylabel('y'); zlabel('t field')
% 
% %% 2D convolution field (with mask)
% Dim = [3,3]; nsubj = 20; resadd = 3; FWHM = 8;
% params = ConvFieldParams([FWHM, FWHM], resadd);
% mu = zeros(Dim); mu(1,1) = 10; mu(1,3) = 10;
% mask = true(Dim); mask(2,1) = 0;
% lat_data = wfield(mask, nsubj);
% tcfield = convfield_t(lat_data, params);
% surf(tcfield.field)
% title('2D convolution t field')
% xlabel('x'); ylabel('y'); zlabel('t field')
% 
% %% 2D convolution field (with mask)
% Dim = [10,10]; nsubj = 20; resadd = 3; FWHM = 8;
% params = ConvFieldParams([FWHM, FWHM], resadd);
% mask = true(Dim); mask(4:7,1:5) = 0;
% lat_data = wfield(mask, nsubj);
% tcfield = convfield_t(lat_data, params);
% surf(tcfield.field)
% title('2D convolution t field')
% xlabel('x'); ylabel('y'); zlabel('t field')
% 
% %% 3D convolution t field
% Dim = [11,11,11]; nsubj = 50; resadd = 3; FWHM = 3;
% params = ConvFieldParams([FWHM, FWHM, FWHM], resadd);
% mask = true(Dim);
% lat_data = wfield(mask, nsubj);
% tcf = convfield_t( lat_data, params );
% surf(tcf.field(:,:,30))
% title('3D convolution t field slice')
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Allow for non field input
if ~isa( lat_data, 'Field' ) && isnumeric(lat_data)
    temp_lat_data = lat_data;
    s_lat_data = size(lat_data);
    s_lat_data = s_lat_data(1:end-1);
    if length(s_lat_data) == 1
        s_lat_data = [s_lat_data, 1];
    end
    lat_data = Field(true(s_lat_data));
    lat_data.field = temp_lat_data;
    clear temp_lat_data;
end

% If params is numeric take it to be the FWHM and choose resadd = 0,
% enlarge = 1 and lat_masked equal to 1.
do_smooth = 1;
if isnumeric(params)
    if params == 0
        do_smooth = 0;
    else
        params = ConvFieldParams( repmat(params,1,lat_data.D), 0 );
    end
end

% Get the dimensions of the data
% Dim = lat_data.masksize;
% D = lat_data.D;

%%  Main function
%--------------------------------------------------------------------------
if do_smooth == 1
    % Obtain the convolution fields for each subject
    cfields  = convfield( lat_data, params, 0 );
    
    % Generate an empty Field with the given mask
    tcfield = Field( cfields.mask );
    
%   Calculate the t-statistic
    tcfield.field = mvtstat( cfields.field );
%     if D > 1
%         tcfield.field = mvtstat( cfields.field, spacep( Dim, params.resadd ) + 2* params.enlarge );
%     else
%         tcfield.field = mvtstat( cfields.field );
%     end
    
    tcfield.xvals = cfields.xvals;
else
    % Generate an empty Field with the given mask
    tcfield = Field( lat_data.mask );
    
    % Calculate the t-statistic
    tcfield.field = mvtstat(lat_data.field);
    
    % Set the xvals to be the xvals of the original data
    tcfield.xvals = lat_data.xvals;
end

% Ensure that the resulting field is masked
tcfield = Mask(tcfield);

end

