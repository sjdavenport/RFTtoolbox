function det_array = vectdet( matrix_array )
% vectdet( matrix_array ) computes the voxelwise determinant of an array
%--------------------------------------------------------------------------
% ARGUMENTS
% matrix_array      a D by D by nvals array
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% %% 1D
% array = zeros(1,2);
% array(1,1) = 2;
% array(1,2) = 4;
% 
% vectdet(array)
% 
% %% 2D 
% array = zeros(2,2,2);
% array(:,:,1) = eye(2);
% 
% vectdet(array)
% 
% %% 2D
% array = zeros(2,2,2);
% array(:,:,1) = [1,2;3,4];
% array(:,:,2) = mvnrnd([0,0], eye(2), 2);
% det(array(:,:,2))
% vectdet(array)
% 
% %% 2D
% array = zeros(2,2,2);
% array(:,:,1) = eye(2);
% array(:,:,2) = mvnrnd([0,0], eye(2), 2);
% det(array(:,:,2))
% vectdet(array)
% 
% %% 3D
% array = zeros(3,3,2);
% array(:,:,1) = eye(3);
% array(:,:,2) = mvnrnd([0,0,0], eye(3), 3);
% det(array(:,:,2))
% vectdet(array)
% 
% %%
% Dim = [3,3];
% lat_data = noisegen([91,109,91], 9, 2);
% reshaped_data = reshape(lat_data, [3,3,91*109*91]);
% reshaped_data_indiv = reshape(lat_data, [3,3,91,109,91]);
% 
% tic
% det_reshaped_data = vectdet(reshaped_data);
% toc
% det_reshaped_data(1)
% det(reshaped_data(:,:,1))
% 
% tic
% det_notparallel = zeros(1,91*109*91);
% for I = 1:91
%     for J = 1:109
%         for K = 1:91
%             det_notparallel = det(reshaped_data_indiv(:,:,I,J,K));        
%         end
%     end
% end
% toc
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
s_array = size(matrix_array);
D = s_array(1);
det_array = zeros(1, s_array(end));

if D == 1
    det_array = matrix_array;
else
    for d = 1:D
        allbutd = setdiff(1:D, d);
        ndbit = squeeze(vectdet(matrix_array(2:D, allbutd, :)));
        if size(ndbit,2) == 1
            ndbit = ndbit';
        end
        det_array = ((-1)^(d-1))*squeeze(matrix_array(1,d,:))'.*ndbit + det_array;
    end
end

end

