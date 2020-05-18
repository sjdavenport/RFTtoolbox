function dEC = ECchange( Z, cc )
% Computes the change in EC by crossing a voxel value from below by
% computing the change in EC from the added cells of a complex.
%  Z:       array of dim 3 x ... x 3 (D-times)
%  cc:      connectivity
%            D=2: 4 or 8
%            D=3: 6 or 26
% Output:
%  dEC change in Euler characteristic, if crossing center voxel from below.
%
%__________________________________________________________________________
% References:
%__________________________________________________________________________
% Author: Fabian Telschow (ftelschow@ucsd.edu)
% Last changes: 2/05/2019
%__________________________________________________________________________

D = length(size(Z));
dEC = 1;

switch D
    case 2
        % use the property that in 2D 4cc and 8cc are related by
        % multiplying the field with a minus sign
        if cc == 8
            Z = -Z;
        end
        
        % subtract number of edges
        dEC = dEC - sum(Z([2 4 6 8]) > Z(5));
        % add number of faces
        if all( Z([1 2 4]) > Z(5) )
            dEC = dEC + 1;
        end
        if all( Z([2 3 6]) > Z(5) )
            dEC = dEC + 1;
        end
        if all( Z([4 7 8]) > Z(5) )
            dEC = dEC + 1;
        end
        if all( Z([6 8 9]) > Z(5) )
            dEC = dEC + 1;
        end
        
    case 3        
        if cc == 6
            % subtract number of edges
            dEC = dEC - sum(Z([5 11 13 15 17 23]) > Z(14));
            
            % add number of faces
            if all(Z([10 11 13]) > Z(14))
                dEC = dEC + 1;
            end
            if all(Z([11 12 15]) > Z(14))
                dEC = dEC + 1;
            end
            if all(Z([15 17 18]) > Z(14))
                dEC = dEC + 1;
            end
            if all(Z([13 16 17]) > Z(14))
                dEC = dEC + 1;
            end
            
            if all(Z([4 5 13]) > Z(14))
                dEC = dEC + 1;
            end
            if all(Z([5 6 15]) > Z(14))
                dEC = dEC + 1;
            end
            if all(Z([13 22 23]) > Z(14))
                dEC = dEC + 1;
            end
            if all(Z([15 23 24]) > Z(14))
                dEC = dEC + 1;
            end
            
            if all(Z([2 5 11]) > Z(14))
                dEC = dEC + 1;
            end
            if all(Z([5 8 17]) > Z(14))
                dEC = dEC + 1;
            end
            if all(Z([11 23 20]) > Z(14))
                dEC = dEC + 1;
            end
            if all(Z([17 23 26]) > Z(14))
                dEC = dEC + 1;
            end
            
            % subtract number of 3D faces
            for shift_up = [0 9]
                if all(Z([1 2 4 5 10 11 13 14]+shift_up) >= Z(14))
                    dEC = dEC - 1;
                end
                if all(Z([2 3 5 6 11 12 14 15]+shift_up) >= Z(14))
                    dEC = dEC - 1;
                end
                if all(Z([5 6 8 9 14 15 17 18]+shift_up) >= Z(14))
                    dEC = dEC - 1;
                end
                if all(Z([4 5 7 8 13 14 16 17]+shift_up) >= Z(14))
                    dEC = dEC - 1;
                end
            end
        elseif cc ==26
            % subtract number of edges
            dEC = dEC - sum(Z(:) > Z(14));
            

            
            % add number of volumes
        else
            error("Please, specify cc to be equal to 6 or 26.")
        end
                
end
dEC = -dEC;
end