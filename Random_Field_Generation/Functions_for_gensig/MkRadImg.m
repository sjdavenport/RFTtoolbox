function Rmap = MkRadImg(Dim,c)
% MKRADIMG calculates the ellipse scaled distance of points in the image
% given by the array: (1:Dim[1], 1:Dim[2], 1:Dim[3]) to the point c.
%--------------------------------------------------------------------------
% USAGE
% Rmap = MKRADIMG(Dim,c)
%--------------------------------------------------------------------------
% ARGUMENTS
% Dim   A 1 by 2 or 1 by 3 vector of the image dimensions, 
%       Dim = [256,256,256] corresponds to a 256*256*256 image.
% c     A 1 by 2 or 1 by 3 vector which specifies the centre that you'd 
%       like to build your ellipse around.
%--------------------------------------------------------------------------
% OUTPUT
% Rmap  A 3d array with dimensions given by Dim such that the value at the 
%       (i,j,k)th coordinate is (i - c(1))^2 + (j - c(2))^2 + (k - c(3))^2 
%       for i in 1:Dim[1], j = 1:Dim[2] and k = 1:Dim[3].
%--------------------------------------------------------------------------
% EXAMPLES
% MkRadImg([2,2,2], [0,0,0])
% %compare to sqrt(4+4+225*4)
%
% Img = double(MkRadImg([256, 256],[128.5,128.5]) <= 20);
% surf(Img);
%
% Img = double(MkRadImg([256, 256, 256],[128.5,128.5, 128.5]) <= 20);
% surf(Img(:,:,128)); %And compare this to:
% surf(Img(:,:,109)); %The radius is smaller because we're intersecting a
%                     %higher part of the sphere.
%--------------------------------------------------------------------------

%Note that if [x,y] = ndgrid(a,b) for vectors a, b of lengths na, nb
% then x is an na by nb matrix such that each every value in the ith row of
% x is a[i]. Ie x is the vector a as a column repeated nb times.
% And y is also an na by nb matrix, but each row of y is equal to b. Ie the
% columns are constant.

%The point is the x doesn't vary in the y coordinate, only in the x
%coordinate.

%Ie if [x,y,z] = ndgrid(a,b,c) 
%then x(i,j,k) = a(i) for all j,k
%and y(i,j,k) = b(j) for all i and k
%and z(i,j,k) = c(k) for all i and j.

%Each of x,y and z is a na by nb by nc array.
%In paricular for vectors starting from 1 up to an integer,
% x(i,j,k) = i, y(i,j,k) = j and z(i,j,k) = k. So x,y, and z are just
% indexing the coordinates. Doing this instead of for loops probably means 
%that things are done much faster.
 if length(Dim) == 2
    [x, y] = ndgrid(1:Dim(1),1:Dim(2));
    Rmap = sqrt((x-c(1)).^2 + (y-c(2)).^2);
 elseif length(Dim) == 3
    [x, y, z] = ndgrid(1:Dim(1),1:Dim(2),1:Dim(3));
    Rmap = sqrt((x-c(1)).^2 + (y-c(2)).^2 + (z-c(3)).^2);
 else
     error('Dim has the wrong dimensions');
 end
 
 %Rmap the values of sqrt(x^2 + y^2 + 225z^2) at each of the coordinates,
 %This is kind of the distance to c, except for the 225 bit! But is
 %effectively some sort of radius of a 3d ellipse.
 

return
