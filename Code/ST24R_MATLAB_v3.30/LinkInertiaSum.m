%% "linkInertiaSum" computes the MANIPULATOR INERTIA MATRIX for the sum of
% two open chain links.
% Use in SE(3).
%
% 	im = LinkInertiaSum(im1,im2)
%
% It gives the MANIPULATOR INERTIA MATRIX (MIM) corresponding to the
% Lagangian's equations of the dynamics of the robot
% formed by this two links. im1 and im2 are the LINK INERTIA MATRIX 
% corresponding to each link. Be awere that you need to repeat this funtion
% by as many links as the robot has for getting the total MIM of a robot.
%
%                           |mI  0|
% IM = Sum(i=1,2)[JstT(th)'*|     |*JstT(th)] = im1+im2   
%                           |0  IT|
% Be careful, because the inertia matrix of each link can have a different
% dimension. Therefore, you must equal
% the size of the links im1 and im2 before to add them.
%
% See also: linkinertia, GeoJacobianT.
%
% Copyright (C) 2001-2019, by Dr. Jose M. Pardos-Gotor.
%
% This file is part of The ST24R "Screw Theory Toolbox for Robotics" MATLAB
% 
% ST24R is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% ST24R is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with ST24R.  If not, see <http://www.gnu.org/licenses/>.
%
% http://www.
%
% CHANGES:
% Revision 1.1  2019/02/11 00:00:01
% General cleanup of code: help comments, see also, copyright
% references, clarification of functions.
%
%% im = linkInertiaSum(im1,im2)
%
function im = LinkInertiaSum(im1,im2)
%
    im = zeros(max(size(im1,1),size(im2,1)),max(size(im1,2),size(im2,2)));
    for i = 1:size(im1,1)
        for j = 1:size(im1,2)
        im(i,j) = im(i,j) + im1(i,j);
        end
    end
    for i = 1:size(im2,1)
        for j = 1:size(im2,2)
        im(i,j) = im(i,j) + im2(i,j);
        end
    end
end
%
    
