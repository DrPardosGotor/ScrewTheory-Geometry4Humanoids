%% "linkInertiaTrans" the LINK TRANSFORMED INERTIA MATRIX for robot's link.
% Use in SE(3).
%
% 	im = linkInertiaTrans(gsi0,m,IT)
%
% It gives the LINK TRANSFORMED INERTIA MATRIX corresponding to the inertia
% of the link into the base frame of the manipulator.
%
% IT is the Inertia Tensor for the link 3x3. m is a the mass of the link.
% the Inertia matrix "im" is 6x6.
%
% "gsi0" = gsi(0), is th location of the link's center of mass at the
% reference configuration.
% "gsi" and "gsi0" are homogeneous matrix 4x4.           
% Adgsli0 = Ad(gsli0^-1)= Adjoint transformation of the inverse of the 
% center of mass ref config.
%                |mI  0|
% im = Adgsli0' *|     |* Adgsli0 = Adgsli0' * M * Adgsli0
%                |0  IT|
% im = Link Transformed Inertia Matrix.
% with M Generalized Inertia Matrix, defined through the diagonal mass and
% the inertia tensor.
%
% See also: LinkInertia, tform2adjoint.
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
%% LinkInertiaTrans(g0,m,it)
%
function im = LinkInertiaTrans(g0,m,it)
%
    gim = [m*eye(3), zeros(3); zeros(3), it];
    adi0 = tform2adjoint(inv(g0));
    im = adi0'*gim*adi0;
end
%   
