%% "LinkInertia" computes the LINK INERTIA MATRIX for a robot's link.
% Use in SE(3).
%
% 	im = LinkInertia([x1...xn],[t1...tn],gsi0,m,IT)
%
% It gives the LINK INERTIA MATRIX corresponding to the Lagangian's
% equations of the dynamics of an Open chain manipulator.
% Beware that this function does compute ONLY the inertia matrix of a link,
% thus for getting the total manipulator Inertia
% matrix, you need to add all links inertia matrices.
%
% The pairs {xi,ti} define the joint twist "xi" and angle "ti" for each
% joint of the manipulator which affects to the link.
% xi is a vector 6x1, ti is a value, so the [x1...xn] is a matrix 6xn and
% [t1...tn] is a transpose vector 1xn.
% IT is the Inertia Tensor for the link 3x3. m is a the mass of the link.
% the Inertia matrix "im" is 6xn.
% "gsi0" = gsi(0), is th location of the link's center of mass at the 
% reference configuration.
% "gsi" and "gsi0" are homogeneous matrix 4x4.           
% gsi(theta)=exp(E1^theta1)*...*exp(En^thetan)*gsi(0) 
%                |mI  0|
% im = JstT(th)'*|     |*JstT(th)=JstT(th)'*M*JstT(th)= Link Inertia Matrix   
%                |0  IT|
% with M Generalized Inertia Matrix, defined through the diagonal mass and
% the inertia tensor.
% and JstT as the manipulator jacobian on the body coordinate frame.
%
% See also: GeosJacobianS, GeosJacobianT.
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
%% im = LinkInertia(x,t,g0,m,it)
%
function im = LinkInertia(x,t,g0,m,it)
%
    gim = [m*eye(3), zeros(3); zeros(3), it];
    jsib = GeoJacobianT(x,t,g0);
    im = jsib'*gim*jsib;
end
%
    
