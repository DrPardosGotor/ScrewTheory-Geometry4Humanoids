%% "ManipulatorInertia" computes the MANIPULATOR INERTIA MATRIX
% for an open chain manipulator.
% Use in SE(3).
%
% 	MI = ManipulatorInertia(x, t, g0, m, it)
%
% Gives the MANIPULATOR INERTIA MATRIX (M) corresponding to the Lagangian's
% equations: 
% M(t)*ddt + C(t,dt)*dt + N(t,dt) = T
% of the dynamics of the robot formed by links on an open chain. 
%
% "x" matrix of the twists which affect the links of the manipulator 6xn.
% "t" matrix of the magnitudes "theta" of the respective twists 1xn.
% "g0" matrix of the homogeneous transformations which represent the links'
% center of masses' locations, at the reference configuration 4x4xn
% (three dimension matrix).
% "m" is the matrix of the masses of each link of the manipulator 1xn.
% "it" is the matrix of the inertia tensors of the links 3x3xn
% (three dimension matrix).
% with a number of links n.
%      |M11...M1n|
% MI = |         |, With Mij = Sum(l=max(i,j),n)[Ei'*Ali'*Ml*Alj*Ej  
%      |Mn1...Mnn|
% 
% With Ml being the link inertia transformed matrix of the link l 6x6,
% incapsulated into a three dimensions matrix 6x6xl.
% With Ei being the twist xi 6x1.
% With Aij being an element 6x6 of the Manipulator adjoint transformation
% A 6x6xixj (four dimensions).
%
% See also: manipulatoradjoint, linkinertiatrans, linkinertia, GeoJacobianT
% linkinertiasum.
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
%% mi = ManipulatorInertia(x, t, g0, m, it)
%
function mi = ManipulatorInertia(x, t, g0, m, it)
%
        a = ManipulatorAdjoint(x,t);
        n = size(x,2);
        mt = zeros(6,6,n);
        mi = zeros(n,n);
        for i=1:n
            mt(:,:,i) = LinkInertiaTrans(g0(:,:,i),m(i),it(:,:,i));    
        end
        for i=1:n
            for j=1:n
                for k=max(i,j):n
                    mi(i,j) = mi(i,j)+x(:,i)'*a(:,:,k,i)'*mt(:,:,k)*a(:,:,k,j)*x(:,j);    
                end
            end
        end
end
%
        
        