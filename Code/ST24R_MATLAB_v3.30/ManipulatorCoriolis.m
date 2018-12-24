%% "ManipulatorCoriolis" computes the MANIPULATOR CORIOLIS MATRIX for an open chain manipulator.
% Use in SE(3).
%
% 	MC = ManipulatorCoriolis(x, t, dt, g0, m, it)
%
% Gives the MANIPULATOR CORIOLIS MATRIX (C) corresponding to the Lagangian's equations: 
% M(t)*ddt + C(t,dt)*dt + N(t,dt) = T
% of the dynamics of the robot formed by links on an open chain. 
%
% "x" is the matrix of the twists which affect the links of the manipulator 6xn.
% "t" is the matrix of the magnitudes "theta" of the respective twists 1xn.
% "dt" is the matrix of the velocities "derive theta" of the respective twists 1xn.
% "g0" is the matrix of the homogeneous transformations which represent the links' center of masses' locations,
% at the reference configuration 4x4xn (three dimension matrix).
% "m" is the matrix of the masses of each link of the manipulator 1xn.
% "it" is the matrix of the inertia tensors of the links 3x3xn
% (three dimension matrix).
% with a number of links n.
%
%     |C11...C1n|
% C = |         |, With Cij = 1/2 * Sum(1,n)[( dMij/dtk + dMik/dtj - dMkj/dti ) * dtk]  
%     |Cn1...Cnn|
% where:
% dMij/dtk = Sum(l=max(i,j),n)[[Ak-1i*Ei,Ek]'*Alk'*Ml*Alj*Ej + Ei'*Ali'*Ml*Alk*[Ak-1j*Ej,Ek]]
% and so on and so forth.
% With Ml being the link inertia transformed matrix of the link l 6x6, 
% incapsulated into a three dimensions matrix 6x6xl.
% With Ei being the twist xi 6x1.
% With Aij being an element 6x6 of the Manipulator adjoint transformation
% A 6x6xixj (four dimensions).
%
% See also: manipulatorinertia, manipulatoradjoint, linkinertiatrans,
% linkinertia.
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
%% cm = ManipulatorCoriolis(x, t, dt, g0, m, it)
%
function cm = ManipulatorCoriolis(x, t, dt, g0, m, it)
%
        ma = ManipulatorAdjoint(x,t);
        n = size(x,2);
        mt = zeros(6,6,n);
        mc = zeros(n,n);
        for i=1:n
            mt(:,:,i) = LinkInertiaTrans(g0(:,:,i),m(i),it(:,:,i));    
        end     
        for i=1:n
            for j=1:n
                for k=1:n
                    a = 0;
                    b = 0;
                    c = 0;
                    for p=max([i j]):n
                            a = a + twistbracket(ma(:,:,k,i)*x(:,i),x(:,k))'*ma(:,:,p,k)'*mt(:,:,p)*ma(:,:,p,j)*x(:,j);
                            a = a + x(:,i)'*ma(:,:,p,i)'*mt(:,:,p)*ma(:,:,p,k)*twistbracket(ma(:,:,k,j)*x(:,j),x(:,k));
                    end
                    for p=max([i k]):n
                            b = b + twistbracket(ma(:,:,j,i)*x(:,i),x(:,j))'*ma(:,:,p,j)'*mt(:,:,p)*ma(:,:,p,k)*x(:,k);
                            b = b + x(:,i)'*ma(:,:,p,i)'*mt(:,:,p)*ma(:,:,p,j)*twistbracket(ma(:,:,j,k)*x(:,k),x(:,j));
                    end
                    for p=max([k j]):n
                            c = c + twistbracket(ma(:,:,i,k)*x(:,k),x(:,i))'*ma(:,:,p,i)'*mt(:,:,p)*ma(:,:,p,j)*x(:,j);
                            c = c + x(:,k)'*ma(:,:,p,k)'*mt(:,:,p)*ma(:,:,p,i)*twistbracket(ma(:,:,i,j)*x(:,j),x(:,i));
                    end
                    mc(i,j) = mc(i,j) + (a+b-c)*dt(k);
                end
            end
        end
        cm = mc * 0.5;
end
%