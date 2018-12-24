%% "ManipulatorAdjoint" Computes the ADJOINT TRANSFORMATION for a list of 
% twists-magnitudes.
% Use in SE(3).
%
% 	A = ManipulatorAdjoint(x,t)
%
% ADJOINT TRANSFORMATION: This is a special notation which gives us a most 
% form of the Adjoint of an open chain manipulator
% We use this notation for an easy calculation of the Manipulator Inertia 
% Matrix and the Manipulator Coriolis Matrix. 
% x is the list of twists 6xn.
% t is the list of magnitudes of the twists Theta 1xn.
%      I                                    if i=j 
% Aij= Ad^-1[(exp(Ej+1,Tj+1)...(exp(Ei,Ti)] if i>j
%      0                                    if i<j
%
% With Aij being an element 6x6 of the Manipulator adjoint transformation
% A 6x6xixj (four dimensions).
% With Ad the adjoint transformation "rigidadjoint". Maps twist vectors 
% to twist vectors. Compute Ad in R^6 from homogeneous 4x4.
%
% See also: manipulatorinertia, manipulatorcoriolis, rigidadjoint.
%
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
%% a = ManipulatorAdjoint(x,t)
%
function a = ManipulatorAdjoint(x,t)
%
    n = size(x,2);
    ax = zeros(4,4,n,n);
    ay = eye(6,6);
    a = zeros(6,6,n,n);
    if n > 1
        for i=2:n
            ax(:,:,i,i-1)=expScrew(x(:,i),t(i));
        end
    end
    if n > 2
        for i=n:-1:3
            for j=i-2:-1:1
            ax(:,:,i,j)=ax(:,:,j+1,j)*ax(:,:,i,j+1);
            end
        end
    end
    for i=1:n
        a(:,:,i,i)=ay;
    end
    if n > 1
        for i=2:n
            for j=1:i-1
            a(:,:,i,j)=inv(tform2adjoint(ax(:,:,i,j)));
            end
        end
    end
end
%
    