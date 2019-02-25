%% "CoriolisAij" CORIOLIS matrix C(t,dt) for an open chain manipulator.
% computation based on the Adjoint Transformation Aij.
% Use in SE(3).
%
% 	Ctdt = CCoriolisAij(TwMag,LiMas,Thetap)
%
% Gives the CORIOLIS MATRIX C(t,dt) for the Lagangian's equations: 
% M(t)*ddt + C(t,dt)*dt + N(t,dt) = T
% of the dynamics of the robot formed by links on an open chain. 
%
% INPUTS:
% TwMag = [Tw1; Mag1, ..., Twn; Magn] (7xn)
% for each rigid body joint (1..n).
% Tw1..Twn: The TWIST components for the joint movement.
% Mag1..Magn: The MAGNITUDE component for the joint SCREW movement.
% LiMas = [CM1; IT1; Mas1, ..., CMn; ITn; Masn] (7xn)
% for each rigid body link (1..n).
% CM1..CMn: Center of Mass x-y-z Position components to S for each Link.
% IT1..ITn: Inertia x-y-z components for each Link refered to its CM.
% Mas1..Masn: The Mass for each Link.
% Thetap = Vector differential Theta dt (1xn) of the joint VELOCITIES.
%
%     |C11...C1n|
% C = |         |, With Cij = 1/2 * Sum(1,n)[( dMij/dtk + dMik/dtj - dMkj/dti ) * dtk]  
%     |Cn1...Cnn|
% where:
% dMij/dtk = Sum(l=max(i,j),n)[[Aki*Ei,Ek]'*Alk'*Ml*Alj*Ej + Ei'*Ali'*Ml*Alk*[Akj*Ej,Ek]]
% 
% With Ml being the link inertia nxn LinkInertiaS. 
% With Ei being the twist xi 6x1.
% With Aij being an element 6x6 of the adjoint transformation Aij2adjoint
%
% See also: LinkInertiaS, Aij2adjoint.
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
%% Ctdt = CCoriolisAij(TwMag,LiMas,Thetap)
%
function Ctdt = CCoriolisAij(TwMag,LiMas,Thetap)
%
    n = size(TwMag,2);
    Ctdt = zeros(n);
    for i = 1:n
        for j = 1:n
            for k = 1:n
                Cosij = DerMtetha(i,j,k)+DerMtetha(i,k,j);
                Cosij = Cosij-DerMtetha(k,j,i);
                Cosij = Cosij*Thetap(k);
                Ctdt(i,j)= Ctdt(i,j) + Cosij;
            end
            Ctdt(i,j) = 1/2 * Ctdt(i,j);
        end
    end
    function dMt = DerMtetha(I,J,K)
        m = max(I,J);
        dMt = 0;
        for l = m:n
            Hsli0 = trvP2tform(LiMas(1:3,l));
            Ii = LiMas(4:6,l);
            mi = LiMas(7,l);
            dMt1 = twistbracket(Aij2adjoint(K,I,TwMag)*TwMag(1:6,I),TwMag(1:6,K))';
            dMt1 = dMt1 * Aij2adjoint(l,K,TwMag)'*LinkInertiaS(Hsli0,Ii,mi);
            dMt1 = dMt1 * Aij2adjoint(l,J,TwMag)*TwMag(1:6,J);
            dMt2 = TwMag(1:6,I)'*Aij2adjoint(l,I,TwMag)';
            dMt2 = dMt2 *LinkInertiaS(Hsli0,Ii,mi)*Aij2adjoint(l,K,TwMag);
            dMt2 = dMt2 *twistbracket(Aij2adjoint(K,J,TwMag)*TwMag(1:6,J),TwMag(1:6,K));
            dMt = dMt + dMt1 + dMt2;
        end
    end
end
%

