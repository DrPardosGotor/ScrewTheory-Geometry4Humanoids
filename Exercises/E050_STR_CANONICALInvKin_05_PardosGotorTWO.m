%% Screw Theory - CANONICAL Subproblem - Exercise Pardos-Gotor TWO (PG2).
%
% Calculate IK for two consecutive SCREWS by PardosGotorTwo function.
% the movements are defined by the SCREWS whose "Twists" parameters
% are: Axis = [Axis1 Axis2], Point = [p1 p2], JointType = 'tra'
% and whose magnitudes are defined by Mag = [Theta1 Theta2].
%
% The magnitude Theta2 is aplied to point pp for moving it to pc
% then the magnitude Theta1 is aplied to pc for moving them to pk.
%
% For checking the PardosGotorTwo function, this exercise has 3 steps:
% STEP1: Apply ForwardKinemats to the Screws for "random" Mag (t2 + t1)
% on pp and then getting a feasible pk.
% STEP2: Calculate the IK solution by PKP2 getting the magnitud
% Theta1Theta2 = [t11 t21] SOLUTION.
% STEP3: Test the solutions got by PKP2 Theta1 & Theta2 applying
% ForwardKinemats to the Screws on pp and checking we get the same pk.
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
%% E050_STR_CANONICALInvKin_05_PardosGotorTWO
%
%
clear
clc
%
pp = [rand*10 rand*10 rand*10]' % for testing various initial points
Mag = [rand rand]; % for testing various magnitudes
%
p1 = [0 0 0]'; p2 = [0 0 0]'; % must imply the intersection of Twists
Point = [p1 p2];
AxisX = [1 0 0]'; AxisY = [0 1 0]'; AxisZ = [0 0 1]';
Axis = [AxisX AxisY]; % whatever for testing the exercise
JointType = ['tra'; 'tra']; % whatever for testing the exercise
% 
% Now we build the TWISTS matrix for the chosen Joints
Twist = joint2twist(Axis(:,1), Point(:,1), JointType(1,:));
for i = 2:size(Point,2)
    Twist = [Twist joint2twist(Axis(:,i), Point(:,i), JointType(i,:))];
end
%
% STEP1: Apply ForwardKinemats to the TWO Screws x2 and then x1 on pp for
% "whatever" Mag (can be% even more than 2pi) for getting a feasible pk.
TwMag1 = [Twist; Mag];
HstR1 = ForwardKinematicsPOE(TwMag1);
pk1h = HstR1*[pp; 1];
pk1 = pk1h(1:3)
%
% STEP2: Calculate the IK solution by PK2 getting the magnitud
% Theta1Theta2 = [t11 t21; t12 t22] DOUBLE SOLUTION.
Th1Th2 = PardosGotorTwo(Twist(:,1), Twist(:,2), pp, pk1)
%
% STEP3: Test the TWO solutions by PG2 Theta1 & Theta2 applying
% ForwardKinemats to the Screws on pp and checking we get the same pk.
TwMag2 = [Twist; Th1Th2];
HstR2 = ForwardKinematicsPOE(TwMag2);
pk2h = HstR2*[pp; 1];
pk2 = pk2h(1:3)
%
% Check that (pk1 = pk2 = pk3) 
%