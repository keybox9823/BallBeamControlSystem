% Copyright 2015 Michael Estes

% This file is part of BallBeam

% BallBeam is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.t
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [ z ] = calcZ( x0, uVal )
% calcZ([4x1 Array], [scalar input value])
% This function accepts two arguments, the first is a four by one array
% representing the current X matrix value, the second is a scalar value of
% u, or the input to the system. The function returns a 4x1 matrix. The
% equations used here were determined from the symbolic xdot matrix defined
% in BallBeamSimulation.m, however, those equations have been simplified here and
% presented using variable precision arithmetic, removing the symbolic
% components, in order to enhance speed of simulation.

z = [x0(3);...
     x0(4);...
     ((- 1440.0*x0(2)*x0(3)^2 - 28000.0*x0(2)*x0(4)*x0(3) + 7000.0*uVal + 981.0*sin(x0(1)) + 137340.0*x0(2)*cos(x0(1)))/(14000.0*x0(2)^2 + 1401.0));...
     ((0.2*(10000.0*x0(2)^2 + 1371.0)*(2.0*x0(2)*x0(3)^2 + 19.62*sin(x0(1))))/(14000.0*x0(2)^2 + 1401.0) - (720.0*(uVal + 2.1582*sin(x0(1)) + 19.62*x0(2)*cos(x0(1)) - 4.0*x0(3)*x0(2)*x0(4)))/(14000.0*x0(2)^2 + 1401.0))...
    ] * 0.01;

end

