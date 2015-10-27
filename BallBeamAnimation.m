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


% Record keeping code for how long it takes to simulate (deprecated)
clear global
clearvars

syms theta thetadot y ydot u s;

% Constants
I2 = 0.05;
I1 = 0.2;
m2 = 2;
r = 0.11;
a = 0.1;
g = 9.81;

% Initial Conditions
x0 = [-22 * pi / 180; 0.9; 0; 0];
x0(1) = input('Enter starting angle (deg): ') * pi / 180;
x0(2) = input('Enter starting position (m): ');

% W matrix, given in handout
W = [I1 + I2 + m2 * (y^2 + r^2) m2 * r + I2 / a; m2 * r + I2 / a m2 + I2 / a^2];

% Q double dot and Q dot matrices (theta double dot, y double dot) (theta
% dot, ydot)
qdd = W^(-1) * ([u; 0] - (m2 * [2*y*thetadot*ydot; -y*thetadot^2] - m2*g*[r*sin(theta) + y*cos(theta); sin(theta)]));
qd = [thetadot; ydot];

% combining qdot, qdouble dot
xdot = [qd; qdd];

% equilibrium conditions for jacobian matrix solution
eq = [0 0 0 0 0 ];

% partial derivatives of theta dot
% diff([equation to differentiate], [variable to differentiate with respect
% to])
% subs([equation], [symbolic variables], [values for the symbolic
% variables])
x3bytheta = subs(diff(qdd(1), theta), {theta, y, thetadot, ydot, u}, eq);
x3byy = subs(diff(qdd(1), y), {theta, y, thetadot, ydot, u}, eq);
x3bythetadot = subs(diff(qdd(1), thetadot), {theta, y, thetadot, ydot, u}, eq);
x3byydot = subs(diff(qdd(1), ydot), {theta, y, thetadot, ydot, u}, eq);

% partial derivatives of y dot
x4bytheta = subs(diff(qdd(2), theta), {theta, y, thetadot, ydot, u}, eq);
x4byy = subs(diff(qdd(2), y), {theta, y, thetadot, ydot, u}, eq);
x4bythetadot = subs(diff(qdd(2), thetadot), {theta, y, thetadot, ydot, u}, eq);
x4byydot = subs(diff(qdd(2), ydot), {theta, y, thetadot, ydot, u}, eq);

% partial derivatives for linearized u matrix
x3byu = subs(diff(qdd(1), u), {theta, y, thetadot, ydot, u}, eq);
x4byu = subs(diff(qdd(2), u), {theta, y, thetadot, ydot, u}, eq);

% A row one: theta dot
% A row two: y dot
% A row three: linearized from above
% A row four: linearized from above
A = [0 0 1 0; 0 0 0 1; x3bytheta x3byy x3bythetadot x3byydot; x4bytheta x4byy x4bythetadot x4byydot];
B = [0; 0; x3byu; x4byu];

% The following does not execute on some computers/versions of matlab
% without the 'double' typecast for matrices A and B, the try/catch block
% should therefore be deprecated, however, in case of any further symengine
% errors there is always a hard-coded backup version of K
try
    K = place(double(A), double(B), [-1, -2, -1.5+2i, -1.5-2i]);
catch
    %K value calculated on a computer which will successfully execute the
    %'place' function without the 'double' typecast
    K = [ 3.6844157929226736566186107470511, 20.512584097859327217125382262997, 1.3826382699868938401048492791612, 1.7673165137614678899082568807339]
end

% allocate space for an output matrix, four rows, 1600 columns
% output = zeros(4, 1600);

uVal = 0;

xdot1 = matlabFunction(simplify(xdot(1)));
xdot2 = matlabFunction(simplify(xdot(2)));
xdot3 = matlabFunction(xdot(3));
xdot4 = matlabFunction(xdot(4));

for i = 1:1600
    
    % if this is the start, use values defined outside of for loop
    % else, redefine x0 and uVal to most recent values
    if i > 1
       x = output(1:4, i - 1);
       u = -K * x;
    else
        x = x0;
        u = 0;
    end
    
    theta = x(1);
    y = x(2);
    thetadot = x(3);
    ydot = x(4);
   
    % Runge Kutte algorithm, see handout
    z1 = double([xdot1(thetadot);...
                 xdot2(ydot);...
                 xdot3(thetadot, theta, u, y, ydot);...
                 xdot4(thetadot, theta, u, y, ydot)] * 0.01);
    
    z2 = double([xdot1(thetadot + z1(3) / 2);...
                 xdot2(ydot + z1(4) / 2);...
                 xdot3(thetadot + z1(3) / 2, theta + z1(1) / 2, u, y + z1(2) / 2, ydot + z1(4) / 2);...
                 xdot4(thetadot + z1(3) / 2, theta + z1(1) / 2, u, y + z1(2) / 2, ydot + z1(4) / 2)] * 0.01);
    
    z3 = double([xdot1(thetadot + z2(3) / 2);...
                 xdot2(ydot + z2(4) / 2);...
                 xdot3(thetadot + z2(3) / 2, theta + z2(1) / 2, u, y + z2(2) / 2, ydot + z2(4) / 2);...
                 xdot4(thetadot + z2(3) / 2, theta + z2(1) / 2, u, y + z2(2) / 2, ydot + z2(4) / 2)] * 0.01);         
    
    z4 = double([xdot1(thetadot + z3(3) / 2);...
                 xdot2(ydot + z3(4) / 2);...
                 xdot3(thetadot + z3(3) / 2, theta + z3(1) / 2, u, y + z3(2) / 2, ydot + z3(4) / 2);...
                 xdot4(thetadot + z3(3) / 2, theta + z3(1) / 2, u, y + z3(2) / 2, ydot + z3(4) / 2)] * 0.01);
    
    % Determine output value and assign to current position in output
    % matrix
    t(:,i) = 0.01 * i;
    output(:,i) = (z1 + 2 * z3 + 2 * z3 + z4) / 6 + x;

end

plot(t(1:1600), output(1,1:1600),t(1:1600), output(2,1:1600));
legend('Angle (rad)','Distance (m)');
legend('show');
xlim([0 16]);
    
drawnow;

%% 

fig=figure('DeleteFcn',@closefigfcn);
axs=axes('Parent',fig);

fill([-0.5, 0.5, 0], [-1, -1, 0],'b','Parent',axs);

beamLength = 2;
ball = rectangle('Position',[(x0(2) - a / 2) (-x0(2) * sin(x0(1))) (2 * a) (2 * a)],...
                 'FaceColor','r',...
                 'Curvature',[1, 1],...
                 'Parent',axs...
                 );
beam = line([beamLength / 2 * cos(x0(1)), -beamLength / 2 * cos(x0(1))], [-beamLength / 2 * sin(x0(1)), beamLength / 2 * sin(x0(1))],...
   'LineWidth',2,...
   'Parent',axs...
   );
axis(axs, [-1.5, 1.5, -0.5, 0.5]);
axis(axs, 'equal');

drawnow;


for i = 1:400
    
   startTime = tic();
   index = i * 4;
   
   ball.Position = [output(2, index) * cos(-output(1, index)) - a - a * sin(-output(1, index)), output(2, index) * sin(-output(1, index)) + a * cos(-output(1, index)) - a, 2 * a, 2 * a]; 
   
   beam.XData = [-beamLength / 2 * cos(output(1, index)), beamLength / 2 * cos(output(1, index))];
   beam.YData = [beamLength / 2 * sin(output(1, index)), -beamLength / 2 * sin(output(1, index))];

   axis(axs, [-1.5, 1.5, -0.5, 0.5]);
   
   if i > 1
       delete(h);
   end
   h = text(-1.25, -0.35, {strcat('Angle (deg): ', num2str(round(output(1, index) * 180 / pi, 3))), strcat('Distance (m): ', num2str(round(output(2, index), 3))), strcat('Time (s): ', num2str(t(index)))}); 
   drawnow;
   pause(0.04 - toc(startTime))
end
%%
% Edit axes of plot
% xlabel('Time (s * 100)')
% ylabel('Y Position (m)')
