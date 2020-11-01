% Solve and model a double pendulum systen

clear all; close all; clc;

syms t theta(t) thetadot thetaddot phi(t) phidot phiddot l0 l1 l2 L1 L2 m1 m2  k phi0 theta0 Nx Ny Ox Oy My g
%T1 T2 theta(t) thetadot thetaddot phi(t) phidot phiddot l1 l2 m1 m2 x1ddot y1ddot x2ddot y2ddot Ax Ay Fgrf Fapp g k l0 rx(t) s Ox Oy Mx My

%set up pendulum vector in the global x-y-z fixed co-ordinate system

% Note: the derivatives of phi and theta are NOT set up as functions of time - this
% will come later
r_21x = (l1 - L1) * cos(theta) + L2 * cos(phi);
r_21y = (l1 - L1) * sin(theta) - L2 * sin(phi);

gamma = atan2(r_21y,r_21x);

Fk = k * (norm([r_21x r_21y]) - l0);

r_1GOx = l1 / 2 * cos(theta); % first rod COM relative ot origin
r_1GOy = l1 / 2 * sin(theta);

r_2GNx = l2 / 2 * cos(phi);
r_2GNy = l2 / 2 - sin(phi);

r_2GOx = 2*r_1GOx + r_2GNx;
r_2GOy = 2*r_1GOy + r_2GNy;

Fkx = Fk * cos(gamma);
Fky = Fk * sin(gamma);

I1 = 1/3 * m1 * l1 ^ 2;
I2 = 1/12 * m2 * l2 ^ 2;

eq1 = m1 * diff(r_1GOx,t,t) == Nx + Ox + Fkx;
eq2 = m1 * diff(r_1GOy,t,t) == Ny + Oy + Fky - m1 * g;
eq3 = I1 * diff(theta,t,t) == - m1 * g * cos(theta) * l1/2 + Ny * l1 * cos(theta) - Nx * l1 * sin(theta) - Fkx * L1 * sin(theta) + Fky * L1 * cos(theta);

eq4 = m2 * diff(r_2GOx,t,t) == -Nx - Fkx;
eq5 = m2 * diff(r_2GOy,t,t) == -Ny - m2 * g + My - Fky;
eq6 = I2 * diff(phi,t,t) == My * cos(phi) * l2/2 + Ny * cos(phi) * l2/2 + Nx * sin(phi) * l2/2 - Fkx * L2 * sin(phi) - Fky * cos(phi) * L2;

eq1m = subs(eq1,[diff(phi,t) diff(phi,t,t) diff(theta,t) diff(theta,t,t)],[phidot phiddot thetadot thetaddot]);
eq2m = subs(eq2,[diff(phi,t) diff(phi,t,t) diff(theta,t) diff(theta,t,t)],[phidot phiddot thetadot thetaddot]);
eq3m = subs(eq3,[diff(phi,t) diff(phi,t,t) diff(theta,t) diff(theta,t,t)],[phidot phiddot thetadot thetaddot]);
eq4m = subs(eq4,[diff(phi,t) diff(phi,t,t) diff(theta,t) diff(theta,t,t)],[phidot phiddot thetadot thetaddot]);
eq5m = subs(eq5,[diff(phi,t) diff(phi,t,t) diff(theta,t) diff(theta,t,t)],[phidot phiddot thetadot thetaddot]);
eq6m = subs(eq6,[diff(phi,t) diff(phi,t,t) diff(theta,t) diff(theta,t,t)],[phidot phiddot thetadot thetaddot]);

% acceleration from kinematics in the x-y fixed co-sys
% x1ddot = -l1 * (thetaddot * sin(theta) + thetadot ^ 2 * cos(theta));
% y1ddot = l1 * (thetaddot * cos(theta) - thetadot ^ 2 * sin(theta));
% x2ddot = x1ddot + -l2 * (phiddot * sin(phi) + phidot ^ 2 * cos(phi));
% y2ddot = y1ddot + -l2 * (phiddot * cos(phi) - phidot ^ 2 * sin(phi));
% 
% rx = l1 * cos(theta) + l2 * cos(phi);
% %write F = ma for first rod
% eq1 =  -Ax + Ox == m1 * x1ddot; 
% eq2 =  -Ay + Oy - Fapp  - m1 * g == m1 * y1ddot; 
% eq3 = (Ox+Ax) * sin(theta) - (Oy - Ay -Fapp)*cos(theta) == m1 * l1 * thetaddot / 6;
% 
% % write F = ma for the second rod (name them eq3 and eq4)
% eq4 = Ax - k * (rx - l0) == m2 * x2ddot;
% eq5 = Ay + My - m2 * g == m2 * y2ddot;
% eq6 = (My - Ay) * cos(phi) - (Ax + k* (rx - l0)) * sin(phi) == m2 * l2 * phiddot / 6;
% 
% eq7 = Oy + My == Fapp;
% %This step isolates phiddot, thetaddot, T1, and T2. 

% Normally this would be done on paper with lots of algebra (e.g. solving one equation for T1 and plugging it
%into another equation to eliminate, and then repeating the process for
%phiddot and thetaddot )
system = solve(eq1m,eq2m,eq3m,eq4m,eq5m,eq6m,phiddot,thetaddot,Nx,Ny,My,Oy);
%%
%Now we need to replace phidot and thetadot with the time derivatives so that
%we can solve the differential equations
eq_thetaddot = subs(system.thetaddot,[phidot thetadot],[diff(phi,t) diff(theta,t)]);
eq_phiddot = subs(system.phiddot,[phidot thetadot],[diff(phi,t) diff(theta,t)]);

%set the appropriate equations equal to the second derivative of theta and
%the second derivative of x
eq_thetaddot = diff(theta,t,t) == eq_thetaddot;
eq_phiddot = diff(phi,t,t) == eq_phiddot;

%This function converts the system of two second order differential
%equations into a series of first order equations
[eqss varss newvarrs] = reduceDifferentialOrder([eq_thetaddot; eq_phiddot],[theta phi])

%This block of code isolates the extra parameters (e.g. mass of the 
%pendulum, pendulum length, gravity) We knew what
%they were anyway but this ensures nothing is missed. 
pDAEs = symvar(eqss);
pDAEvars = symvar(varss)
extraParams = setdiff(pDAEs,pDAEvars)

%convert symbolic experessions to numerical
f = daeFunction(eqss,varss,L1, L2, Ox, g, k, l0, l1, l2, m1, m2)

%%
clearvars -except f
g = 9.81; l1 = 0.1;  m1 = 0.02; m2 = 0.04;
L1= 0.1*l1; 
Ox = 0; 
phi0 = 35 * pi / 180; theta0 = 25*pi/180;
% phi_l0 = 13 * pi / 180; theta_l0 = 14*pi/180;
Fapp = 100; k = 100; 
h = sin(theta0) * l1;
l2 = h/sin(phi0);
L2 = l2 * 0.8; 

% make functions of the location of each mass
rx1 = @(theta) l1*cos(theta);
ry1 = @(theta) l1*sin(theta);
rx2 = @(theta,phi) rx1(theta) + l2 * cos(phi);
ry2 = @(theta,phi) ry1(theta) - l2 * sin(phi);

% length of the spring, determine resting length
l0x = @(theta,phi) (l1 - L1) * cos(theta) + L2 * cos(phi);
l0y = @(theta,phi) (l1 - L1) * sin(theta) - L2 * sin(phi);
l0 = norm([l0x(theta0,phi0) l0y(theta0,phi0)]);

% l0  = rx2(theta0,phi0);

F = @(t,Y,Yp) f(t, Y, Yp,L1, L2, Ox, g, k, l0, l1, l2, m1, m2);
opt = odeset('RelTol', 10.0^(-8), 'AbsTol' , 10.0^(-8));

y0 = [theta0, phi0 , 0, 0];
yp0 = zeros(1,4);

%use decic to estimate initial guesses
[y0 yp0] = decic(F,0,y0,[],yp0,[],opt)
dt = 0.001 %specify timestep (runs in adaptive mode (variable timestep) if not specified)
[tSol ySol] = ode15i(F,[0:dt:1],y0, yp0,opt)
figure;
subplot(2,1,1)
plot(tSol,ySol(:,1:2),'-o')
legend('theta','phi')
subplot(2,1,2)
plot(tSol,ySol(:,3:4),'-o')
legend('vel_theta','vel_phi')
%visualize 
   figure;
   plot(0,0,'ko')
   ind = 1;
   for iii = 1
for i = 1:length(tSol) %skipping every other frame
   cla reset;
   hold on
   rx1i = rx1(ySol(i,1));
   ry1i = ry1(ySol(i,1));
   rx2i = rx2(ySol(i,1),ySol(i,2));
   ry2i = ry2(ySol(i,1),ySol(i,2));
   
   plot([0,rx1i],[0, ry1i]) % plot the first rope
   plot([rx1i,rx2i], [ry1i, ry2i]) % plot the second rope
%    plot(rx1(ySol(i,1)),ry1(ySol(i,1)),'.','MarkerSize',5*m1) % plot the mass proportional to the magnitude for the first mass
%    plot(rx2(ySol(i,1),ySol(i,2)),ry2(ySol(i,1),ySol(i,2)),'.','MarkerSize',5*m2) % plot the second mass

   %double check the function by testing that the lengths of the strings are
%consistent
   L1s(ind,1) = norm([rx1i,ry1i]);
   L2s(ind,1) = norm([rx1i,ry1i]-[rx2i,ry2i]);
 
   % sin waves to represent 
   ncoils = 10;
   lrxi = L1 * cos(ySol(i,1));
   lryi = L1 * sin(ySol(i,1));
   lxi = l0x(ySol(i,1),ySol(i,2));
   lyi = l0y(ySol(i,1),ySol(i,2));
   
   xx1e = linspace(lrxi,lrxi+lxi,100);
   yy1 = linspace(lryi,lryi+lyi,100);
   yy1e = 0.005*sin((2*pi/((rx1i+lxi)/ncoils))*xx1e)+yy1;
   plot(xx1e, yy1e, 'r')
   % set limits on the plots 
   xlim([-0.25 0.5])
   ylim([-0.1 0.2])
   drawnow
   pause(0.05)
   hold off
   
   ind = ind+1;
end

   end




