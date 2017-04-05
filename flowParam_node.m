function [P, rho, u, M, T, c, e] = flowParam_node(S, Q, j)
%   This subroutine will calculate the updated values of the flow
%   parameters at each node of the grid.  

% Define constants
R = 287; %specific gas constant
gam = 1.4; %specific heat ratio

P = ((gam-1)/S(j))*(Q(3*j)-0.5*(Q(3*j-1)^2)/Q(3*j-2));%pressure
rho = Q(3*j-2)/S(j);%density
u = Q(3*j-1)/Q(3*j-2);%velocity
c = sqrt(gam*P*S(j)/Q(3*j-2));% sound speed
M = u/c;%Mach
T = P/(R*rho);%temperature
e = Q(3*j)/S(j);%internal energy 'e'