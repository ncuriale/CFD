function [P, rho, u, M, T, c, e] = flowParam(S, Q, nodes)
%   This subroutine will calculate the updated values of the flow
%   parameters at each node of the grid.  

% Define constants
R = 287; %specific gas constant
gam = 1.4; %specific heat ratio

for j=1:nodes
    P(j) = ((gam-1)/S(j))*(Q(3*j)-0.5*(Q(3*j-1)^2)/Q(3*j-2));%pressure
    rho(j) = Q(3*j-2)/S(j);%density
    u(j) = Q(3*j-1)/Q(3*j-2);%velocity
    c(j) = sqrt(gam*P(j)*S(j)/Q(3*j-2));% sound speed
    M(j) = u(j)/c(j);%Mach
    T(j) = P(j)/(R*rho(j));%temperature
    e(j) = Q(3*j)/S(j);%internal energy 'e'
end