function [FD] = FluxDiff(j, Q, S, nodes, gam, dx, bcFlag)
% Calculates the value of the 2nd order centered flux difference at node 'j'.

% Include boundary conditions
[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR] = BCs(bcFlag);

if (j==1)
    FD(1,1) = Q(3*(j+1)-1) - QL(2);
    FD(2,1) = [(gam-1)*Q(3*(j+1))-((gam-3)/2)*(Q(3*(j+1)-1)^2)/Q(3*(j+1)-2)] - [(gam-1)*QL(3)-((gam-3)/2)*(QL(2)^2)/QL(1)];
    FD(3,1) = [gam*Q(3*(j+1))*Q(3*(j+1)-1)/Q(3*(j+1)-2)-((gam-1)/2)*(Q(3*(j+1)-1)^3)/(Q(3*(j+1)-2)^2)]  - [gam*QL(3)*QL(2)/QL(1)-((gam-1)/2)*(QL(2)^3)/(QL(1)^2)];
    FD = FD/(2*dx); 
    
elseif (j<nodes)
    FD(1,1) = Q(3*(j+1)-1) - Q(3*(j-1)-1);
    FD(2,1) = [(gam-1)*Q(3*(j+1))-((gam-3)/2)*(Q(3*(j+1)-1)^2)/Q(3*(j+1)-2)] - [(gam-1)*Q(3*(j-1))-((gam-3)/2)*(Q(3*(j-1)-1)^2)/Q(3*(j-1)-2)];
    FD(3,1) = [gam*Q(3*(j+1))*Q(3*(j+1)-1)/Q(3*(j+1)-2)-((gam-1)/2)*(Q(3*(j+1)-1)^3)/(Q(3*(j+1)-2)^2)] - [gam*Q(3*(j-1))*Q(3*(j-1)-1)/Q(3*(j-1)-2)-((gam-1)/2)*(Q(3*(j-1)-1)^3)/(Q(3*(j-1)-2)^2)];
    FD = FD/(2*dx); 
 
else
    FD(1,1) = QR(2) - Q(3*(j-1)-1);
    FD(2,1) = [(gam-1)*QR(3)-((gam-3)/2)*(QR(2)^2)/QR(1)] - [(gam-1)*Q(3*(j-1))-((gam-3)/2)*(Q(3*(j-1)-1)^2)/Q(3*(j-1)-2)];
    FD(3,1) = [gam*QR(3)*QR(2)/QR(1)-((gam-1)/2)*(QR(2)^3)/(QR(1)^2)] - [gam*Q(3*(j-1))*Q(3*(j-1)-1)/Q(3*(j-1)-2)-((gam-1)/2)*(Q(3*(j-1)-1)^3)/(Q(3*(j-1)-2)^2)];
    FD = FD/(2*dx); 
    
end