function [D] = Dissipation(j, Q, S, gam, nodes, dx, bcFlag)
% Calculates the dissipation vector by summing the 2nd and 4th order 
% dissipation terms

% Use boundary conditions
[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR] = BCs(bcFlag);

% Dissipation terms at node j
% [F2,F2L,F4,F4L] = dissipation_terms(S,Q,dx,nodes,bcFlag);
[c4, c2] = dissTerms_node(j, dx, nodes, S, gam, bcFlag, Q);
% % Flow parameters at current value of Q
% [P, rho, u, M, T, c, e] = flowParam(S, Q, nodes);
% 
% % Dissipation terms at node j
% [c4, c2] = dissTerms(j, u, c, P, dx, nodes, S, gam, bcFlag);

% RHS dissipation at node 'j'
Qi=Q(3*j-2:3*j);
if (j>2 && j<nodes-1)
    Qm2=Q(3*(j-2)-2:3*(j-2));
    Qm1=Q(3*(j-1)-2:3*(j-1));
    Qp1=Q(3*(j+1)-2:3*(j+1));
    Qp2=Q(3*(j+2)-2:3*(j+2));
elseif (j==1)
    Qm1=QL;
    Qp1=Q(3*(j+1)-2:3*(j+1));
    Qp2=Q(3*(j+2)-2:3*(j+2));
elseif (j==2)
    Qm2=QL;
    Qm1=Q(3*(j-1)-2:3*(j-1));
    Qp1=Q(3*(j+1)-2:3*(j+1));
    Qp2=Q(3*(j+2)-2:3*(j+2));
elseif (j==nodes-1)
    Qm2=Q(3*(j-2)-2:3*(j-2));
    Qm1=Q(3*(j-1)-2:3*(j-1));
    Qp1=Q(3*(j+1)-2:3*(j+1));
    Qp2=QR;
elseif (j==nodes)
    Qm2=Q(3*(j-2)-2:3*(j-2));
    Qm1=Q(3*(j-1)-2:3*(j-1));
    Qp1=QR;
end

%fourth order term
if (j>1 && j<nodes)
%     c4=[F4(j-1),F4(j)];
%     c2=[F2(j-1),F2(j)];
    fot = c4(2)*(Qm1-3*Qi+3*Qp1-Qp2) - c4(1)*(Qm2-3*Qm1+3*Qi-Qp1);
elseif (j==1)
%     c4=[F4L,F4(j)];
%     c2=[F2L,F2(j)];
    fot = c4(2)*(Qm1-3*Qi+3*Qp1-Qp2) - c4(1)*(-Qm1+2*Qi-Qp1);
elseif (j==nodes)
%     c4=[F4(j-1),F4(j)];
%     c2=[F2(j-1),F2(j)];
    fot = c4(2)*(Qm1-2*Qi+Qp1) - c4(1)*(Qm2-3*Qm1+3*Qi-Qp1);
end

%2nd order term
sot = c2(2)*(Qp1-Qi) - c2(1)*(Qi-Qm1);

D = sot + fot; %Dissipation vector
 
