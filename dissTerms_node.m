function [c4, c2] = dissTerms_node(j, dx, nodes, S, gam, bcFlag, Q)

% Include boundary conditions
[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR] = BCs(bcFlag);
    
% Define constants
k2 = 0.5; % Tuning parameter for 2nd order dissipation terms
k4 = 1/32; % Tuning parameter for 4th order dissipation terms

%Flow parameters at node j
[Pi, rhoi, ui, Mi, Ti, ci, ei] = flowParam_node(S, Q, j);
if (j>3 && j<nodes-2)
    [Pm3, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j-3);
    [Pm2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j-2);
    [Pm1, ~, um, ~, ~, cm, ~] = flowParam_node(S, Q, j-1);
    [Pp1, ~, up, ~, ~, cp, ~] = flowParam_node(S, Q, j+1);
    [Pp2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j+2);
    [Pp3, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j+3);
elseif(j==1)
    Pm1=PL;
    um=uL;
    cm=cL;
    [Pp1, ~, up, ~, ~, cp, ~] = flowParam_node(S, Q, j+1);
    [Pp2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j+2);
    [Pp3, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j+3);
elseif(j==2)
    Pm2=PL;
    [Pm1, ~, um, ~, ~, cm, ~] = flowParam_node(S, Q, j-1);
    [Pp1, ~, up, ~, ~, cp, ~] = flowParam_node(S, Q, j+1);
    [Pp2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j+2);
    [Pp3, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j+3);
elseif(j==3)
    Pm3=PL;
    [Pm2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j-2);
    [Pm1, ~, um, ~, ~, cm, ~] = flowParam_node(S, Q, j-1);
    [Pp1, ~, up, ~, ~, cp, ~] = flowParam_node(S, Q, j+1);
    [Pp2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j+2);
    [Pp3, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j+3);
elseif (j==nodes-2)
    [Pm3, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j-3);
    [Pm2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j-2);
    [Pm1, ~, um, ~, ~, cm, ~] = flowParam_node(S, Q, j-1);
    [Pp1, ~, up, ~, ~, cp, ~] = flowParam_node(S, Q, j+1);
    [Pp2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j+2);
    Pp3=PR;
elseif (j==nodes-1)
    [Pm3, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j-3);
    [Pm2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j-2);
    [Pm1, ~, um, ~, ~, cm, ~] = flowParam_node(S, Q, j-1);
    [Pp1, ~, up, ~, ~, cp, ~] = flowParam_node(S, Q, j+1);
    Pp2=PR;
elseif (j==nodes)
    [Pm3, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j-3);
    [Pm2, ~, ~, ~, ~, ~, ~] = flowParam_node(S, Q, j-2);
    [Pm1, ~, um, ~, ~, cm, ~] = flowParam_node(S, Q, j-1);
    Pp1=PR;
    up=uR;
    cp=cR;
end
     
% Upsilon
ups_i = abs(Pp1 - 2*Pi + Pm1)/abs(Pp1 + 2*Pi + Pm1);
if (j>2 && j<nodes-1)
    ups_m2 = abs(Pm1 - 2*Pm2 + Pm3)/abs(Pm1 + 2*Pm2 + Pm3);
    ups_m1 = abs(Pi - 2*Pm1 + Pm2)/abs(Pi + 2*Pm1 + Pm2);
    ups_p1 = abs(Pp2 - 2*Pp1 + Pi)/abs(Pp2 + 2*Pp1 + Pi);
    ups_p2 = abs(Pp3 - 2*Pp2 + Pp1)/abs(Pp3 + 2*Pp2 + Pp1);
elseif (j==1)
    ups_p1 = abs(Pp2 - 2*Pp1 + Pi)/abs(Pp2 + 2*Pp1 + Pi);
    ups_p2 = abs(Pp3 - 2*Pp2 + Pp1)/abs(Pp3 + 2*Pp2 + Pp1);
elseif (j==2)
    ups_m1 = abs(Pi - 2*Pm1 + Pm2)/abs(Pi + 2*Pm1 + Pm2);
    ups_p1 = abs(Pp2 - 2*Pp1 + Pi)/abs(Pp2 + 2*Pp1 + Pi);
    ups_p2 = abs(Pp3 - 2*Pp2 + Pp1)/abs(Pp3 + 2*Pp2 + Pp1);
elseif (j==nodes-1) 
    ups_m2 = abs(Pm1 - 2*Pm2 + Pm3)/abs(Pm1 + 2*Pm2 + Pm3);
    ups_m1 = abs(Pi - 2*Pm1 + Pm2)/abs(Pi + 2*Pm1 + Pm2);
    ups_p1 = abs(Pp2 - 2*Pp1 + Pi)/abs(Pp2 + 2*Pp1 + Pi);
elseif (j==nodes)
    ups_m2 = abs(Pm1 - 2*Pm2 + Pm3)/abs(Pm1 + 2*Pm2 + Pm3);
    ups_m1 = abs(Pi - 2*Pm1 + Pm2)/abs(Pi + 2*Pm1 + Pm2);
end
     
% 2nd order dissipation 'e2'
if (j>2 && j<nodes-1)
    e2_m = k2*max(max(ups_i, ups_m1), ups_m2);
    e2_i = k2*max(max(ups_p1, ups_i), ups_m1);
    e2_p = k2*max(max(ups_p2, ups_p1), ups_i);
elseif (j==1)
    e2_i = k2*max(ups_p1, ups_i); 
    e2L = e2_i;
    e2_p = k2*max(max(ups_p2, ups_p1), ups_i);
elseif (j==2)
    e2_m = k2*max(ups_p1, ups_i); 
    e2_i = k2*max(max(ups_p1, ups_i), ups_m1);
    e2_p = k2*max(max(ups_p2, ups_p1), ups_i);
elseif (j==nodes-1)
    e2_m = k2*max(max(ups_i, ups_m1), ups_m2);
    e2_i = k2*max(max(ups_p1, ups_i), ups_m1);
    e2_p = k2*max(ups_i, ups_m1); 
elseif (j==nodes)
    e2_m = k2*max(max(ups_i, ups_m1), ups_m2);
    e2_i = k2*max(ups_i, ups_m1); 
    e2R = e2_i;
end

% 4th order dissipation 'e4'
% Spectral radius 'sigma'
e4_i = max(0, (k4 - e2_i));
sigma_i = (abs(ui) + ci)/dx;
if (j>1 && j<nodes)
    e4_m = max(0, (k4 - e2_m));
    e4_p = max(0, (k4 - e2_p));
    sigma_m = (abs(um) + cm)/dx;
    sigma_p = (abs(up) + cp)/dx;
elseif (j==1)
    e4L = max(0, (k4 - e2L));
    e4_p = max(0, (k4 - e2_p));
    sigma_m = (abs(uL) + cL)/dx;
    sigma_p = (abs(up) + cp)/dx;
elseif (j==nodes)
    e4_m = max(0, (k4 - e2_m));
    e4R = max(0, (k4 - e2R));
    sigma_m = (abs(um) + cm)/dx;
    sigma_p = (abs(uR) + cR)/dx;
end

% Calculate values of dissipation terms at j
if (j>1 && j<nodes)
    c4(2) = (sigma_i*e4_i + sigma_p*e4_p)/2;
    c4(1) = (sigma_i*e4_i + sigma_m*e4_m)/2;   
    c2(2) = (sigma_i*e2_i + sigma_p*e2_p)/2;
    c2(1) = (sigma_i*e2_i + sigma_m*e2_m)/2;
    
elseif (j==1)
    c4(2) = (sigma_i*e4_i + sigma_p*e4_p)/2;
    c4(1) = (sigma_i*e4_i + sigma_m*e4L)/2;
    c2(2) = (sigma_i*e2_i + sigma_p*e2_p)/2;
    c2(1) = (sigma_i*e2_i + sigma_m*e2L)/2;
    
elseif (j==nodes)
    c4(2) = (sigma_i*e4_i + sigma_p*e4R)/2;
    c4(1) = (sigma_i*e4_i + sigma_m*e4_m)/2;  
    c2(2) = (sigma_i*e2_i + sigma_p*e2R)/2;
    c2(1) = (sigma_i*e2_i + sigma_m*e2_m)/2;
end

    
end

    












