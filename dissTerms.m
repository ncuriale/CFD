function [c4, c2] = dissTerms(j, u, c, P, dx, nodes, S, gam, bcFlag)
% Calculates the values of the 4th order and 2nd order dissipation terms.

% Include boundary conditions
[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR] = BCs(bcFlag);

% Define constants
k2 = 0.5; % Tuning parameter for 2nd order dissipation terms
k4 = 1/50; % Tuning parameter for 4th order dissipation terms

% Upsilon
ups(1) = abs(P(2) - 2*P(1) + PL)/abs(P(2) + 2*P(1) + PL); % inflow boundary
for k=2:(nodes-1);
    ups(k) = abs(P(k+1) - 2*P(k) + P(k-1))/abs(P(k+1) + 2*P(k) + P(k-1));
end
ups(nodes) = abs(PR - 2*P(nodes) + P(nodes-1))/abs(PR + 2*P(nodes) + P(nodes-1)); % outflow boundary

% 2nd order dissipation 'e2'
e2(1) = k2*max(ups(2), ups(1)); % inflow boundary
e2L = e2(1);
for k=2:(nodes-1)
    e2(k) = k2*max(max(ups(k+1), ups(k)), ups(k-1)); 
end
e2(nodes) = k2*max(ups(nodes), ups(nodes-1)); % outflow boundary
e2R = e2(nodes);

% 4th order dissipation 'e4'
e4L = max(0, (k4 - e2L));
for k=1:nodes;
    e4(k) = max(0, (k4 - e2(k)));
end
e4R = max(0, (k4 - e2R));

% Spectral radius 'sigma'
sigmaL = (abs(uL) + cL)/dx;
for k=1:nodes;
    sigma(k) = (abs(u(k)) + c(k))/dx;
end
sigmaR = (abs(uR) + cR)/dx;

% Calculate values of dissipation terms at j
if (j==1)
    c4(2) = (sigma(j)*e4(j) + sigma(j+1)*e4(j+1))/2;
    c4(1) = (sigma(j)*e4(j) + sigmaL*e4L)/2;
    c2(2) = (sigma(j)*e2(j) + sigma(j+1)*e2(j+1))/2;
    c2(1) = (sigma(j)*e2(j) + sigmaL*e2L)/2;
    
elseif (j==nodes)
    c4(2) = (sigma(j)*e4(j) + sigmaR*e4R)/2;
    c4(1) = (sigma(j)*e4(j) + sigma(j-1)*e4(j-1))/2;  
    c2(2) = (sigma(j)*e2(j) + sigmaR*e2R)/2;
    c2(1) = (sigma(j)*e2(j) + sigma(j-1)*e2(j-1))/2;
    
else
    c4(2) = (sigma(j)*e4(j) + sigma(j+1)*e4(j+1))/2;
    c4(1) = (sigma(j)*e4(j) + sigma(j-1)*e4(j-1))/2;   
    c2(2) = (sigma(j)*e2(j) + sigma(j+1)*e2(j+1))/2;
    c2(1) = (sigma(j)*e2(j) + sigma(j-1)*e2(j-1))/2;  
    
end


% tol=1e-5;
% if (c4(1)<tol) 
%     c4(1)=0.0;
% elseif(c4(2)<tol)
%     c4(2)=0.0;
% elseif(c2(1)<tol)
%     c2(1)=0.0;
% elseif(c2(2)<tol)
%     c2(2)=0.0;
% end
% 

    












