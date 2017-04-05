function [F2,F2L,F4,F4L] = dissipation_terms(S,Q,dx,N,bcFlag)

gamma =1.4;
k2 = 0.5;
k4 = 1/32;

% Include boundary conditions
[uL, PL, rhoL, eL, SL, QL, cL, uR, PR, rhoR, eR, SR, QR, cR] = BCs(bcFlag);
% Flow parameters at current value of Q
[P, rho, u, M, T, c, e] = flowParam(S, Q, N);

% Compute wave_speeds
lambda = abs(u) + c;
lambdaBC = abs([uL,uR]) + [cL,cR];

% Compute Gammas where Gamma_j = [eps(i)*lambda]_j
p_check = zeros(N,1);
p_check(1) = abs((P(2)-2*P(1)+PL)/(P(2)+2*P(1)+PL));
for j = 2:N-1
    p_check(j) = abs((P(j+1)-2*P(j)+P(j-1))/(P(j+1)+2*P(j)+P(j-1)));
end
p_check(N) = abs((PR-2*P(N)+P(N-1))/(PR+2*P(N)+P(N-1)));

eps2 = zeros(N,1);eps4=eps2;G2=eps2;G4=eps2; % Allocation
for j=1:N
    if j ==1
        eps2(j) = k2*max([p_check(j+1),p_check(j)]);
    elseif j < N
        eps2(j) = k2*max([p_check(j+1),p_check(j),p_check(j-1)]);
    else
        eps2(j) = k2*max([p_check(j),p_check(j-1)]);            
    end
    eps4(j) = max(0,k4-eps2(j));
    G2(j)   = eps2(j)*lambda(j);
    G4(j)   = eps4(j)*lambda(j);
end

% Get terms for BC
eps2BC(1) = k2*max(p_check(1),p_check(2));
eps2BC(2) = k2*max(p_check(N-1),p_check(N));
eps4BC(1) = max(0,k4-eps2BC(1));
eps4BC(2) = max(0,k4-eps2BC(2));

G2BC(1)   = eps2BC(1)*lambdaBC(1);
G2BC(2)   = eps2BC(2)*lambdaBC(2);
G4BC(1)   = eps4BC(1)*lambdaBC(1);
G4BC(2)   = eps4BC(2)*lambdaBC(2);

% Get j+1/2 Gamma values

F2L = 1/2/dx*(G2BC(1)+G2(1));
F4L = 1/2/dx*(G4BC(1)+G4(1));
F2 = zeros(N,1);
F4 = zeros(N,1);

for i = 1:N
  if i<N
    F2(i) = 1/2/dx*(G2(i)+G2(i+1));
    F4(i) = 1/2/dx*(G4(i)+G4(i+1));
  else
    F2(i) = 1/2/dx*(G2(i)+G2BC(2));
    F4(i) = 1/2/dx*(G4(i)+G4BC(2));
  end
end